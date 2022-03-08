#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Union
from io import BytesIO

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.morphology import white_tophat, disk
from scipy.ndimage import maximum_filter
from imctools.io.mcd.mcdparser import McdParser
from imctools.io.mcd.mcdxmlparser import McdXmlParser

from stardist.models import StarDist2D
from csbdeep.utils import Path, normalize

from skimage.measure import find_contours, grid_points_in_poly
from skimage.filters import gaussian
from skimage import morphology
from skimage import exposure


_curr_loc = Path(__file__).parent
SPILLMAT_CSV = (_curr_loc / ".." / ".." / "databases" / "spillover.csv").resolve()
#SPILLMAT_CSV = Path("/projects/robson-lab/research/imc/data/spillover.csv")


def cross_dist(u, v, w=None):
    if w is None:
        w = np.zeros_like(u)
    return np.sign(np.cross(u - w, v - w))


def close_contour(contour, xmax, ymax):
    if len(contour) < 5:
        return contour
    
    c1, c0 = contour[[-1, 0]]

    if np.equal(c1[0], c0[0]):
        return contour
    elif np.equal(c1[1], c0[1]):
        return contour
    
    corners = np.array([[0,0], [0,ymax], [xmax,ymax], [xmax,0]])
    crosses = cross_dist(c1, c0, corners)
    
    return np.vstack((
        contour,
        corners[crosses < 0]
    ))


def get_img_by_key(adata, by):
    n_obs = adata.shape[0]
    shape = adata.obsm["X_spatial"].max(axis=0)[::-1].astype(int) + 1
    if np.prod(shape) != n_obs:
        shape = adata.obsm["X_spatial_lowres"].max(axis=0)[::-1].astype(int) + 1
    img = adata.obs[by].values.reshape(shape)
    return img


def paint_tissue(adata, by="log1p_total_intensity", key_added="in_tissue"):
    img = get_img_by_key(adata, by)
    
    contours = find_contours(img, fully_connected="high")
    mask = np.ones_like(img, dtype=bool)
    print(img.shape)
    for contour in contours:
        c = close_contour(contour, *img.shape)
        mask[grid_points_in_poly(img.shape, c)] = False
    
    adata.obs[key_added] = mask.ravel()
    adata.obs[key_added] = adata.obs[key_added].astype("category")

    
def paint_tissue_fast(adata, by="log1p_total_intensity", key_added_in="in_tissue", key_added_out="background", sigma=3, threshold=0.3):
    img = get_img_by_key(adata, by)
    rmin, rmax = img.min(), img.max()
    if (rmax - rmin) > 100:
        print("Data range is very large, going to sublinear transform")
        img = np.sqrt(img)
        print(img.min(), img.max())
        
    blurred = gaussian(img, sigma=sigma)
    p = np.percentile(blurred, q=threshold*100)
    mask = morphology.remove_small_holes(
        morphology.remove_small_objects(
            blurred > p, 500
        ), 500
    )
    adata.obs[key_added_in] = mask.ravel().astype(int)
    adata.obs[key_added_out] = (~mask.ravel()).astype(int)


def shape_from_xy(arr):
    assert len(arr.shape) > 1
    return tuple(arr.max(axis=np.argmax(arr.shape))[::-1] + 1)


def detect_nuclei(data_1ch):
    model = StarDist2D.from_pretrained("2D_versatile_fluo")
    labels, details = model.predict_instances(normalize(data_1ch))
    return np.vstack(np.nonzero(labels)).T


def load_spillmat(infile=None):
    if not infile:
        infile = SPILLMAT_CSV
    return pd.read_csv(infile, index_col=0)


def align_spillmat(*, spillmat=None, input_metals=None):
    if spillmat is None:
        spillmat = load_spillmat()

    if input_metals is None:
        input_metals = set(spillmat.index.union(spillmat.columns))

    sm = spillmat.reindex(index=input_metals, columns=input_metals, fill_value=0)
    filled = sm.values
    np.fill_diagonal(filled, 1.0)
    return pd.DataFrame(filled, index=sm.index, columns=sm.columns)


def compensate_long(long, spillmat):
    comp_ = long @ np.linalg.inv(spillmat.T)
    comp_ = np.clip(comp_, 0, comp_.max()).astype(np.float32)
    return comp_


def _load_imc_acquisition(mcd_file, acq_label, get_metadata=True, preview=False):
    parser = McdParser(mcd_file)
    
    if acq_label is None:
        acq_label = parser.session.acquisition_ids[0]
    
    if isinstance(acq_label, int):
        aid = acq_label
        acq = parser.session.acquisitions.get(aid)

    else:
        for aid, acq in parser.session.acquisitions.items():
            if acq.description == acq_label:
                print(f"Found {acq_label} with ID {aid}")
                break
        else:
            raise ValueError(f"Label {acq_label} not found in {mcd_file}.")
        
    preab_image = None
    if acq.has_before_ablation_image:
        preab_image = imread(BytesIO(
            parser.get_before_ablation_image(aid)
        ))

    if get_metadata:
        xml = McdXmlParser(parser.get_mcd_xml(), mcd_file).metadata
        
    #raw = parser._get_acquisition_raw_data(acq)
    raw = parser.get_acquisition_data(aid)
    if preview:
        plt.imshow(raw.get_image_by_label("DNA3"), cmap="Blues")
    return acq, raw.image_data, preab_image, xml


def _coarse_grain(data, window_size=10):
    from skimage.util.shape import view_as_blocks
    
    window = (window_size, window_size, 1)
    pad_x = window_size - (data.shape[0] % window_size)
    pad_y = window_size - (data.shape[1] % window_size)
    
    padded = np.pad(data, ((0, pad_x), (0, pad_y), (0, 0)))
    
    blocks = view_as_blocks(padded, window)
    # expected to be x, y, c because data has been rolled prior to this function
    ndim = blocks.ndim
    
    return blocks.sum(axis=tuple(range(ndim-3, ndim)))


def clip_hot_pixels(img, hp_filter_shape, hp_threshold=None):
    """
    See
    https://github.com/BodenmillerGroup/ImcPluginsCP/blob/a53bb7e1dea60b859d57677ea9a15281fa84d493/plugins/smoothmultichannel.py#L417
    """
    if hp_filter_shape[0] % 2 != 1 or hp_filter_shape[1] % 2 != 1:
        raise ValueError(
            "Invalid hot pixel filter shape: %s" % str(hp_filter_shape)
        )
    if hp_threshold is None:
        hp_threshold = np.percentile(img, q=98)

    hp_filter_footprint = np.ones(hp_filter_shape)
    hp_filter_footprint[
        int(hp_filter_shape[0] / 2), int(hp_filter_shape[1] / 2)
    ] = 0

    max_img = maximum_filter(
        img, footprint=hp_filter_footprint, mode="reflect"
    )
    hp_mask = img - max_img > hp_threshold
    img = img.copy()
    img[hp_mask] = max_img[hp_mask]
    return img


def clean_channel(img, hp_filter_shape=(9, 9), noise_blob_radius=3):
    # this should be done on total intensity, not per channel
    #wh = white_tophat(img, selem=disk(noise_blob_radius))
    #cleaned = img - wh
    #return clip_hot_pixels(cleaned, hp_filter_shape)
    return clip_hot_pixels(img, hp_filter_shape)


def mcd_to_anndata(
    mcd_file: Union[str, Path],
    library_id: str,
    *,
    acquisition_label: Union[str, int, None] = None,
    compensate: bool = True,
    remove_hot_pixels: bool = True,
    remove_unused_channels: bool = True,
    preview: bool = False
):
    """
    Converts an MCD acqusition to an AnnData object
    
    """
    mcd_file = Path(mcd_file)
    assert mcd_file.exists()
    
    acquisition, raw_data, preab_image, metadata = _load_imc_acquisition(
        mcd_file, acquisition_label, get_metadata=True, preview=preview
    )
    raw_data = np.moveaxis(raw_data, 0, 2)
    
    # get long data
    # compensate
    compensated = None
    if compensate:
        print("Compensating")
        n_channels = acquisition.n_channels
        long = raw_data.reshape(-1, n_channels)
        spillmat = align_spillmat(input_metals=acquisition.channel_names)
        compensated = compensate_long(long, spillmat).reshape(raw_data.shape)
        
    # clean
    cleaned = None
    if remove_hot_pixels:
        print("Removing small blobs and hot pixels")
        x_ = raw_data
        if compensated is not None:
            x_ = compensated
        # only doing maximum filter --- not sure if the conway filter is needed/worth it
        cleaned = np.zeros_like(x_)
        for k in range(x_.shape[-1]):
            cleaned[..., k] = clean_channel(x_[...,k])
        tot_int = cleaned.sum(axis=-1)
        wh = white_tophat(tot_int, selem=disk(3))
        mask = np.where((tot_int - wh) == 0)
        cleaned[mask] = 0

                      
    var_info = pd.DataFrame({
        "ab_mass": acquisition.channel_masses,
        "ab_label": acquisition.channel_labels,
        "ab_name": acquisition.channel_names
    })
    var_info["ab_label"] = var_info.ab_label.fillna(var_info.ab_name)
    var_info.set_index("ab_label", inplace=True)

    print("Finding nuclei")
    nuc_inds = np.where(var_info.index.str.startswith("DNA"))
    nuc_x = raw_data
    if cleaned is not None:
        nuc_x = cleaned
    elif compensated is not None:
        nuc_x = compensated
    nuc_channel = nuc_x[..., nuc_inds].sum(axis=-1)
    nuc_points = detect_nuclei(nuc_channel)
    
    # filter out garbage channels:
    drop_inds = np.logical_or(
        var_info.ab_name.str.extract("([A-z]+)(\d+)", expand=True).agg(lambda r: "".join(r[::-1]), axis=1) == var_info.index,
        var_info.ab_name == var_info.index
    )
    var_info = var_info.loc[~drop_inds, :]
    keep_inds = np.where(~drop_inds)
    raw_data = raw_data[..., keep_inds]
    
    adata = sc.AnnData(
        raw_data.reshape(-1, len(var_info)),
        var = var_info,
        obs = pd.DataFrame(index=pd.RangeIndex(stop=np.prod(raw_data.shape[:-1])))
    )
    sc.pp.calculate_qc_metrics(adata, expr_type="intensity", var_type="antibodies", percent_top=None, inplace=True)
    
    # We need the layers:
    if compensated is not None:
        adata.layers["compensated"] = compensated[..., keep_inds].reshape(adata.shape)
    if cleaned is not None:
        adata.layers["cleaned"] = cleaned[..., keep_inds].reshape(adata.shape)

    # add mcd metadata
    adata.uns["mcd"] = metadata
    
    # Add in spatial data
    ys, xs = np.meshgrid(np.arange(raw_data.shape[1]), np.arange(raw_data.shape[0]))
    adata.obsm["X_spatial"] = np.vstack((ys.ravel(), xs.ravel())).T
    # We consider the pre-ablation image the 'hi-res' image
    adata.uns["spatial"] = {
        library_id: dict(
            images=dict(hires=preab_image),
            scalefactors=dict(spot_diameter_fullres=1, tissue_hires_scalef=1)
        )
    }

    nuclei_counts = np.zeros(raw_data.shape[:-1])
    nuclei_counts[nuc_points[:,0], nuc_points[:,1]] = 1
    
    adata.obs["nuclei_counts"] = nuclei_counts.ravel()

    # ID tissue
    paint_tissue_fast(adata)

    return adata


def coarse_grain_imc_adata(
    adata: sc.AnnData,
    *,
    layer: Union[str, None] = None,
    window_size: int = 10,
):
    img_shape = shape_from_xy(adata.obsm["X_spatial"])
    img_shape += (adata.shape[1],)
    if layer is not None:
        x_ = adata.layers[layer].reshape(img_shape)
    else:
        x_ = np.asarray(adata.X).reshape(img_shape)
    
    coarse = _coarse_grain(x_, window_size)
    nuc_counts = _coarse_grain(adata.obs["nuclei_counts"].values.reshape(img_shape[:-1])[...,None], window_size)
    
    coarse_adata = sc.AnnData(
        coarse.reshape(-1, adata.shape[1]),
        var = adata.var,
        obs = pd.DataFrame(index=pd.RangeIndex(stop=np.prod(coarse.shape[:-1]))),
        uns = adata.uns.copy()
    )
    sc.pp.calculate_qc_metrics(coarse_adata, expr_type="intensity", var_type="antibodies", percent_top=None, inplace=True)

    # Add in spatial data
    ys, xs = np.meshgrid(np.arange(coarse.shape[1]), np.arange(coarse.shape[0]))
    coarse_adata.obsm["X_spatial"] = (np.vstack((ys.ravel(), xs.ravel())).T * window_size) + window_size/2
    coarse_adata.obsm["X_spatial_lowres"] = np.vstack((ys.ravel(), xs.ravel())).T
    # We consider the pre-ablation image the 'hi-res' image
    lib_id = list(coarse_adata.uns["spatial"].keys())[0]

    paint_tissue_fast(coarse_adata)
    
    coarse_adata.uns["spatial"][lib_id]["scalefactors"] = dict(
        spot_diameter_fullres=window_size, tissue_hires_scalef=1, tissue_lowres_scalef=1/window_size
    )
    
    coarse_adata.obs["nuclei_counts"] = nuc_counts.ravel()
    
    return coarse_adata


def preprocess_imc_data(adata, vst_cofactor=5.):
    adata_ = adata[adata.obs.in_tissue.astype(bool), :].copy()
    sc.pp.filter_cells(adata_, min_counts=50)
    
    adata_.layers["cleaned"] = adata_.X.copy()
    adata_.X = np.arcsinh(adata_.X / vst_cofactor)
    adata_.layers["vst"] = adata_.X.copy()
    
    maxs = np.percentile(adata_.X, q=98, axis=0)
    adata_.X = np.clip(adata_.X, 0, maxs)
    adata_.layers["vst-clipped"] = adata_.X.copy()
    
    print("filtering done")
    redux_vars = adata_.var_names[~adata_.var_names.str.startswith("DNA")]
    adata_.obsm["X_redux"] = adata_[:, redux_vars].X.copy()
    sc.pp.neighbors(adata_, n_neighbors=15, use_rep="X_redux", metric="correlation")
    print("neighbors done")
    sc.tl.umap(adata_, min_dist=0.5)
    sc.tl.leiden(adata_, resolution=1)
    return adata_


def convert_imc_adata_to_text(
    adata, 
    outpath=None,
    cluster_key="leiden", 
    cluster_prefix="Cluster_", 
    layer=None
):
    df = sc.get.obs_df(adata, adata.var_names.tolist(), layer=layer)
    
    coords = pd.DataFrame(adata.obsm["X_spatial_lowres"], columns=list("XY"))
    coords["Z"] = 0 
    
    unknown = pd.DataFrame(
        np.zeros_like(coords.values, dtype=np.uint8),
        columns=["Start_push", "End_push", "Pushes_duration"],
        index=df.index
    )
    
    coords.index = df.index
    clusters = pd.get_dummies(adata.obs[cluster_key], prefix=cluster_prefix, prefix_sep="")
    res = pd.concat((unknown, coords, df, clusters), axis=1)
        
    # now we need to make this rectangular
    shape = res[["X","Y"]].max(axis=0).values + 1
    size = np.prod(shape)
    rect = np.ones(shape)
    rect[res.X, res.Y] = 0
    nonzeros = np.nonzero(rect)
    
    empty = pd.DataFrame(np.zeros((len(nonzeros[0]), len(res.columns))), columns=res.columns)
    empty["X"] = nonzeros[0]
    empty["Y"] = nonzeros[1]
    final = pd.concat((res, empty), axis=0)
    #final = final[unknown.columns.tolist() + ["Y", "X"] + final.columns[5:].tolist()]
    final = final.sort_values(["Y", "X"])
    final["Z"] = np.arange(len(final), dtype=np.uint16)
    
    for col in final.columns:
        final[col] = pd.to_numeric(final[col], downcast="unsigned")
    
    if outpath is not None:
        final.to_csv(outpath, index=False, sep="\t")
    return final
