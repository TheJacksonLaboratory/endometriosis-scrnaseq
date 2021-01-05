import numpy as np
from matplotlib.text import Annotation
from matplotlib.legend_handler import HandlerBase
from matplotlib import patheffects
from matplotlib.colors import CSS4_COLORS, hex2color


class AnnotationHandler(HandlerBase):
    def create_artists(self, legend, artist, xdescent, ydescent,
                       width, height, fontsize, trans):
        a = Annotation(
            artist.get_text(), [width/2, height/2],
            color=artist.get_color(),
            xycoords=trans, fontsize=fontsize*0.5,
            ha="center", va="center",
            bbox=dict(boxstyle="circle", facecolor=artist.get_bbox_patch().get_facecolor())
        )
        #a.set_path_effects(artist.get_path_effects())
        a.set_label(artist.get_label())
        return [a]


def add_numbered_circles_to_umap(umap_ax, circle_prefix="", circle_bg=None, circle_kwargs={}, min_cells=400, text_color="black"):
    """
    Works on a single axes at a time
    : circle_bg : If `None`, uses the color of the groups.  Otherwise, expects a color, like `'white'`.
    : circle_kwargs : A dict to change the circle params (like `linewidth` or `edgecolor`)
    : text_color : Can be one of `None`, `"white"` or `"black"`.
    """
    allowed_colors = ["black", "white", None]
    assert text_color in allowed_colors, f"Text color needs to be in {allowed_colors}. You gave {text_color}"
    if isinstance(circle_bg, str):
        circle_bg = hex2color(CSS4_COLORS[circle_bg])

    # get labels
    leg = umap_ax.get_legend()
    labels = [_.get_label() for _ in leg.legendHandles]
    facecolors = np.array([_.get_facecolors()[0] for _ in leg.legendHandles])[:,:3]
    leg.remove()
    pcol = None
    for pc in umap_ax.collections:
        if pc.get_offsets().shape[0] == pc.get_facecolors().shape[0]:
            pcol = pc
            break
    else:
        raise Exception("Can't figure out which points to use")

    color_vector = pcol.get_facecolors()[:,:3]
    x, y = pcol.get_offsets().T
    circle_params = dict(edgecolor="black", linewidth=0.5)
    circle_params.update(**circle_kwargs)

    new_handles = []
    for k, (label, group) in enumerate(zip(labels, facecolors), start=1):
        inds = (color_vector == group).all(axis=1)
        mx, my = np.median(x[inds]), np.median(y[inds])
        bg = group if circle_bg is None else circle_bg
        if not text_color:
            fontcolor = "black" if bg[:3].dot([0.299, 0.587, 0.114]) > 150/256 else "white"
        else:
            fontcolor = text_color
        stroke_color = allowed_colors[int(fontcolor == "black")]
        dx = -0.5 if inds.sum() < min_cells else 0
        dy = 0.8 if inds.sum() < min_cells else 0

        text = f"{circle_prefix}{k}"
        #text = f"{text:^4}"
        h = umap_ax.annotate(text, [mx+dx, my+dy], xycoords="data", fontsize=12, ha="center", va="center", color=fontcolor,
               bbox=dict(boxstyle="circle", facecolor=bg, **circle_params))
        #hw = umap_ax.annotate(text, [mx+dx, my+dy], xycoords="data",fontsize=8, ha="center", va="center", color=fontcolor,
        #       bbox=dict(boxstyle="circle", facecolor="white", **circle_params))
        h.set_label(label)
        #hw.set_label(label)
        h.set_path_effects([patheffects.withStroke(linewidth=1, foreground=stroke_color)])
        #hw.set_path_effects([patheffects.withStroke(linewidth=1, foreground=stroke_color)])
        new_handles.append(h)

    umap_ax.legend(new_handles, labels, bbox_to_anchor=(1, 0.5), loc="center left", frameon=False, fontsize="small",
                   handler_map={Annotation: AnnotationHandler()}, handletextpad=0.5)


def strip_legend(ax):
    leg = ax.get_legend()
    h = leg.legendHandles
    labs = [_.get_label() for _ in h]
    leg.remove()
    return h, labs
