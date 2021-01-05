import numpy as np
from matplotlib.colors import Normalize

from scanpy.plotting._anndata import _prepare_dataframe


def poor_mans_matrixplot(ax, adata, var_names, groupby, standard_scale=None, use_raw=False, log=False, layer=None, swap_axes=False, **kwds):
    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw,
        log,
        layer=layer,
    )

    if standard_scale == 'group':
        obs_tidy = obs_tidy.sub(obs_tidy.min(1), axis=0)
        obs_tidy = obs_tidy.div(obs_tidy.max(1), axis=0).fillna(0)
    elif standard_scale == 'var':
        obs_tidy -= obs_tidy.min(0)
        obs_tidy = (obs_tidy / obs_tidy.max(0)).fillna(0)
    elif standard_scale is None:
        pass

    mean_obs = obs_tidy.groupby(level=0).mean()
    if swap_axes:
        mean_obs = mean_obs.T
        yticklabels = mean_obs.index
        xticklabels = categories
        ylabel = ""
        xlabel = groupby
    else:
        xticklabels = mean_obs.columns
        yticklabels = categories
        xlabel = ""
        ylabel = groupby
    nx = mean_obs.shape[1]
    ny = mean_obs.shape[0]

    normalize = Normalize(
        vmin=kwds.get('vmin'), vmax=kwds.get('vmax')
    )

    pc = ax.pcolor(mean_obs, edgecolor='gray', **kwds)

    x_ticks = np.arange(nx) + 0.5
    y_ticks = np.arange(ny) + 0.5
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticklabels(xticklabels, rotation=90)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.tick_params(axis='both', labelsize='small')
    ax.grid(False)

    ax.set_xlim(0, nx)
    ax.set_ylim(ny, 0)
