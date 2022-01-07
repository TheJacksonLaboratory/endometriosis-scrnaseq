import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cmocean as cmo

red_colormap = color_map = sns.blend_palette(
    ["0.9", sns.xkcd_rgb["bright red"]], as_cmap=True
)
tab20a = sns.mpl_palette("tab20_r", 20)
tab20b = sns.mpl_palette("tab20b_r", 20)
tab20c = sns.mpl_palette("tab20c_r", 20)
main_palette = (
)

heatmap_cmap = cmo.tools.cmap(cmo.cm.balance(np.linspace(0, 1, 256), 0.9))
