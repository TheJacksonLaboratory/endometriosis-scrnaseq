import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cmocean as cmo

from matplotlib.colors import hex2color

red_colormap = color_map = sns.blend_palette(
    ["0.9", sns.xkcd_rgb["bright red"]], as_cmap=True
)

main_palette = [
	'#828282', '#8da8ea', '#9dbf98', '#fbb2ba', '#e97367'
]

mye_palette = [
	'#6b6ecf', '#fde0dd', '#a494eb', '#fbb2ba', '#f883ab', '#ea4d9c',
	'#b22fad', '#d79fd5', '#ad017e', '#600070', '#d6616b', '#ead5e9',
	'#a55194', '#843c39', '#c49c94'
]

lym_palette = [
	'#327925', '#4c8a41', '#669c5d', '#81ad79', '#9dbf98', '#b8d1b4',
	'#d2e2d0', '#b7ced3', '#7da9b3', '#6297a3', '#468593', '#2a7483',
	'#7fcdb6', '#637939'
]

epi_palette = [
	'#256fa8', '#5890bb', '#8cb1ce', '#9ecae1', '#e5eff9', '#bfcdee',
	'#8da8ea', '#5b84e6', '#08488e', '#758da3'
]

stro_palette = [

	'#e53524', '#e97367', '#ec948b', '#f2f1f1', '#ebd8cf', '#dba589',
	'#c55a23', '#843c39', '#feb852', '#e6550d', '#fd8d3c', '#fdd0a2'
]

endo_palette = [
	'#978a84', '#2b2b2b', '#828282', '#f0d2cf', '#c8aca9', '#5c5c5c',
	'#516572', '#758da3'
]

org_palette = [
	'#ededed', '#256fa8', '#5890bb', '#8cb1ce', '#9ecae1', '#e5eff9',
	'#bfcdee', '#8da8ea', '#5b84e6', '#08488e', '#758da3'
]

heatmap_cmap = cmo.tools.cmap(cmo.cm.balance(np.linspace(0, 1, 256), 0.9))

sample_type_color_mapping = {
    "Ctrl": "#9dcbb9",
    "EcP": "#2c387a",
    "EcPA": "#4f699e",
    "EcO": "#757097",
    "EuE": "#34693D",
}

celltype_color_mapping = {
    "epithelial": "#8da8ea", 
    "stromal": "#e97367", 
    "endothelial": "#828282",
    "myeloid": "#fbb2ba", 
    "lymphocytes": "#9dbf98",
}

sample_type_palette = [hex2color(v) for v in sample_type_color_mapping.values()]
patient_palette = sns.color_palette("husl", 18, desat=0.8)
celltype_palette = [hex2color(v) for v in celltype_color_mapping.values()]

heatmap_color = sns.diverging_palette(14, 215, as_cmap=True)
