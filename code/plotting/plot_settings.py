import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import cmocean as cmo
import scanpy as sc


matplotlib.rcParams["figure.facecolor"] = "white"
matplotlib.rcParams["figure.dpi"] = 300
matplotlib.rcParams["figure.subplot.left"] = 0
matplotlib.rcParams["figure.subplot.right"] = 1
matplotlib.rcParams["legend.fontsize"] = 8.88
matplotlib.rcParams["savefig.format"] = "pdf"
matplotlib.rcParams["xtick.labelsize"] = 8.0
matplotlib.rcParams["ytick.labelsize"] = 8.0
matplotlib.rcParams["axes.titlesize"] = 10.0
matplotlib.rcParams["axes.labelsize"] = 10.0
matplotlib.rcParams["savefig.pad_inches"] = 0
matplotlib.rcParams["font.family"] = ["Arial"]
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

sc.settings._vector_friendly = False
