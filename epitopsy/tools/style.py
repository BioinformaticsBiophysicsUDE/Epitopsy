import os
import numpy as np
#import pylab as py

style = {}

#### plot parameters
style["fontsize"] = 25
style["ticksize"] = 15
style["linewidth"] = 2.5
style["markeredgewidth"] = 2.5
style["dpi"] = 300
style["default_pymol_color"] = "gray80"

## DI plot parameters
eukaryotes_name = "eukaryote"
viral_name = "viral"
bacterial_name = "bacterial"
style["DI_protein_types"] = [ bacterial_name, eukaryotes_name, viral_name]
style["subplot_bottom"] = 0.15
style["{0}_color".format(bacterial_name)] = "green"
style["{0}_alpha".format(bacterial_name)] = 1
style["{0}_color".format(viral_name)] = "red"
style["{0}_alpha".format(viral_name)] = 0.75
style["{0}_color".format(eukaryotes_name)] = "blue"
style["{0}_alpha".format(eukaryotes_name)] = 0.75
style["pw_color"] = "black"
style["msa_color"] = "red"
style["pure_di_color"] = "green"
style["pure_mi_color"] = "orange"#"purple"
style["pw_msa_alpha"] = 0.5

#### path settings
struct_folder = os.path.join(os.getenv("HOME"),"Dropbox","my_thesis",
    "data", "structures")
style["struct_folder"] = struct_folder
style["ligand_dict"] = {"HS":os.path.join(struct_folder,"HS_3E7J","HS_3E7J.pqr"),
        "H":os.path.join(struct_folder,"heparin_1GMO","heparin_1GMO_A.pqr"),
        "HSsulf":os.path.join(struct_folder,"HS","HSsulf.pqr"),
        "HSunsulf":os.path.join(struct_folder,"HS","HSunsulf.pqr"),
        "CA":os.path.join(struct_folder,"CA","ca.pqr"),
        "SO4":os.path.join(struct_folder,"sulfate","sulfate.pqr"),
        "H2O":os.path.join(struct_folder,"water","water.pqr"),
        "H2Osmall":os.path.join(struct_folder,"water","small_water.pqr")}
style["pymol_align"] = os.path.join(os.getenv("HOME"),"workspace",
        "epitopsyPython","src","epitopsy","tools","pymol_align.py")
style["energy_box_file_path"] = "gibbs_free_energy.dx"
style["counter_box_file_path"] = "counter_matrix.dx.gz"
# ligand accessible surface
style["LAS_box_file_path"] = "gibbs_free_energy_LAS.dx"
style["plot_box_file_path"] = style["LAS_box_file_path"]

## ala scan analysis
style["energy_func"] = "normal_fav_surface_score"
style["analysis_funcs"] = ["DG",
    "LAS_fav_surface_score",
    "LAS_fav_volume_score",
    "LAS_fav_surface",
    "LAS_fav_volume",
    "normal_fav_surface_score",
    "normal_fav_volume_score",
    "normal_fav_surface",
    "normal_fav_volume"]
style["alascan_ylabel"] = {"LAS_fav_surface":r"$A_{+}$"+ r"$<-1k_{B}T$",
        "LAS_fav_surface_score":r"$wA_{+}$"+ r"$<-1k_{B}T$",
        "LAS_fav_volume":r"$V_{+}$"+ r"$<-1k_{B}T$",
        "LAS_fav_volume_score":r"$wV_{+}$"+ r"$<-1k_{B}T$",
        "normal_fav_surface":r"$A_{+}$"+ r"$<-1k_{B}T$",
        "normal_fav_surface_score":r"$wA_{+}$"+ r"$<-1k_{B}T$",
        "normal_fav_volume":r"$V_{+}$"+ r"$<-1k_{B}T$",
        "normal_fav_volume_score":r"$wV_{+}$"+ r"$<-1k_{B}T$",
        "DG":r"$\Delta G_{\mathit{bind}}$",
        "k_D":r"$k_{D}$"}


## hedgehog path settings
style["ca_color"] = "gray70"
style["zn_color"] = "black"
style["hh_save_folder"] = os.path.join(os.getenv("HOME"),"projects",
        "Project3-Electrostatics", "Docking", "esp_based", "sonic_hedgehog",
        "publication", "plots")
style["hh_full_sulfur_energy_level"] = 1
style["hedgehog_pdb_id_HIP"] = {"1vhh": [135],
        "2wfq": [135],
        "2wfx": [135],
        "2wg4": [135],
        "3d1m": [135],
        "3ho5": [134],
        "3m1n": [134],
        "3m1n-1": [134],
        "3m1n-2": [134],
        "3mxw": [134],
        "3n1r": [135]}
# shh
style["shh_color"] = "blue"
style["shh_pymol_color"] = "blue"
style["shh_id_shift"] = 0
# ihh
style["ihh_color"] = "Brown"
style["ihh_pymol_color"] = "brown"
style["ihh_id_shift"] = 3
# dhh
style["dhh_color"] = "GoldenRod"
style["dhh_pymol_color"] = "sand"
style["dhh_id_shift"] = -2
# drosoh
style["hh_color"] = "gray"
style["hh_alpha"] = 0.75
style["drosoh_pymol_color"] = "gray"
style["drosoh_id_shift"] = 60

## chemokine path settings
style["chemokine_save_folder"] = os.path.join(os.getenv("HOME"), "projects",
        "Project3-Electrostatics", "Docking", "esp_based", "Mip",
        "publication", "plots")
style["cc_full_sulfur_energy_level"] = 2
style["CCL3_color"] = "cyan"
style["CCL4_color"] = "red"

# verify settings
style["verify_energy_cutoff"] = -2

## DI path settings
if os.path.split(os.getenv("HOME"))[1] == "wilms":
    style["data_path"] = os.path.join("/data", "wilms", "projects",
        "Project6-Direct_Information", "protein_data")
    style["depreciated_data_path"] = os.path.join("/data", "wilms", "projects",
        "Project6-Direct_Information", "protein_data_depreciated")
    style["{0}_path".format(bacterial_name)] = os.path.join("/data",
        "wilms", "projects", "Project6-Direct_Information",
        "bacterial_proteins")
    style["{0}_path".format(viral_name)] = os.path.join("/data",
          "wilms", "projects", "Project6-Direct_Information",
          "viral_proteins")
    style["{0}_path".format(eukaryotes_name)] = os.path.join("/data",
          "wilms", "projects", "Project6-Direct_Information",
          "eukaryote_proteins")
    style["interesting_path"] = os.path.join("/data",
          "wilms", "projects", "Project6-Direct_Information",
          "interesting_proteins")
    style["di_plot_path"] = os.path.join("/data", "wilms", "projects",
        "Project6-Direct_Information", "plots")
    style["additional_data"] = os.path.join("/data", "wilms", "projects",
        "Project6-Direct_Information", "additional_data")
else:
    style["data_path"] = os.path.join(os.getenv("HOME"), "projects",
        "Project6-Direct_Information", "protein_data")
    style["depreciated_data_path"] = os.path.join(os.getenv("HOME"), "projects",
        "Project6-Direct_Information", "protein_data_depreciated")
    style["{0}_path".format(bacterial_name)] = os.path.join(
        os.getenv("HOME"), "projects", "Project6-Direct_Information",
        "bacterial_proteins")
    style["{0}_path".format(viral_name)] = os.path.join(os.getenv("HOME"),
          "projects", "Project6-Direct_Information", "viral_proteins")
    style["{0}_path".format(eukaryotes_name)] = os.path.join(
        os.getenv("HOME"), "projects", "Project6-Direct_Information",
        "eukaryote_proteins")
    style["interesting_path"] = os.path.join(os.getenv("HOME"),
          "projects", "Project6-Direct_Information",
          "interesting_proteins")
    style["di_plot_path"] = os.path.join(os.getenv("HOME"), "projects",
        "Project6-Direct_Information", "plots")
    style["additional_data"] = os.path.join(os.getenv("HOME"), "projects",
        "Project6-Direct_Information", "additional_data")

#### various settings

## DI performance settings
style["dca_new_best_theta"] = 0.4
style["dca_old_best_theta"] = 0.1
style["mi_new_best_theta"] = 0.4
style["mi_old_best_theta"] = 0.1
style["dca_new_best_identity"] = 0.6
style["dca_old_best_identity"] = 0.9
style["mi_new_best_identity"] = 0.6
style["mi_old_best_identity"] = 0.9
style["seq_dist"] = 5  # min sequence separation between coevoling residues
style["max_dist"] = 8  # max distance between coevolving residues


#def color_boxplot(box_data, color, alpha, linewidth):
#    black = "black"
#    red = "red"
#    blue = "blue"
#    medians_color = "black"
#    if medians_color == color:
#        raise ValueError("Colors are identical! medians: {0} == box: {1}".format(medians_color, color))
#    py.setp(box_data["boxes"],color=color, alpha=alpha, edgecolor=black)
#    py.setp(box_data["caps"], linewidth=linewidth,color=color, alpha=alpha)
#    py.setp(box_data["whiskers"], linewidth=linewidth,color=color, alpha=alpha)
#    py.setp(box_data["medians"], linewidth=linewidth, color=medians_color)
#    py.setp(box_data["fliers"], markersize=10, markeredgewidth=2, color=black)
