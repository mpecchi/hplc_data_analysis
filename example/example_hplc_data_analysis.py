# %%
import pathlib as plib
from hplc_data_analysis.hplc import Project

# this is the relative path to the folder where the data is
# which is the same as where the script is __file__
folder_path = plib.Path(__file__).resolve().parent
# if running as a Jupyter notebook, use absolute paths
# folder_path = r"absolute path to folder"
folder_path = r"/Users/matteo/Projects/hplc_data_analysis/example/data"
hplc = Project(folder_path)
# %%
# create files info (dilution and total sample concentration information need to be added manually)
files_info_created = hplc.create_files_info()
# %% load the files info with dilution data and concentration
files_info = hplc.load_files_info(update_saved_files_info=False)
# %%
# create replicates and samples info
replicates_info = hplc.create_replicates_info()
samples_info = hplc.create_samples_info()
# create the samples (intended as multiple replicates (multiple files))
hplc.create_samples()
# %%
# to get a list of all files, replicates and samples
list_of_all_filenames = files_info.index.tolist()
list_of_all_replicatenames = replicates_info.index.tolist()
list_of_all_samplenames = samples_info.index.tolist()
# %%
# this creates the compounds properties for all files
compounds_properties = hplc.create_compounds_properties(update_saved_files_info=False)
# %%
# once the compounds properties are created and saved, they can be loaded to save time
compounds_properties = hplc.load_compounds_properties()
# %%

# plot the chosen param for all files, names_to_keep is a list of the files to plot
mf = hplc.plot_report(
    files_replicates_or_samples="files",
    param="conc_vial_mg_L",
    names_to_keep=[
        "210_FW250C1h1_1",
        # "210_FW250C1h1_2",
        "210_FWCP250C1h1_1",
        "210_FWCP250C1h1_2",
        "254_FW250C1h1_1",
        "254_FW250C1h1_2",
        "254_FWCP250C1h1_1",
        "254_FWCP250C1h1_2",
        # "254_FWCPMn250C1h1_1",
    ],
    height=5,
    width=10,
    x_ticklabels_rotation=20,
    legend_ncols=2,
)
# %%
# plot replicates, which are the combination of files with
# different methods but same vial
mf = hplc.plot_report(
    files_replicates_or_samples="replicates",
    param="conc_vial_mg_L",
    names_to_keep=[
        "FW250C1h1_1",
        "FW250C1h1_2",
        "FWCP250C1h1_1",
        "FWCP250C1h1_2",
    ],
    y_axis_min_threshold=10000,
    height=5,
    width=10,
    x_ticklabels_rotation=20,
    legend_ncols=1,
)
# %%
# plot samples, which are the combination of replicates for the same
# material
mf = hplc.plot_report(
    files_replicates_or_samples="samples",
    param="conc_vial_mg_L",
    names_to_keep=["FW250C1h1", "FWCP250C1h1", "FWCPMn250C1h1"],
    y_axis_min_threshold=10000,
    x_ticklabels_rotation=20,
    legend_ncols=1,
)
# %%
mf = hplc.plot_report(
    files_replicates_or_samples="samples",
    param="fraction_of_sample_fr",
    names_to_keep=["FW250C1h1", "FWCP250C1h1", "FWCPMn250C1h1"],
    y_axis_min_threshold=10000,
    x_ticklabels_rotation=20,
    legend_ncols=1,
)
# %%
mf = hplc.plot_report(
    files_replicates_or_samples="samples",
    param="fraction_of_sample_fr",
    report_or_aggrrep="aggrrep",
    names_to_keep=["FW250C1h1", "FWCP250C1h1", "FWCPMn250C1h1"],
    y_axis_min_threshold=10000,
    x_ticklabels_rotation=20,
    legend_ncols=1,
)


# uncomment the following if you want to save everything as an excel file in the output folder
# hplc.save_files_samples_reports()

# %%
