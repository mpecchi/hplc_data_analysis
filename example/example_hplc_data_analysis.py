# %%
from hplc_data_analysis.hplc import Project

# this has to be changed to where the _example/data folder is
folder_path = r"/Users/matteo/Projects/hplc_data_analysis/example/data"
# folder_path = plib.Path(plib.Path(__file__).cwd(), 'data')
# change this to folder_path = plib.Path(r"C:\Path\To\Your\Data") for your project
# class methods need to be called at the beginning to influence all instances

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


for level, param in zip(["files", "repliactes", "samples"], hplc.acceptable_params):
    hplc.plot_report(r)
