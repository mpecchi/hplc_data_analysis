# hplc_data_analysis

## A Python tool to manage multiple HPLC qualitative tables and automatically split chemicals into functional groups.

An open-source Python tool that can automatically:
- handle multiple HPLC semi-quantitative data tables (obtained with different methods)
- duild a database of all identified compounds and their relevant properties using PubChemPy
- split each compound into its functional groups using a published fragmentation algorithm
- produce single file report, single replicate (intended as the sum of more methods applied to one vial), comprehensive multi-sample (average and deviation of replicates) reports and aggregated reports based on functional group mass fractions in the samples
- provides plotting capabilities

## Naming convention for samples

To ensure the code handles replicates of the same sample correctly, names have to follow the convention:
*method_name-of-sample-with-dashes-only_replicatenumber*
(for now, method can either be 210 or 254, as the wavelenght analyzed using a UVvis detector

Examples that are *correctly* handled:
- 210_Bio-oil-foodwaste-250C_1
- 254_Bio-oil-foodwaste-250C_1
- 210_FW_2
- 254_FW_2

Examples of *NON-ACCEPTABLE* names are
- 210-bio_oil_1
- 254-FW1

## Example

A comprehensive example is provided on the GitHub repository to show how inputs should be formatted.
To test the module, install the `hplc_data_analysis` module, download the example folder given in the repository, and run the example_hplc_data_analysis.py. The folder_path needs to be set to where your data folder is.

The example code is shown here for convenience:
<!-- EXAMPLE_START -->
```python
# example of hplc_data_analysis
import pathlib as plib
from hplc_data_analysis import Project
# this has to be changed to where the _example/data folder is
folder_path = plib.Path(plib.Path.cwd().parent, 'example/data')
# folder_path = plib.Path(plib.Path(__file__).cwd(), 'data')
# change this to folder_path = plib.Path(r"C:\Path\To\Your\Data") for your project
# class methods need to be called at the beginning to influence all instances
Project.set_folder_path(folder_path)
Project.set_plot_grid(False)
Project.set_plot_font('Sans')  # ('Times New Roman')

p = Project(rebuild_compounds_properties=False)
# %%
# p.create_compounds_properties()
# get files, replicates, and samples info
files_info = p.files_info
replicates_info = p.replicates_info
samples_info = p.samples_info
# create files, replicates and samples dfs
p.create_files_replicates_samples()
#%%
# example of file
FWCP250C1h1_1_210 = p.files['210_FWCP250C1h1_1']
# examples of replicates (merged files with diff. methods)
FWCP250C1h1_1 = p.replicates['FWCP250C1h1_1']
FWCP250C1h1_2 = p.replicates['FWCP250C1h1_2']
# examples of samples (average of replicates)
FWCP250C1h1 = p.samples['FWCP250C1h1']
FWCP250C1h1_std = p.samples_std['FWCP250C1h1']
# add all statistics to info dfs
p.update_all_info_statistics()
# create report and aggrreports for different parameters
# for files, replicates and samples
files_report_conc = \
    p.create_files_param_report(param='conc_vial_if_undiluted_mg_L')
replicates_report_conc = \
    p.create_replicates_param_report(param='conc_vial_if_undiluted_mg_L')
replicates_report_fr = \
    p.create_replicates_param_report(param='fraction_of_sample_fr')
replicates_aggrrep_fr = \
    p.create_replicates_param_aggrrep(param='fraction_of_sample_fr')
samples_report_fr, samples_report_fr_std = \
    p.create_samples_param_report(param='fraction_of_sample_fr')
samples_aggrrep_fr, samples_report_fr_std = \
    p.create_samples_param_aggrrep(param='fraction_of_sample_fr')
p.save_files_info()
p.save_replicates_info()
p.save_samples_param_report()
p.save_samples_param_aggrrep()
# %%
p.plot_ave_std(param='conc_vial_mg_L', aggr=True, min_y_thresh=0,
    y_lim=[0, 50000],legend_location='outside',
    color_palette='Set2')
# %%
p.plot_ave_std(param='conc_vial_if_undiluted_mg_L', replicates_or_samples='replicates',
    min_y_thresh=0, legend_location='outside', xlab_rot=20,
    # only_samples_to_plot=['FW250C1h1', 'FWCP250C1h1'],
    y_lim=[0, 50000], annotate_outliers=True)
# %%
p.plot_ave_std(replicates_or_samples='samples', param='fraction_of_sample_fr',
    aggr=True, min_y_thresh=2000, legend_location='outside', xlab_rot=20,
    # only_samples_to_plot=['FW250C1h1', 'FWCP250C1h1'],
    y_lim=[0, 50000], annotate_outliers=True)
```
<!-- EXAMPLE_END -->
