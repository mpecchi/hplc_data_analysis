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
folder_path = plib.Path(plib.Path.cwd(), 'data')
# change this to folder_path = plib.Path(r"C:\Path\To\Your\Data") for your project
# class methods need to be called at the beginning to influence all instances
Project.set_folder_path(folder_path)
Project.set_plot_grid(False)
Project.set_plot_font('Sans')  # ('Times New Roman')

p = Project(rebuild_compounds_properties=False)
# %%
# p.create_compounds_properties()
p.create_files_replicates_samples()
a = p.files
aa = p.replicates[0]
aaa = p.replicates[1]
b = p.samples[0]
bb = p.samples_std[0]

p.create_all_statistics()
c, d = p.create_param_report()
p.save_files_report()
p.save_param_report()
p.save_param_aggrrep()
zz, zzstd = p.return_param_aggrrep(param='f_sample')

p.plot_ave_std(param='f_sample', aggr=True, min_y_thresh=0,
    color_palette='Set2')
p.plot_ave_std(min_y_thresh=200, legend_location='outside',
                only_samples_to_plot=['FW250C1h1', 'FWCP250C1h1'],
                )
```
<!-- EXAMPLE_END -->
