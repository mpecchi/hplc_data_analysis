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