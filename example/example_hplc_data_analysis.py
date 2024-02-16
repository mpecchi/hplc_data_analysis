# example of hplc_data_analysis
import pathlib as plib
from hplc_data_analysis import Project
# this has to be changed to where the _example/data folder is
folder_path = plib.Path(plib.Path.cwd().parent, 'example/data')
# folder_path = plib.Path(plib.Path(__file__).cwd(), 'data')
# change this to folder_path = plib.Path(r"C:\Path\To\Your\Data") for your project
# class methods need to be called at the beginning to influence all instances
# class methods need to be called at the beginning to influence all instances
Project.set_folder_path(folder_path)
Project.set_plot_grid(False)
Project.set_plot_font('Sans')  # ('Times New Roman')

p = Project(rebuild_compounds_properties=False)
# %%
# p.create_compounds_properties()
p.create_files_replicates_samples()
files_info = p.files_info
replicates_info = p.replicates_info
samples_info = p.samples_info
#%%
a = p.files['210_FWCP250C1h1_1']
aa = p.replicates['FWCP250C1h1_1']
aaa = p.replicates['FWCP250C1h1_2']
b = p.samples['FWCP250C1h1']
bb = p.samples_std['FWCP250C1h1']

p.create_all_statistics()
files_report = p.create_files_param_report()
replicates_report = p.create_replicates_param_report()
replicates_aggrrep = p.create_replicates_param_aggrrep()
c, d = p.create_samples_param_report()
p.save_files_report()
p.save_replicates_report()
p.save_samples_param_report()
p.save_samples_param_aggrrep()
zz, zzstd = p.return_param_aggrrep(param='f_sample')
# %%
p.plot_ave_std(param='conc_vial_mg_L', aggr=True, min_y_thresh=0,
    y_lim=[0, 50000],legend_location='outside',
    color_palette='Set2')
# %%
p.plot_ave_std(replicates_or_samples='replicates', param='conc_vial_mg_L',
    min_y_thresh=2000, legend_location='outside', xlab_rot=20,
    # only_samples_to_plot=['FW250C1h1', 'FWCP250C1h1'],
    y_lim=[0, 50000], annotate_outliers=True)
# %%
p.plot_ave_std(replicates_or_samples='replicates', param='conc_vial_mg_L',
    aggr=True, min_y_thresh=2000, legend_location='outside', xlab_rot=20,
    # only_samples_to_plot=['FW250C1h1', 'FWCP250C1h1'],
    y_lim=[0, 50000], annotate_outliers=True)
# %%
