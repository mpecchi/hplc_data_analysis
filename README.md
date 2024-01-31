This example demonstrates how to use the hplc_data_analysis package for basic High-Performance Liquid Chromatography (HPLC) data analysis and visualization. It covers the initial setup, data preparation, and plotting of average and standard deviation for a specific parameter. The script sets up the project folder, initializes a Project instance, and processes HPLC data files. It then creates a parameter report, assigns colors to compounds for visualization purposes, and finally plots the data. This example is ideal for getting started with analyzing HPLC data using our package.

```python
# Import necessary libraries
import pandas as pd
import seaborn as sns
import pathlib as plib

# Import the Project class from hplc_data_analysis package
from hplc_data_analysis import Project

# Set the folder path for the project
folder_path = plib.Path(r"C:\Path\To\Your\Data")
Project.set_folder_path(folder_path)

# Initialize a Project instance
p = Project(rebuild_compounds_properties=False)

# Generate files, replicates, and samples
p.create_files_replicates_samples()

# Generate statistical data
p.create_all_statistics()

# Create and retrieve parameter reports
param_report, _ = p.create_param_report()

# Define compounds for visualization
compounds = param_report.index.tolist()

# Map compounds to colors for visualization
compounds_to_colors = pd.DataFrame({
    'item': compounds,
    'clr': sns.color_palette('Set2', len(compounds))
}).set_index('item')

# Plotting average and standard deviation for a specific parameter
p.plot_ave_std(filename='example_plot', param='f_sample', aggr=False,
               min_y_thresh=0.02, y_lim=[0, 0.25],
               item_to_color_to_hatch=compounds_to_colors)