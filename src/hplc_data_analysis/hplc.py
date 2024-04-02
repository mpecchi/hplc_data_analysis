# %%
from __future__ import annotations
import pathlib as plib
import numpy as np
import pandas as pd
from typing import Any
from scipy.signal import savgol_filter
from lmfit.models import GaussianModel, LinearModel
from typing import Literal
from hplc_data_analysis.pubchem import (
    name_to_properties,
    get_compound_from_pubchempy,
    report_difference,
)
from hplc_data_analysis.myfigure import MyFigure, clrs, lnstls, htchs, mrkrs


class Project:

    def __init__(
        self,
        folder_path: plib.Path,
        name: str | None = None,
        auto_save_reports: bool = True,
        load_delimiter: str = "\t",
        load_skiprows: int = 18,
        files_info_defauls_columns: list[str] | None = None,
        columns_to_keep_in_files: list[str] | None = None,
        columns_to_rename_in_files: dict[str, str] | None = None,
        compounds_to_rename: dict[str, str] | None = None,
        param_to_axis_label: dict[str, str] | None = None,
        plot_font: Literal["Dejavu Sans", "Times New Roman"] = "Dejavu Sans",
        plot_grid: bool = False,
    ):
        self.folder_path = folder_path
        self.out_path = plib.Path(folder_path, "output")
        if name is None:
            self.name = self.folder_path.parts[-1]
        else:
            self.name = name
        self.plot_font = plot_font
        self.plot_grid = plot_grid
        self.auto_save_reports = auto_save_reports
        self.load_delimiter = load_delimiter
        self.load_skiprows = load_skiprows
        if files_info_defauls_columns is None:
            self.files_info_defauls_columns = [
                "dilution_factor",
                "total_sample_conc_in_vial_mg_L",
                "sample_yield_on_feedstock_basis_fr",
            ]
        else:
            self.files_info_defauls_columns = files_info_defauls_columns

        if columns_to_rename_in_files is None:
            self.columns_to_rename_in_files = {
                "Name": "comp_name",
                "R.Time": "retention_time",
                "Height": "height",
                "Area": "area",
                "Conc.": "conc_vial_mg_L",
            }
        else:
            self.columns_to_rename_in_files = columns_to_rename_in_files

        self.columns_to_keep_in_files = self.columns_to_rename_in_files.values()

        if compounds_to_rename is None:
            self.compounds_to_rename = {
                "3-methyl-(2H)-furan-5-one": "4-methyl-2H-furan-5-one",
                "4-methyl-(2H)-furan-5-one": "4-methyl-2H-furan-5-one",
                "2,3-pentanedione": "pentane-2,3-dione",
                "(2R,3S,4R,5R)-2,3,4,5,6-pentahydroxyhexanoic acid": "gluconic acid",
                "5-(hydroxymethyl)furan-2-carbaldehyde": "5-HMF",
            }
        else:
            self.compounds_to_rename = compounds_to_rename
        if param_to_axis_label is None:
            self.param_to_axis_label = {
                "AdjArea": "Peak Area [-]",
                "conc_vial_mg_L": "vial conc. [mg/L] (ppm)",
                "conc_vial_if_undiluted_mg_L": "vial conc. [mg/L] (ppm)",
                "fraction_of_sample_fr": "mass fraction [g/g$_{sample}$]",
                "fraction_of_feedstock_fr": "mass fraction [g/g$_{feedstock}$]",
            }
        else:
            self.param_to_axis_label = param_to_axis_label
        self.files_info = None
        self.replicates_info = None
        self.samples_info = None

        self.files_info_created = False
        self.replicates_info_created = False
        self.samples_info_created = False
        self.files = {}
        self.replicates = {}
        self.samples = {}
        self.samples_std = {}
        self.files_reports = {}
        self.replicates_reports = {}
        self.samples_reports = {}
        self.samples_reports_std = {}
        self.files_aggrreps = {}
        self.replicates_aggrreps = {}
        self.samples_aggrreps = {}
        self.samples_aggrreps_std = {}

        self.files_replicates_samples_created = False
        self.list_of_files_param_reports = []
        self.list_of_replicates_param_reports = []
        self.list_of_samples_param_reports = []
        self.list_of_files_param_aggrreps = []
        self.list_of_replicates_param_aggrreps = []
        self.list_of_samples_param_aggrreps = []

        self.samples: dict[str, Sample] = {}
        self.samplenames: list[str] = []

        self.multireports: dict[str, pd.DataFrame] = {}
        self.multireport_types_computed: list[str] = []

    def load_files_info(self):
        """Attempts to load the 'files_info.xlsx' file containing metadata about GCMS
        files. If the file is not found, it creates a new 'files_info' DataFrame with
        default values based on the GCMS files present in the project's input path and
        saves it to 'files_info.xlsx'. This method ensures 'files_info' is loaded with
        necessary defaults and updates the class attribute 'files_info_created' to True."""
        try:
            self.files_info = pd.read_excel(
                plib.Path(self.folder_path, "files_info.xlsx"),
                engine="openpyxl",
                index_col="filename",
            )
            self._add_default_to_files_info()
            print("Info: files_info loaded")
        except FileNotFoundError:
            print("Info: files_info not found")
            self.files_info = pd.DataFrame()
            self.create_files_info()
            self._add_default_to_files_info()
        self.files_info.to_excel(plib.Path(self.folder_path, "files_info.xlsx"))
        self.files_info_created = True
        return self.files_info

    def create_files_info(self):
        """ """
        filename = [a.parts[-1].split(".")[0] for a in list(self.folder_path.glob("**/*.txt"))]
        hplc_method = [f.split("_")[0] for f in filename]
        samplename = [f.split("_")[1] for f in filename]
        replicatenumber = [f.split("_")[2] for f in filename]
        replicatename = [s + "_" + r for s, r in zip(samplename, replicatenumber)]
        self.files_info = pd.DataFrame(
            {
                "filename": filename,
                "hplc_method": hplc_method,
                "samplename": samplename,
                "replicatename": replicatename,
            }
        )
        self.files_info.set_index("filename", drop=True, inplace=True)

    def _add_default_to_files_info(self):
        """ """
        for col in self.files_info_defauls_columns:
            if col not in list(self.files_info):
                self.files_info[col] = 1

    def add_sample(self, samplename: str, sample: Sample):
        """
        Add a sample to the project.

        :param samplename: The name of the sample to add.
        :type samplename: str
        :param sample: The sample object to add.
        :type sample: Sample
        """
        if samplename not in self.samplenames:
            self.samplenames.append(samplename)
            self.samples[samplename] = sample

        else:
            print(f"{samplename = } already present in project. Sample not added.")

    def create_replicates_info(self):
        """Creates a summary 'replicates_info' DataFrame from 'files_info',
        aggregating data for each replicate, and updates the 'replicates_info'
        attribute with this summarized data."""
        if not self.files_info_created:
            self.load_files_info()
        _replicates_info = self.files_info.reset_index().groupby("replicatename").agg(list)
        _replicates_info["samplename"] = [sn[0] for sn in _replicates_info["samplename"]]
        _replicates_info.reset_index(inplace=True)
        _replicates_info.set_index("replicatename", drop=True, inplace=True)
        self.replicates_info = _replicates_info
        self.replicates_info_created = True
        print("Info: create_replicates_info: replicates_info created")
        return self.replicates_info

    def create_samples_info(self):
        """Creates a summary 'samples_info' DataFrame from 'files_info',
        aggregating data for each sample, and updates the 'samples_info'
        attribute with this summarized data."""
        if not self.replicates_info_created:
            self.create_replicates_info()
        _samples_info = self.files_info.reset_index().groupby("samplename").agg(list)
        _samples_info.reset_index(inplace=True)
        _samples_info.set_index("samplename", drop=True, inplace=True)
        self.samples_info = _samples_info
        self.samples_info_created = True
        print("Info: create_samples_info: samples_info created")
        return self.samples_info

    def load_samples(self):
        if not self.samples_info_created:
            self.create_samples_info()
        for samplename in self.samples_info.index.tolist():
            sample_info = files_info.loc[files_info["samplename"] == samplename, :]
            self.samples[samplename] = Sample(self, samplename, sample_info)
        self.samples_created = True

    # def update_all_info_statistics(self):
    #     """ """
    #     if not self.files_replicates_samples_created:
    #         self.create_files_replicates_samples()
    #     for file in self.files.values():
    #         self._update_info_statistics(file, self.files_info)
    #     self.save_files_info()
    #     for replicate in self.replicates.values():
    #         self._update_info_statistics(replicate, self.replicates_info)
    #     self.save_replicates_info()
    #     for sample in self.samples.values():
    #         self._update_info_statistics(sample, self.samples_info)
    #     self.save_samples_info()

    def _update_info_statistics(self, df, info):
        """ """
        name = df.index.name
        # max values
        info.loc[name, "max_height"] = df["height"].max()
        info.loc[name, "max_area"] = df["area"].max()
        info.loc[name, "max_conc_vial_mg_L"] = df["conc_vial_mg_L"].max()
        info.loc[name, "max_conc_vial_if_undiluted_mg_L"] = df["conc_vial_if_undiluted_mg_L"].max()
        info.loc[name, "max_fraction_of_sample_fr"] = df["fraction_of_sample_fr"].max()
        info.loc[name, "max_fraction_of_feedstock_fr"] = df["fraction_of_feedstock_fr"].max()
        # total values
        info.loc[name, "total_conc_vial_mg_L"] = df["conc_vial_mg_L"].sum()
        info.loc[name, "total_conc_vial_if_undiluted_mg_L"] = df[
            "conc_vial_if_undiluted_mg_L"
        ].sum()
        info.loc[name, "total_fraction_of_sample_fr"] = df["fraction_of_sample_fr"].sum()
        info.loc[name, "total_fraction_of_feedstock_fr"] = df["fraction_of_feedstock_fr"].sum()
        info.loc[name, "compound_with_max_conc"] = df[
            df["conc_vial_mg_L"] == df["conc_vial_mg_L"].max()
        ].index[0]

    def load_class_code_frac(self):
        """Loads the 'classifications_codes_fractions.xlsx' file containing information
        on SMARTS classifications. It first searches in the project's input path, then
        in the shared path. It logs the status and returns the DataFrame containing
        classification codes and fractions."""
        try:  # first try to find the file in the folder
            self.class_code_frac = pd.read_excel(
                plib.Path(self.folder_path, "classifications_codes_fractions.xlsx")
            )
            print("Info: load_class_code_frac: classifications_codes_fractions loaded")

        except FileNotFoundError:
            raise FileNotFoundError('the file "classifications_codes_fractions.xlsx" was not found')
        all_classes = self.class_code_frac.classes.tolist()
        codes = self.class_code_frac.codes.tolist()  # list of code for each class
        mfs = self.class_code_frac.mfs.tolist()  # list of mass fraction of each class
        self.dict_classes_to_codes = dict(zip(all_classes, codes))  # dictionaries
        self.dict_classes_to_mass_fractions = dict(zip(all_classes, mfs))  # dictionaries
        return self.class_code_frac

    def load_compounds_properties(self):
        """Attempts to load the 'compounds_properties.xlsx' file containing physical
        and chemical properties of compounds. If not found, it creates a new properties
        DataFrame and updates the 'compounds_properties_created' attribute."""
        compounds_properties_path = plib.Path(self.folder_path, "compounds_properties.xlsx")
        if compounds_properties_path.exists():
            cpdf = pd.read_excel(
                compounds_properties_path,
                index_col="comp_name",
            )
            # cpdf = _order_columns_in_compounds_properties(cpdf)
            # cpdf = cpdf.fillna(0)
            self.compounds_properties = cpdf
            self.compounds_properties_created = True
            print("Info: compounds_properties loaded")
        else:
            print("Warning: compounds_properties.xlsx not found, creating it")
            cpdf = self.create_compounds_properties()
        return self.compounds_properties

    def create_compounds_properties(self):
        """ """
        if not self.files_replicates_samples_created:
            self.create_files_replicates_samples()
        self.class_code_frac = self.load_class_code_frac()

        all_compounds = pd.concat([df for df in self.samples.values()])
        unique_compounds = pd.Index(all_compounds.index.unique())
        cpdf = pd.DataFrame()
        for name in unique_compounds:
            cpdf = name_to_properties(
                comp_name=name,
                dict_classes_to_codes=self.dict_classes_to_codes,
                dict_classes_to_mass_fractions=self.dict_classes_to_mass_fractions,
                df=cpdf,
            )
        cpdf.index.name = "comp_name"
        self.compounds_properties = cpdf
        self.compounds_properties_created = True
        # save db in the project folder in the input
        cpdf.to_excel(plib.Path(self.folder_path, "compounds_properties.xlsx"))
        print("Info: create_compounds_properties: compounds_properties created and saved")
        return self.compounds_properties

    def create_files_param_report(self, param="conc_vial_mg_L"):
        """ """
        if not self.files_replicates_samples_created:
            self.load_samples()
        rep = pd.DataFrame(
            index=self.compounds_properties.index, columns=self.files_info.index, dtype="float"
        )
        rep.index.name = param

        for comp in rep.index.tolist():  # add conc values
            for filename in list(rep):
                try:
                    rep.loc[comp, filename] = self.files[filename].loc[comp, param]
                except KeyError:
                    rep.loc[comp, filename] = 0

        rep = rep.sort_index(key=rep.max(1).get, ascending=False)
        rep = rep.loc[:, rep.any(axis=0)]  # drop columns with only 0s
        self.files_reports[param] = rep
        self.list_of_files_param_reports.append(param)
        return self.files_reports[param]


class Sample:

    def __init__(
        self,
        project: Project,
        samplename: str,
        sample_info: pd.DataFrame,
    ):
        # store the sample in the project
        self.project_name = project.name

        # prject defaults unless specified

        self.load_skiprows = project.load_skiprows
        self.load_delimiter = project.load_delimiter
        self.files_info_defauls_columns = project.files_info_defauls_columns
        self.columns_to_keep_in_files = project.columns_to_keep_in_files
        self.columns_to_rename_in_files = project.columns_to_rename_in_files
        self.compounds_to_rename = project.compounds_to_rename
        self.auto_save_reports = project.auto_save_reports

        self.folder_path = project.folder_path

        self.samplename = samplename

        self.sample_info = sample_info
        self.files: dict[str, pd.DataFrame] = {}
        self.replicates: dict[str, pd.DataFrame] = {}
        self.replicate_files: dict[str, pd.DataFrame] = {}
        for replicatename in self.sample_info["replicatename"].tolist():
            replicate_info = self.sample_info.loc[
                self.sample_info["replicatename"] == replicatename, :
            ]
            _files = []
            for filename in replicate_info.index.tolist():
                file = self.load_single_file(filename)
                self.files[filename] = file
                _files.append(file)
            self.replicates[replicatename] = self.create_replicate_from_files(_files, replicatename)

        self.create_ave_std_from_replicates(list(self.replicates.values()))

    def load_single_file(self, filename: str) -> pd.DataFrame:

        file: pd.DataFrame = pd.read_csv(
            plib.Path(self.folder_path, filename + ".txt"),
            delimiter=self.load_delimiter,
            index_col=0,
            skiprows=self.load_skiprows,
        )
        file.rename(self.columns_to_rename_in_files, inplace=True, axis="columns")
        file = file.loc[file["comp_name"].notna(), self.columns_to_keep_in_files]
        file.set_index("comp_name", inplace=True)
        file.rename(self.compounds_to_rename, inplace=True)
        if any(file.index.duplicated(keep="first")):
            duplicates = file[file.index.duplicated(keep=False)]
            file = file[~file.index.duplicated(keep="first")]
            print(f"WARNING: duplicates in {filename = }")
            print(f"{duplicates = }, first instance has been kept")
        file = file.loc[file["conc_vial_mg_L"] > 0, :]
        file["conc_vial_if_undiluted_mg_L"] = (
            file["conc_vial_mg_L"] * self.sample_info.loc[filename, "dilution_factor"]
        )
        file["fraction_of_sample_fr"] = (
            file["conc_vial_mg_L"]
            / self.sample_info.loc[filename, "total_sample_conc_in_vial_mg_L"]
        )
        file["fraction_of_feedstock_fr"] = (
            file["fraction_of_sample_fr"]
            * self.sample_info.loc[filename, "sample_yield_on_feedstock_basis_fr"]
        )
        file.index.name = filename
        return file

    def create_replicate_from_files(self, files_to_merge, replicatename):
        """ """
        replicate = pd.concat(files_to_merge, join="outer")
        replicate = replicate.groupby(replicate.index).max()
        replicate.index.name = replicatename
        return replicate

    def create_ave_std_from_replicates(self, replicates) -> None:
        # Align indices and columns of all DataFrames to the first replicate
        aligned_dfs = [df.align(replicates[0], join="outer", axis=0)[0] for df in replicates]
        aligned_dfs = [df.align(replicates[0], join="outer", axis=1)[0] for df in aligned_dfs]

        # Fill missing values with 0 in each DataFrame
        filled_dfs = [df.fillna(0) for df in aligned_dfs]

        # Calculate the average and standard deviation
        self.ave = pd.concat(filled_dfs).groupby(level=0).mean()
        self.std = pd.concat(filled_dfs).groupby(level=0).std()

        return


# %%
# if __file__ == "main":
folder_path = plib.Path(r"C:\Users\mp933\OneDrive - Cornell University\Python\HPLC\SALLE")
hplc = Project(folder_path)
files_info = hplc.load_files_info()
# replicates_info = hplc.create_replicates_info()
samples_info = hplc.create_samples_info()
hplc.load_samples()
# %%
# hplc.create_compounds_properties()
# %%


# _samples_info = self.files_info.reset_index().groupby("samplename").agg(list)
# _samples_info.reset_index(inplace=True)
# _samples_info.set_index("samplename", drop=True, inplace=True)
# %%
# df = samples_info.loc[
#     a, :
# ]  # Convert the DataFrame such that each element of the list becomes a row
# # %%
# df_exploded = df.apply(pd.Series.explode).reset_index()
# # %%
# # Group by the original index to get individual dataframes
# grouped = df_exploded.groupby(df_exploded.index)

# # Create a dictionary or list to hold the individual dataframes
# dfs = {k: v.drop("index", axis=1) for k, v in grouped}
# # %%

# %%
