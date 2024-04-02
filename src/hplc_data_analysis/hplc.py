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

        if columns_to_keep_in_files is None:
            self.columns_to_keep_in_files = ["R.Time", "Height", "Area", "Conc."]
        else:
            self.columns_to_keep_in_files = columns_to_keep_in_files

        if columns_to_rename_in_files is None:
            self.columns_to_rename_in_files = {
                "R.Time": "retention_time",
                "Height": "height",
                "Area": "area",
                "Conc.": "conc_vial_mg_L",
            }
        else:
            self.columns_to_rename_in_files = columns_to_rename_in_files

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

    def update_all_info_statistics(self):
        """ """
        if not self.files_replicates_samples_created:
            self.create_files_replicates_samples()
        for file in self.files.values():
            self._update_info_statistics(file, self.files_info)
        self.save_files_info()
        for replicate in self.replicates.values():
            self._update_info_statistics(replicate, self.replicates_info)
        self.save_replicates_info()
        for sample in self.samples.values():
            self._update_info_statistics(sample, self.samples_info)
        self.save_samples_info()

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
                plib.Path(Project.in_path, "classifications_codes_fractions.xlsx")
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
        """Retrieves and organizes properties for underivatized compounds using pubchempy,
        updating the 'compounds_properties' attribute and saving the properties
        to 'compounds_properties.xlsx'."""
        print("Info: create_compounds_properties: started")

        if not self.class_code_frac_loaded:
            self.load_class_code_frac()
        if not self.list_of_all_compounds_created:
            self.create_list_of_all_compounds()
        # cpdf = pd.DataFrame(index=pd.Index(self.list_of_all_compounds))
        #
        cpdf = pd.DataFrame()
        print("Info: create_compounds_properties: looping over names")
        for name in self.list_of_all_compounds:
            cpdf = name_to_properties(
                comp_name=name,
                dict_classes_to_codes=self.dict_classes_to_codes,
                dict_classes_to_mass_fractions=self.dict_classes_to_mass_fractions,
                df=cpdf,
            )
        # cpdf = self._order_columns_in_compounds_properties(cpdf)
        # cpdf = cpdf.fillna(0)
        cpdf.index.name = "comp_name"
        self.compounds_properties = cpdf
        self.compounds_properties_created = True
        # save db in the project folder in the input
        cpdf.to_excel(plib.Path(Project.in_path, "compounds_properties.xlsx"))
        print("Info: create_compounds_properties: compounds_properties created and saved")
        return self.compounds_properties


class Sample:

    def __init__(
        self,
        project: Project,
        name: str,
        filenames: list[str],
        folder_path: plib.Path | None = None,
        label: str | None = None,
        column_name_mapping: dict[str:str] | None = None,
        load_skiprows: int = 0,
        time_moist: float = 38.0,
        time_vm: float = 147,
        heating_rate_deg_min: float | None = None,
        temp_i_temp_b_threshold: float | None = None,
    ):
        # store the sample in the project
        self.project_name = project.name
        project.add_sample(name, self)
        # prject defaults unless specified

        self.out_path = project.out_path
        self.plot_font = project.plot_font
        self.plot_grid = project.plot_grid
        self.auto_save_reports = project.auto_save_reports

        if folder_path is None:
            self.folder_path = project.folder_path
        else:
            self.folder_path = folder_path
        if load_skiprows is None:
            self.load_skiprows = project.load_skiprows
        else:
            self.load_skiprows = load_skiprows
        # sample default
        self.name = name
        self.filenames = filenames
        self.n_repl = len(self.filenames)
        if not label:
            self.label = name
        else:
            self.label = label


# %%
# if __file__ == "main":
folder_path = plib.Path(r"C:\Users\mp933\OneDrive - Cornell University\Python\HPLC\SALLE")
hplc = Project(folder_path)
# hplc.create_files_info()
# %%
