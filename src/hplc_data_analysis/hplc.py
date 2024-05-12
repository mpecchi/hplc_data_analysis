# %%
from __future__ import annotations
import pathlib as plib
import pandas as pd
from typing import Literal
import matplotlib.patches as mpatches
from matplotlib.axes import Axes
from hplc_data_analysis.pubchem import name_to_properties
from myfigure.myfigure import MyFigure, colors, hatches


class Project:

    def __init__(
        self,
        folder_path: plib.Path,
        name: str | None = None,
        auto_save_reports: bool = True,
        load_delimiter: str = "\t",
        load_skiprows: int = 18,
        files_info_defauls_columns: list[str] | None = None,
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
        self.acceptable_params: list[str] = list(self.param_to_axis_label.keys())
        self.files_info: pd.DataFrame | None = None
        self.replicates_info: pd.DataFrame | None = None
        self.samples_info: pd.DataFrame | None = None
        self.list_of_unique_compounds: list[str] | None = None
        self.class_code_frac: pd.DataFrame | None = None
        self.dict_classes_to_codes: dict[str, str] | None = None
        self.dict_classes_to_mass_fractions: dict[str, float] | None = None
        self.compounds_properties: pd.DataFrame | None = None

        self.samples: dict[str, Sample] = {}
        self.file_dfs: dict[str, pd.DataFrame] = {}
        self.replicate_dfs: dict[str, pd.DataFrame] = {}
        self.sample_dfs: dict[str, pd.DataFrame] = {}
        self.sample_dfs_std: dict[str, pd.DataFrame] = {}
        self.samplenames: list[str] = []
        self.files_reports: dict[str, pd.DataFrame] = {}
        self.replicates_reports: dict[str, pd.DataFrame] = {}
        self.samples_reports: dict[str, pd.DataFrame] = {}
        self.samples_reports_std: dict[str, pd.DataFrame] = {}
        self.files_aggrreps: dict[str, pd.DataFrame] = {}
        self.replicates_aggrreps: dict[str, pd.DataFrame] = {}
        self.samples_aggrreps: dict[str, pd.DataFrame] = {}
        self.samples_aggrreps_std: dict[str, pd.DataFrame] = {}

        self.list_of_files_param_reports = []
        self.list_of_replicates_param_reports = []
        self.list_of_samples_param_reports = []
        self.list_of_files_param_aggrreps = []
        self.list_of_replicates_param_aggrreps = []
        self.list_of_samples_param_aggrreps = []

    def load_files_info(self) -> pd.DataFrame:
        """ """
        files_info_path = plib.Path(self.folder_path, "files_info.xlsx")
        if files_info_path.exists():
            self.files_info = pd.read_excel(
                files_info_path, engine="openpyxl", index_col="filename"
            )
            print("Info: files_info loaded")
        else:
            print("Info: files_info not found")
            self.create_files_info()
        self._add_default_to_files_info()
        self.files_info.to_excel(plib.Path(self.folder_path, "files_info.xlsx"))
        return self.files_info

    def create_files_info(self) -> pd.DataFrame:
        """ """
        filename: int = [a.parts[-1].split(".")[0] for a in list(self.folder_path.glob("**/*.txt"))]
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
        return self.files_info

    def _add_default_to_files_info(self):
        """ """
        for col in self.files_info_defauls_columns:
            if col not in list(self.files_info):
                self.files_info[col] = 1

    def create_replicates_info(self):
        """Creates a summary 'replicates_info' DataFrame from 'files_info',
        aggregating data for each replicate, and updates the 'replicates_info'
        attribute with this summarized data."""
        if self.files_info is None:
            _ = self.load_files_info()
        self.replicates_info = self.files_info.reset_index().groupby("replicatename").agg(list)

        self.replicates_info.reset_index(inplace=True)
        self.replicates_info.set_index("replicatename", drop=True, inplace=True)
        print("Info: create_replicates_info: replicates_info created")
        return self.replicates_info

    def create_samples_info(self):
        """Creates a summary 'samples_info' DataFrame from 'files_info',
        aggregating data for each sample, and updates the 'samples_info'
        attribute with this summarized data."""
        if self.replicates_info is None:
            _ = self.create_replicates_info()
        self.samples_info = self.files_info.reset_index().groupby("samplename").agg(list)
        self.replicates_info["samplename"] = [a[0] for a in self.replicates_info["samplename"]]
        self.samples_info.reset_index(inplace=True)
        self.samples_info.set_index("samplename", drop=True, inplace=True)
        print("Info: create_samples_info: samples_info created")
        return self.samples_info

    def create_samples(self):
        if self.samples_info is None:
            self.create_samples_info()
        for samplename in self.samples_info.index.tolist():
            sample_info = self.files_info.loc[self.files_info["samplename"] == samplename, :]
            self.samples[samplename] = Sample(self, samplename, sample_info)

    def create_files_param_report(self, param="conc_vial_mg_L"):
        """
        Create a report that consolidates the values of a specified parameter from different DataFrames,
        using the union of all indices found in the individual DataFrames.

        :param param: The parameter to extract from each DataFrame. Defaults to "conc_vial_mg_L".
        :return: A DataFrame containing the consolidated report.
        """
        if not self.sample_dfs:
            self.create_samples()
        # Create a dictionary of Series, each Series named after the file and containing the 'param' values
        series_dict = {
            filename: self.file_dfs[filename][param].rename(filename)
            for filename in self.files_info.index
            if param in self.file_dfs[filename].columns
        }
        # Get the union of all indices from the individual DataFrames
        rep = pd.concat(series_dict.values(), axis=1, keys=series_dict.keys(), join="outer")
        # Reindex the DataFrame to include all unique indices, filling missing values with 0
        rep = rep.sort_index(key=rep.max(axis=1).get, ascending=False)
        rep = rep.loc[:, rep.any(axis=0)]
        # Save and return the report
        self.files_reports[param] = rep.fillna(0)
        self.list_of_files_param_reports.append(param)
        return self.files_reports[param]

    def create_replicates_param_report(self, param="conc_vial_mg_L"):
        """
        Create a report that consolidates the values of a specified parameter from different DataFrames,
        using the union of all indices found in the individual DataFrames.

        :param param: The parameter to extract from each DataFrame. Defaults to "conc_vial_mg_L".
        :return: A DataFrame containing the consolidated report.
        """
        if not self.sample_dfs:
            self.create_samples()
        # Create a dictionary of Series, each Series named after the replicate and containing the 'param' values
        series_dict = {
            replicatename: self.replicate_dfs[replicatename][param].rename(replicatename)
            for replicatename in self.replicates_info.index
            if param in self.replicate_dfs[replicatename].columns
        }
        # Get the union of all indices from the individual DataFrames
        rep = pd.concat(series_dict.values(), axis=1, keys=series_dict.keys(), join="outer")
        # Sort by the max value in each row, then filter out columns that only contain 0s
        rep = rep.sort_index(key=rep.max(axis=1).get, ascending=False)
        rep = rep.loc[:, rep.any(axis=0)]
        # Save and return the report
        self.replicates_reports[param] = rep.fillna(0)
        self.list_of_replicates_param_reports.append(param)
        return self.replicates_reports[param]

    def create_samples_param_report(self, param="conc_vial_mg_L"):
        """
        Create two reports that consolidate the average and standard deviation of a specified parameter
        from different sample DataFrames, assuming both sets of DataFrames share the same indices.

        :param param: The parameter to extract from each DataFrame. Defaults to "conc_vial_mg_L".
        :return: A tuple of two DataFrames containing the consolidated averages and standard deviations.
        """
        if not self.sample_dfs:
            self.create_samples()

        series_dict = {
            samplename: self.sample_dfs[samplename][param].rename(samplename)
            for samplename in self.samples_info.index
            if param in self.sample_dfs[samplename].columns
        }
        series_dict_std = {
            samplename: self.sample_dfs_std[samplename][param].rename(samplename)
            for samplename in self.samples_info.index
            if param in self.sample_dfs_std[samplename].columns
        }
        # Get the union of all indices from the individual sample DataFrames (assuming indices are the same for std and avg)
        rep = pd.concat(series_dict.values(), axis=1, keys=series_dict.keys(), join="outer")
        rep_std = pd.concat(
            series_dict_std.values(), axis=1, keys=series_dict_std.keys(), join="outer"
        )
        # Populate the DataFrames with values

        # Sort by the max value in each row and filter out columns that only contain 0s in the average report
        rep = rep.sort_index(key=rep.max(axis=1).get, ascending=False)
        rep = rep.loc[:, rep.any(axis=0)]
        # Ensure the standard deviation DataFrame aligns with the average DataFrame
        rep_std = rep_std.reindex_like(rep)

        # Save and return the reports
        self.samples_reports[param] = rep.fillna(0)
        self.samples_reports_std[param] = rep_std
        self.list_of_samples_param_reports.append(param)

        return self.samples_reports[param], self.samples_reports_std[param]

    def load_class_code_frac(self):
        """ """
        class_code_frac_path = plib.Path(self.folder_path, "classifications_codes_fractions.xlsx")
        if class_code_frac_path.exists():
            self.class_code_frac = pd.read_excel(class_code_frac_path)
        else:
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
            self.compounds_properties = pd.read_excel(
                compounds_properties_path, index_col="comp_name"
            )
            print("Info: compounds_properties loaded")
        else:
            print("Warning: compounds_properties.xlsx not found, creating it")
            self.compounds_properties = self.create_compounds_properties()

        return self.compounds_properties

    def create_compounds_properties(self):
        """ """
        if self.list_of_unique_compounds is None:
            self.create_list_of_unique_compounds()
        if self.class_code_frac is None:
            self.class_code_frac = self.load_class_code_frac()

        self.compounds_properties = pd.DataFrame()
        for name in self.list_of_unique_compounds:
            self.compounds_properties = name_to_properties(
                comp_name=name,
                dict_classes_to_codes=self.dict_classes_to_codes,
                dict_classes_to_mass_fractions=self.dict_classes_to_mass_fractions,
                df=self.compounds_properties,
            )
        self.compounds_properties.index.name = "comp_name"
        # save db in the project folder in the input
        self.compounds_properties.to_excel(plib.Path(self.folder_path, "compounds_properties.xlsx"))
        print("Info: create_compounds_properties: compounds_properties created and saved")
        return self.compounds_properties

    def create_list_of_unique_compounds(self) -> list[str]:
        if len(self.sample_dfs) == 0:
            self.create_samples()

        all_compounds = pd.concat([df for df in self.sample_dfs.values()])
        self.list_of_unique_compounds = pd.Index(all_compounds.index.unique())
        return self.list_of_unique_compounds

    def create_files_param_aggrrep(self, param="conc_vial_mg_L"):
        """Aggregates compound concentration data by functional group for each
        parameter across all FILES, providing a summarized view of functional
        group concentrations. This aggregation facilitates the understanding
        of functional group distribution across FILES."""
        print("Info: create_param_aggrrep: ", param)
        if param not in self.acceptable_params:
            raise ValueError(f"{param = } is not an acceptable param")
        if param not in self.list_of_files_param_reports:
            self.create_files_param_report(param)
        if self.compounds_properties is None:
            self.load_compounds_properties()
        # fg = functional groups, mf = mass fraction
        filenames = self.files_info.index.tolist()
        _all_comps = self.files_reports[param].index.tolist()
        cols_with_fg_mf_labs = list(self.compounds_properties)
        fg_mf_labs = [
            c for c in cols_with_fg_mf_labs if c.startswith("fg_mf_") if c != "fg_mf_total"
        ]
        fg_labs = [c[6:] for c in fg_mf_labs]
        # create a df with iupac name index and fg_mf columns (underiv and deriv)
        all_comps_df = self.compounds_properties
        all_comps_df = all_comps_df[~all_comps_df.index.duplicated(keep="first")]
        fg_mf_all = pd.DataFrame(index=_all_comps, columns=fg_mf_labs)
        for idx in fg_mf_all.index.tolist():
            fg_mf_all.loc[idx, fg_mf_labs] = all_comps_df.loc[idx, fg_mf_labs]
        # create the aggregated dataframes and compute aggregated results
        aggrrep = pd.DataFrame(columns=filenames, index=fg_labs, dtype="float")
        aggrrep.index.name = param  # is the parameter
        for col in filenames:
            list_iupac = self.files_reports[param].index
            signal = self.files_reports[param].loc[:, col].values
            for fg, fg_mf in zip(fg_labs, fg_mf_labs):
                # each compound contributes to the cumulative sum of each
                # functional group for the based on the mass fraction it has
                # of that functional group (fg_mf act as weights)
                # if fg_mf in subrep: multiply signal for weight and sum
                # to get aggregated
                weights = fg_mf_all.loc[list_iupac, fg_mf].astype(signal.dtype)

                aggrrep.loc[fg, col] = (signal * weights).sum()
        aggrrep = aggrrep.loc[(aggrrep != 0).any(axis=1), :]  # drop rows with only 0
        aggrrep = aggrrep.sort_index(key=aggrrep[filenames].max(1).get, ascending=False)
        self.files_aggrreps[param] = aggrrep
        self.list_of_files_param_aggrreps.append(param)
        return aggrrep

    def create_replicates_param_aggrrep(self, param="conc_vial_mg_L"):
        """Aggregates compound concentration data by functional group for each
        parameter across all FILES, providing a summarized view of functional
        group concentrations. This aggregation facilitates the understanding
        of functional group distribution across FILES."""
        print("Info: create_param_aggrrep: ", param)
        if param not in self.acceptable_params:
            raise ValueError(f"{param = } is not an acceptable param")
        if param not in self.list_of_replicates_param_reports:
            self.create_replicates_param_report(param)
        if self.compounds_properties is None:
            self.load_compounds_properties()
        # fg = functional groups, mf = mass fraction
        replicatenames = self.replicates_info.index.tolist()
        _all_comps = self.replicates_reports[param].index.tolist()
        cols_with_fg_mf_labs = list(self.compounds_properties)
        fg_mf_labs = [
            c for c in cols_with_fg_mf_labs if c.startswith("fg_mf_") if c != "fg_mf_total"
        ]
        fg_labs = [c[6:] for c in fg_mf_labs]
        # create a df with iupac name index and fg_mf columns (underiv and deriv)
        comps_df = self.compounds_properties
        all_comps_df = comps_df
        all_comps_df = all_comps_df[~all_comps_df.index.duplicated(keep="first")]
        fg_mf_all = pd.DataFrame(index=_all_comps, columns=fg_mf_labs)
        for idx in fg_mf_all.index.tolist():
            fg_mf_all.loc[idx, fg_mf_labs] = all_comps_df.loc[idx, fg_mf_labs]
        # create the aggregated dataframes and compute aggregated results
        aggrrep = pd.DataFrame(columns=replicatenames, index=fg_labs, dtype="float")
        aggrrep.index.name = param  # is the parameter
        aggrrep.fillna(0, inplace=True)
        for col in replicatenames:
            list_iupac = self.replicates_reports[param].index
            signal = self.replicates_reports[param].loc[:, col].values
            for fg, fg_mf in zip(fg_labs, fg_mf_labs):
                # each compound contributes to the cumulative sum of each
                # functional group for the based on the mass fraction it has
                # of that functional group (fg_mf act as weights)
                # if fg_mf in subrep: multiply signal for weight and sum
                # to get aggregated
                weights = fg_mf_all.loc[list_iupac, fg_mf].astype(signal.dtype)

                aggrrep.loc[fg, col] = (signal * weights).sum()
        aggrrep = aggrrep.loc[(aggrrep != 0).any(axis=1), :]  # drop rows with only 0
        aggrrep = aggrrep.sort_index(key=aggrrep[replicatenames].max(1).get, ascending=False)
        self.replicates_aggrreps[param] = aggrrep
        self.list_of_replicates_param_aggrreps.append(param)
        return aggrrep

    def create_samples_param_aggrrep(self, param: str = "conc_vial_mg_L"):
        print(f"Info: create_samples_param_aggrrep: {param = }")
        if param not in self.acceptable_params:
            raise ValueError(f"{param = } is not an acceptable param")
        if param not in self.list_of_replicates_param_aggrreps:
            self.create_replicates_param_aggrrep(param)
        replicate_to_sample_rename = dict(
            zip(self.replicates_info.index.tolist(), self.replicates_info["samplename"])
        )
        replicateagg = self.replicates_aggrreps[param].copy()
        replicateagg.rename(columns=replicate_to_sample_rename, inplace=True)
        self.samples_aggrreps[param] = replicateagg.T.groupby(by=replicateagg.columns).mean().T
        self.samples_aggrreps_std[param] = replicateagg.T.groupby(by=replicateagg.columns).std().T
        self.list_of_samples_param_aggrreps.append(param)
        return self.samples_aggrreps[param], self.samples_aggrreps_std[param]

    def plot_report(
        self,
        report_or_aggrrep: Literal["report", "aggrrep"] = "report",
        files_replicates_or_samples: Literal["files", "replicates", "samples"] = "samples",
        param: str = "conc_vial_mg_L",
        names_to_keep: list[str] | None = None,
        labels: list[str] | None = None,
        show_total_in_twinx: bool = False,
        y_axis_min_threshold: float | None = None,
        item_to_color_to_hatch: pd.DataFrame | None = None,
        yt_sum_label: str = "total\n(right axis)",
        remove_insignificant_values: bool = False,
        **kwargs,
    ) -> MyFigure:
        """ """
        if param not in self.acceptable_params:
            raise ValueError(f"{param = } is not an acceptable param")

        out_path = plib.Path(self.out_path, "plots")
        out_path.mkdir(parents=True, exist_ok=True)

        if report_or_aggrrep == "report":  # then use compounds reports
            if files_replicates_or_samples == "files":
                if param not in self.list_of_files_param_reports:
                    self.create_files_param_report(param)
                df_ave = self.files_reports[param].T
                df_std = pd.DataFrame()
            if files_replicates_or_samples == "replicates":
                if param not in self.list_of_replicates_param_reports:
                    self.create_replicates_param_report(param)
                df_ave = self.replicates_reports[param].T
                df_std = pd.DataFrame()
            elif files_replicates_or_samples == "samples":
                if param not in self.list_of_samples_param_reports:
                    self.create_samples_param_report(param)
                df_ave = self.samples_reports[param].T
                df_std = self.samples_reports_std[param].T
        else:  # use aggregated reports
            if files_replicates_or_samples == "files":
                if param not in self.list_of_files_param_aggrreps:
                    self.create_files_param_aggrrep(param)
                df_ave = self.files_aggrreps[param].T
                df_std = pd.DataFrame()
            if files_replicates_or_samples == "replicates":
                if param not in self.list_of_replicates_param_aggrreps:
                    self.create_replicates_param_aggrrep(param)
                df_ave = self.replicates_aggrreps[param].T
                df_std = pd.DataFrame()
            elif files_replicates_or_samples == "samples":
                if param not in self.list_of_samples_param_aggrreps:
                    self.create_samples_param_aggrrep(param)
                df_ave = self.samples_aggrreps[param].T
                df_std = self.samples_aggrreps_std[param].T

        if names_to_keep is not None:
            df_ave = df_ave.loc[names_to_keep, :].copy()
            if files_replicates_or_samples == "samples":
                df_std = df_std.loc[names_to_keep, :].copy()

        if labels is not None:
            df_ave.index = labels
            if files_replicates_or_samples == "samples":
                df_std.index = labels

        if y_axis_min_threshold is not None:
            df_ave = df_ave.loc[:, (df_ave > y_axis_min_threshold).any(axis=0)].copy()
            if files_replicates_or_samples == "samples":
                df_std = df_std.loc[:, df_ave.columns].copy()

        if item_to_color_to_hatch is not None:  # specific color and hatches to each fg
            plot_colors = [item_to_color_to_hatch.loc[item, "clr"] for item in df_ave.columns]
            plot_hatches = [item_to_color_to_hatch.loc[item, "htch"] for item in df_ave.columns]
        else:  # no specific colors and hatches specified
            plot_colors = colors
            plot_hatches = hatches

        if df_std.isna().all().all() or df_std.empty:
            std_available = False
        else:
            std_available = True

        if remove_insignificant_values:
            if std_available:
                mask = (df_ave.abs() > df_std.abs()) | df_std.isna()
                df_ave = df_ave[mask]
                df_std = df_std[mask]

        default_kwargs = {
            "filename": report_or_aggrrep + files_replicates_or_samples + param,
            "out_path": out_path,
            "height": 4,
            "width": 4,
            "grid": self.plot_grid,
            "text_font": self.plot_font,
            "y_lab": self.param_to_axis_label[param],
            "yt_lab": self.param_to_axis_label[param],
            "twinx": True if show_total_in_twinx else False,
            "auto_apply_hatches_to_bars": False,
        }
        # Update kwargs with the default key-value pairs if the key is not present in kwargs
        kwargs = {**default_kwargs, **kwargs}
        myfig = MyFigure(rows=1, cols=1, **kwargs)
        df_ave.plot(
            ax=myfig.axs[0],
            kind="bar",
            width=0.9,
            edgecolor="k",
            legend=False,
            capsize=3,
            color=plot_colors,
            yerr=df_std if std_available else None,
        )

        apply_hatches_to_ax(myfig.axs[0], plot_hatches)

        if show_total_in_twinx:
            myfig.axts[0].scatter(
                df_ave.index,
                df_ave.sum(axis=1).values,
                color="k",
                linestyle="None",
                edgecolor="k",
                facecolor="grey",
                s=100,
                label=yt_sum_label,
                alpha=0.5,
            )
            if std_available:
                myfig.axts[0].errorbar(
                    df_ave.index,
                    df_ave.sum(axis=1).values,
                    df_std.sum(axis=1).values,
                    capsize=3,
                    linestyle="None",
                    color="grey",
                    ecolor="k",
                    label="_nolegend_",
                )

        # Identify new patches added by the DataFrame plot
        myfig.save_figure()
        return myfig


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
                project.file_dfs[filename] = self.files[filename]
                _files.append(file)
            self.replicates[replicatename] = self.create_replicate_from_files(_files, replicatename)
            project.replicate_dfs[replicatename] = self.replicates[replicatename]

        ave, std = self.create_ave_std_from_replicates(list(self.replicates.values()))
        project.sample_dfs[self.samplename] = ave
        project.sample_dfs_std[self.samplename] = std

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

    def create_ave_std_from_replicates(self, replicates) -> tuple[pd.DataFrame]:
        # Align indices and columns of all DataFrames to the first replicate
        aligned_dfs = [df.align(replicates[0], join="outer", axis=0)[0] for df in replicates]
        aligned_dfs = [df.align(replicates[0], join="outer", axis=1)[0] for df in aligned_dfs]

        # Fill missing values with 0 in each DataFrame
        filled_dfs = [df.fillna(0) for df in aligned_dfs]

        # Calculate the average and standard deviation
        self.ave = pd.concat(filled_dfs).groupby(level=0).mean()
        self.std = pd.concat(filled_dfs).groupby(level=0).std()

        return self.ave, self.std


def apply_hatches_to_ax(ax: Axes, hatches_list: list[str]) -> None:
    """
    Apply hatch patterns to bars in the bar plots of each subplot.

    This method iterates over all subplots and applies predefined hatch patterns to each bar,
    enhancing the visual distinction between bars, especially in black and white printouts.
    """
    # Check if the plot is a bar plot
    bars = [b for b in ax.patches if isinstance(b, mpatches.Rectangle)]
    # If there are no bars, return immediately
    if not bars:
        return
    num_groups = len(ax.get_xticks(minor=False))
    # Determine the number of bars in each group
    bars_in_group = len(bars) // num_groups
    patterns = hatches_list[:bars_in_group]  # set hatch patterns in correct order
    plot_hatches_list = []  # list for hatches in the order of the bars
    for h in patterns:  # loop over patterns to create bar-ordered hatches
        for _ in range(int(len(bars) / len(patterns))):
            plot_hatches_list.append(h)
    # loop over bars and hatches to set hatches in correct order
    for b, hatch in zip(bars, plot_hatches_list):
        b.set_hatch(hatch)
        b.set_edgecolor("k")


# %%
# if __file__ == "main":
