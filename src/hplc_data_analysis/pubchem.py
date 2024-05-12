import pathlib as plib
import numpy as np
import pandas as pd
import ele
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdmolops
from rdkit.Chem.AllChem import (  # pylint: disable=no-name-in-module
    GetMorganFingerprintAsBitVect,
)
from hplc_data_analysis.fragmenter import Fragmenter


def get_compound_from_pubchempy(comp_name: str) -> pcp.Compound:
    if not isinstance(comp_name, str) or comp_name.isspace():
        print(f"WARNING get_compound_from_pubchempy got an invalid {comp_name =}")
        return None
    cond = True
    while cond:  # to deal with HTML issues on server sides (timeouts)
        try:
            # comp contains all info about the chemical from pubchem
            try:
                comp_inside_list = pcp.get_compounds(comp_name, "name")
            except ValueError:
                print(f"{comp_name = }")
                return None
            if comp_inside_list:
                comp = comp_inside_list[0]
            else:
                print(
                    f"WARNING: name_to_properties {comp_name=} does not find an entry in pcp",
                )
                return None
            cond = False
        except pcp.PubChemHTTPError:  # timeout error, simply try again
            print("Caught: pcp.PubChemHTTPError (keep trying)")
    return comp


def _order_columns_in_compounds_properties(
    unsorted_df: pd.DataFrame | None,
) -> pd.DataFrame | None:
    if unsorted_df is None:
        return None
    priority_cols: list[str] = [
        "iupac_name",
        "underiv_comp_name",
        "molecular_formula",
        "canonical_smiles",
        "molecular_weight",
        "xlogp",
    ]

    # Define a custom sort key function
    def sort_key(col):
        if col in priority_cols:
            return (-1, priority_cols.index(col))
        if col.startswith("el_mf"):
            return (2, col)
        elif col.startswith("el_"):
            return (1, col)
        elif col.startswith("fg_mf_unclassified"):
            return (5, col)
        elif col.startswith("fg_mf"):
            return (4, col)
        elif col.startswith("fg_"):
            return (3, col)
        else:
            return (0, col)

    # Sort columns using the custom key
    sorted_columns = sorted(unsorted_df.columns, key=sort_key)
    sorted_df = unsorted_df.reindex(sorted_columns, axis=1)
    sorted_df.index.name = "comp_name"
    # Reindex the DataFrame with the sorted columns
    return sorted_df


def name_to_properties(
    comp_name: str,
    dict_classes_to_codes: dict[str:str],
    dict_classes_to_mass_fractions: dict[str:float],
    df: pd.DataFrame = pd.DataFrame(),
    precision_sum_elements: float = 0.05,
    precision_sum_functional_group: float = 0.05,
) -> pd.DataFrame:
    """
    used to retrieve chemical properties of the compound indicated by the
    comp_name and to store those properties in the df

    Parameters
    ----------
    GCname : str
        name from GC, used as a unique key.
    search_name : str
        name to be used to search on pubchem.
    df : pd.DataFrame
        that contains all searched compounds.
    df_class_code_frac : pd.DataFrame
        contains the list of functional group names, codes to be searched
        and the weight fraction of each one to automatically calculate the
        mass fraction of each compounds for each functional group.
        Classes are given as smarts and are looked into the smiles of the comp.

    Returns
    -------
    df : pd.DataFrame
        updated dataframe with the searched compound.
    CompNotFound : str
        if GCname did not yield anything CompNotFound=GCname.

    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("The argument df must be a pd.DataFrame.")

    if not isinstance(comp_name, str) or comp_name.isspace():
        return _order_columns_in_compounds_properties(df)

    if comp_name in df.index.tolist():
        return _order_columns_in_compounds_properties(df)

    comp = get_compound_from_pubchempy(comp_name)

    if comp is None:
        df.loc[comp_name, "iupac_name"] = "unidentified"
        return _order_columns_in_compounds_properties(df)

    try:
        valid_iupac_name = comp.iupac_name.lower()
    except AttributeError:  # iupac_name not give
        valid_iupac_name = comp_name.lower()

    df.loc[comp_name, "iupac_name"] = valid_iupac_name
    df.loc[comp_name, "molecular_formula"] = comp.molecular_formula
    df.loc[comp_name, "canonical_smiles"] = comp.canonical_smiles
    df.loc[comp_name, "molecular_weight"] = float(comp.molecular_weight)

    try:
        df.loc[comp_name, "xlogp"] = float(comp.xlogp)
    except TypeError:  # float() argument must be a string or a real number, not 'NoneType'
        df.loc[comp_name, "xlogp"] = np.nan
    elements = set(comp.to_dict()["elements"])
    el_dict = {}
    el_mf_dict = {}

    for el in elements:
        el_count = comp.to_dict()["elements"].count(el)
        el_mass = ele.element_from_symbol(el).mass

        # Using similar logic as in the fg_dict example
        if el not in el_dict:
            el_dict[el] = 0
            el_mf_dict[el] = 0.0

        el_dict[el] += int(el_count)
        el_mf_dict[el] += float(el_count) * float(el_mass) / float(comp.molecular_weight)
    # Now, update the DataFrame in a similar way to the fg_dict example
    for key, value in el_dict.items():
        df.at[comp_name, f"el_{key}"] = int(value)

    for key, value in el_mf_dict.items():
        df.at[comp_name, f"el_mf_{key}"] = float(value)
    cols_el_mf = [col for col in df.columns if col.startswith("el_mf_")]
    residual_els = df.loc[comp_name, cols_el_mf].sum() - 1
    # check element sum
    try:
        assert residual_els <= precision_sum_elements
    except AssertionError:
        print(f"the total mass fraction of elements in {comp_name =} is > 0.001")
    # apply fragmentation using the Fragmenter class (thanks simonmb)
    frg = Fragmenter(
        dict_classes_to_codes,
        fragmentation_scheme_order=dict_classes_to_codes.keys(),
        algorithm="simple",
    )
    fragmentation, _, _ = frg.fragment(comp.canonical_smiles)
    fg_dict = {}
    fg_mf_dict = {}
    # Iterate over each item in the dictionary
    for key, value in fragmentation.items():
        # Determine the root key (the part before an underscore, if present)
        root_key = key.split("_")[0]
        # if root_key in hetero_atoms:
        #     pass
        # Check if the root key is in the sum_dict; if not, initialize it
        if root_key not in fg_dict:
            fg_dict[root_key] = 0
            fg_mf_dict[root_key] = 0
        # Add the value to the corresponding root key in the sum_dict
        fg_dict[root_key] += int(fragmentation[key])
        fg_mf_dict[root_key] += (
            float(fragmentation[key])
            * float(dict_classes_to_mass_fractions[key])
            / df.loc[comp_name, "molecular_weight"].astype(float)
        )  # mass fraction of total

    # Update df with fg_dict
    for key, value in fg_dict.items():
        df.at[comp_name, f"fg_{key}"] = int(value)  # Update the cell
    # Update df with fg_mf_dict
    for key, value in fg_mf_dict.items():
        df.at[comp_name, f"fg_mf_{key}"] = float(value)  # Update the cell
    cols_fg_mf = [col for col in df.columns if col.startswith("fg_mf")]
    residual_fgs = df.loc[comp_name, cols_fg_mf].sum() - 1
    try:
        assert residual_fgs <= precision_sum_functional_group
    except AssertionError:
        print(f"{df.loc[comp_name, cols_fg_mf].sum()=}")
        print(f"the total mass fraction of functional groups in {comp_name =} is > 0.05")
    if residual_fgs < -precision_sum_functional_group:
        df.at[comp_name, "fg_mf_unclassified"] = abs(residual_fgs)
    df.loc[df["iupac_name"] != "unidentified"] = df.loc[df["iupac_name"] != "unidentified"].fillna(
        0
    )
    return _order_columns_in_compounds_properties(df)


# %%
def get_iupac_from_pcp(comp_name: str) -> str:
    """get iupac name for compound using pubchempy, needs internet connection

    :param comp_name: _description_
    :type comp_name: str
    :return: lowercase iupac name for the compound
    :rtype: str
    """
    cond = True
    while cond:  # to deal with HTML issues on server sides (timeouts)
        try:
            # comp contains all info about the chemical from pubchem
            try:
                comp = pcp.get_compounds(comp_name, "name")[0]

            except ValueError:
                print(f"Calibration iupac addition: compund {comp_name} did not work")
            try:
                iup: str = comp.iupac_name
            except AttributeError:  # iupac_name not give
                print(
                    f"Calibration iupac addition: compund {comp_name} does not have a iupac entry"
                )
            cond = False
        except pcp.PubChemHTTPError:  # timeout error, simply try again
            print("Caught: pcp.PubChemHTTPError")
    return iup.lower()


def report_difference(rep1, rep2, diff_type="absolute"):
    """
    calculates the ave, std and p percentage of the differnece between
    two reports where columns and index are the same.
    Replicates (indicated as XX_1, XX_2) are used for std.

    Parameters
    ----------
    rep1 : pd.DataFrame
        report that is conisdered the reference to compute differences from.
    rep2 : pd.DataFrame
        report with the data to compute the difference.
    diff_type : str, optional
        type of difference, absolute vs relative (to rep1)
        . The default is 'absolute'.

    Returns
    -------
    dif_ave : pd.DataFrame
        contains the average difference.
    dif_std : pd.DataFrame
        contains the std, same units as dif_ave.
    dif_stdp : pd.DataFrame
        contains the percentage std compared to ref1.

    """
    idx_name = rep1.index.name
    rep1 = rep1.transpose()
    rep2 = rep2.transpose()

    # put the exact same name on files (by removing the '_#' at end)
    repl_idx1 = [i if "_" not in i else i.split("_")[0] for i in rep1.index.tolist()]
    repl_idx2 = [i if "_" not in i else i.split("_")[0] for i in rep2.index.tolist()]
    rep1.loc[:, idx_name] = repl_idx1
    rep2.loc[:, idx_name] = repl_idx2
    # compute files and std of files and update the index
    rep_ave1 = rep1.groupby(idx_name, sort=False).mean().reset_index()
    rep_std1 = rep1.groupby(idx_name, sort=False).std().reset_index()
    rep_ave1.set_index(idx_name, inplace=True)
    rep_std1.set_index(idx_name, inplace=True)
    rep_ave2 = rep2.groupby(idx_name, sort=False).mean().reset_index()
    rep_std2 = rep2.groupby(idx_name, sort=False).std().reset_index()
    rep_ave2.set_index(idx_name, inplace=True)
    rep_std2.set_index(idx_name, inplace=True)

    if diff_type == "absolute":
        dif_ave = rep_ave1 - rep_ave2
        dif_std = np.sqrt(rep_std1**2 + rep_std2**2)
        dif_stdp = np.sqrt(rep_std1**2 + rep_std2**2) / dif_ave * 100
    if diff_type == "relative":
        dif_ave = (rep_ave1 - rep_ave2) / rep_ave1
        dif_std = np.sqrt(rep_std1**2 + rep_std2**2) / rep_ave1
        dif_stdp = np.sqrt(rep_std1**2 + rep_std2**2) / rep_ave1 / dif_ave * 100

    return dif_ave, dif_std, dif_stdp
