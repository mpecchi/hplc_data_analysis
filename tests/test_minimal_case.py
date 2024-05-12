# %%
import pytest


folder_path = r"/Users/matteo/Projects/hplc_data_analysis/tests/data_minimal_case"


def test_values(hplcproject):
    files_info_created = hplcproject.create_files_info()
    # load the files info with dilution data and concentration
    files_info = hplcproject.load_files_info(update_saved_files_info=False)
    # create replicates and samples info
    replicates_info = hplcproject.create_replicates_info()
    samples_info = hplcproject.create_samples_info()

    # to get a list of all files, replicates and samples
    list_of_all_filenames = files_info.index.tolist()
    list_of_all_replicatenames = replicates_info.index.tolist()
    list_of_all_samplenames = samples_info.index.tolist()

    assert list_of_all_filenames == [
        "210_FW250C1h1_1",
        "210_FW250C1h1_2",
        "254_FW250C1h1_1",
        "254_FW250C1h1_2",
        "254_FWCPMn250C1h1_1",
    ]
    assert list_of_all_replicatenames == [
        "FW250C1h1_1",
        "FW250C1h1_2",
        "FWCPMn250C1h1_1",
    ]
    assert list_of_all_samplenames == [
        "FW250C1h1",
        "FWCPMn250C1h1",
    ]
    hplcproject.create_samples()
    params = [
        "area",
        "conc_vial_mg_L",
        "conc_vial_if_undiluted_mg_L",
        "fraction_of_sample_fr",
        "fraction_of_feedstock_fr",
    ]
    for param, results in zip(
        params,
        [
            [150.0, 200.5, 350.5],
            [10.0, 10.5, 30.5],
            [10.0, 16.0, 107.0],
            [0.1, 0.1, 0.24423076923076925],
            [0.05, 0.055, 0.1828846153846154],
        ],
    ):
        assert hplcproject.samples["FW250C1h1"].ave[param].tolist() == results


# %%
