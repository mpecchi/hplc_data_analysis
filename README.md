# HPLC data analysis in Python

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

![Python 3.10](https://img.shields.io/badge/python-3.10%2B-blue)

[![Testing (CI)](https://github.com/mpecchi/hplc_data_analysis/actions/workflows/continuous_integration.yaml/badge.svg)](https://github.com/mpecchi/hplc_data_analysis/actions/workflows/continuous_integration.yaml)

[![Publishing (CI)](https://github.com/mpecchi/hplc_data_analysis/actions/workflows/python-publish.yaml/badge.svg)](https://github.com/mpecchi/hplc_data_analysis/actions/workflows/python-publish.yaml)


The `hplc_data_analysis` tool automates the typical analysis of HPLC data, saving time, avoiding human error, and increasing comparability of results from different groups. 

Some key features:
- handle multiple HPLC semi-quantitative data tables (obtained with different methods)
- duild a database of all identified compounds and their relevant properties using PubChemPy
- split each compound into its functional groups using a published fragmentation algorithm
- produce single file report, single replicate (intended as the sum of more methods applied to one vial), comprehensive multi-sample (average and deviation of replicates) reports and aggregated reports based on functional group mass fractions in the samples
- provides plotting capabilities

## Framework

### File

A **.txt** or **.csv** file located in the project folder that contains **time**, **area**, and **concentration** information for many **compunds** for a measure.

Depending on the instrument export parameter, the structure of ``Files`` can differ slightly. ``Project-Parameters`` ensure that the loading process can lead to the same data structure to allow to perform all downstream computations.

A good naming convention for ``Files`` ensures the code handles replicates of the same sample correctly. Filenames have to follow the convention:
*method_name-of-sample-with-dashes-only_replicatenumber*

Examples that are *correctly* handled:
- 210_Bio-oil-foodwaste-250C_1
- 254_Bio-oil-foodwaste-250C_1
- 210_FW_2
- 254_FW_2

Examples of *NON-ACCEPTABLE* names are
- 210-bio_oil_1
- 254-FW1

### Replicate

If more ``Files`` belong to the same material (``Sample``, see below) but represent different methods that see different compounds (for example, different wavelengths are used in the detector), they can be merged into the same ``Replicate``.

A ``Replicate`` is the union of files with different methods that are complementary in the analysis of a material.


### Sample

A collection of ``Replicates`` that replicate the same measure and allow to assess  reproducibility.


### Project

The ``folder path`` indicates where the ``Files`` are located and where the ``output`` folder will be created.

The ``Project-Parameters`` are valid for each ``Sample``.


The ``Project`` can generate ``Reports`` and ``Plots`` for all ``Files``, ``Replicates``, or ``Sample`` or only for some of them.

### Reports 

Reports contain the results for one ``parameter`` (abbreviated as ``param``) for all ``Files``, ``Replicates``, or ``Sample``.

There are two types of reports:

#### Reports (simple-reports or compound-reports)
These report report the ``param`` value for each compound in each ``Files``, ``Replicates``, or ``Sample``.

Example: the values of ``conc_vial_mg_L`` for each compound in each ``File`` are collected in a single pandas dataframe (and saved as excel worksheet) for an easy comparison.

#### Aggrreps (aggregated reports)
These report report the ``param`` value for each aggregated functional group in each ``Files``, ``Replicates``, or ``Sample``.

The results of componds are **aggregated** by functional group (see this [paper](https://doi.org/10.1039/D3SU00345K) for details).

### Plots

Each report can be plotted using the ``plot_report`` method of the ``Project`` class.


## Documentation

Check out the [documentation](https://hplc-data-analysis.readthedocs.io/).

## Installation

You can install the package from [PyPI](https://pypi.org/project/hplc_data_analysis/):

## Examples

Each example is available as a folder in the ``examples`` folder and contains the code and the necessary input data.
To run examples:
1. Install ``hplc_data_analysis`` in your Python environment
2. Download the folder that contains the example
3. Run the code 
4. If you run the scripts as Jupyter Notebooks, replace the relative path at the beginning of each example with the absolute path to the folder where the code is located 

## Plotting with myfigure

Plots rely on the package ``myfigure``, a package to simplify scientific plotting in data analysis packages.
Check out its [documentation](https://myfigure.readthedocs.io/) and 
[GitHub](https://github.com/mpecchi/myfigure/).