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


