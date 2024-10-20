[![DOI](https://zenodo.org/badge/616543497.svg)](https://zenodo.org/badge/latestdoi/616543497)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Welcome to PyParse!
===================================

Authors: [Joe Mason](https://github.com/thatchemistryguy), Francesco Rianjongdee, Harry Wilders, [David Fallon](https://github.com/fallonda)


Description
--------------- 

This script will read Liquid Chromatography Mass Spectrometry (LCMS) data in the Waters OpenLynx™ browser report (.rpt)
or Shimadzu .daml file formats, and assign peaks to compounds specified in  a .csv platemap. This assignment is then used to generate heatmaps and 
other visualisations to compare and contrast different LCMS runs. It was designed specifically for the analysis of data generated from 
high-throughput chemistry, and is suitable for reaction optimisations, parallel synthesis, library validation experiments and direct-to-biology. 


Example Usage 
---------------

Check out the [Github Pages](https://thatchemistryguy.github.io/PyParse/index.html) site! Here you'll find full documentation, including a walkthrough using the example data set provided!

You can also find our published user-guide for chemists at our peer-reviewed article at [Digital Discovery](https://doi.org/10.1039/D3DD00167A)

```

	python PyParse.py example_rpt.rpt example_platemap.csv -o new_output_directory

```
(Saves all output tables, data and visualisations to "new_output_directory".)

New!! About This Branch
-------------------

This branch, genAnalyticalTable, generates a .csv file containing all possible information about every analyte for every well in the plate. It is designed to encompass everything, in a machine-readable format ready
for upload to a database or similar. 

It is most useful when used in conjunction with [PyParse_designer](https://github.com/thatchemistryguy/PyParse_designer), which a simple and lightweight tool to help the user generate fully detailed platemaps
with minimal effort. When used, the user gets a full end-to-end workflow to design the plate, analyse the data, and prepare this data for long term storage in a database. 

Workflow in a Nutshell:

1. PyParse_designer is used to define every parameter for every well in the plate
2. This platemap is then combined with the reaction UPLC-MS data in PyParse, which looks to assign peaks to analytes.
3. PyParse generates an analytical table, which includes information about both the inputs (catalyst amount, ID, temperature, etc) and the outputs (product percentage area, retention time observed, m/z values, etc)
4. The user uploads this analytical table to a global data lake for long-term storage, for future interrogation, and to facilitate the creation of machine-learning models.

Citation
-----------

Publications which make use of PyParse to aid analysis of high-throughput LC-MS data should cite the peer-reviewed article:


Mason J., Wilders H., Fallon D.J., Thomas R.P., Bush J.T., Tomkinson N.C.O., Rianjongdee, F.; Automated LC-MS Analysis and Data Extraction for High-Throughput Chemistry; Digital Discovery (**2023**), 2, 1894 - 1899; https://doi.org/10.1039/D3DD00167A

(An earlier version of this manuscript was published on ChemRxiv (**2023**), https://doi.org/10.26434/chemrxiv-2023-1x288)
		
License
---------------

Apache 2.0



