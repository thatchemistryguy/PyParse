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
high-throughput chemistry, and is suitable for reaction optimisations, parallel synthesis
and library validation experiments. 


Example Usage 
---------------

Check out the [Github Pages](https://thatchemistryguy.github.io/PyParse/index.html) site! Here you'll find full documentation, including a walkthrough using the example data set provided!

You can also find our published user-guide for chemists at our peer-reviewed article at [Digital Discovery](https://doi.org/10.1039/D3DD00167A)

```

	python PyParse.py example_rpt.rpt example_platemap.csv -o new_output_directory

```
(Saves all output tables, data and visualisations to "new_output_directory".)

Citation
-----------

Publications which make use of PyParse to aid analysis of high-throughput LC-MS data should cite the peer-reviewed article:


Mason J., Wilders H., Fallon D.J., Thomas R.P., Bush J.T., Tomkinson N.C.O., Rianjongdee, F.; Automated LC-MS Analysis and Data Extraction for High-Throughput Chemistry; Digital Discovery (**2023**), 2, 1894 - 1899; https://doi.org/10.1039/D3DD00167A

(An earlier version of this manuscript was published on ChemRxiv (**2023**), https://doi.org/10.26434/chemrxiv-2023-1x288)
		
License
---------------

Apache 2.0



