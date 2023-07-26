
Welcome to PyParse!
===================================

Authors: [Joe Mason](https://github.com/thatchemistryguy), Francesco Rianjongdee, Harry Wilders, [David Fallon](https://github.com/fallonda)


Description
--------------- 

This script will read Liquid Chromatography Mass Spectrometry (LCMS) data in the Waters OpenLynx™ browser report (.rpt) file
format, and assign peaks to compounds specified in  a .csv platemap. This assignment is then used to generate heatmaps and 
other visualisations to compare and contrast different LCMS runs. It was designed specifically for the analysis of data generated from 
high-throughput chemistry, and is suitable for reaction optimisations, parallel synthesis
and library validation experiments. 

See the [Waters OpenLynx™](https://www.waters.com/nextgen/ie/en/library/application-notes/2007/openlynx-open-access-and-software-tools-for-managing-an-open-access-laboratory-environment.html) page for further information about the OpenLynx™ file format. 

Example Usage 
---------------

Check out the [Github Pages](https://thatchemistryguy.github.io/PyParse/index.html) site! Here you'll find full documentation, including a walkthrough using the example data set provided!

You can also find our published user-guide for chemists at [ChemRxiv](https://doi.org/10.26434/chemrxiv-2023-1x288).

```

	python PyParse.py example_rpt.rpt example_platemap.csv -o new_output_directory -r 8 -c 12 -pt corrP/STD

```
(Used for an LCMS plate with 8 rows and 12 columns, to generate visualisations based 
on the corrected ratio of product to internal standard. Saves all output tables, data and visualisations
to "new_output_directory".)

Citation
-----------

Publications which make use of PyParse to aid analysis of high-throughput LC-MS data should cite the above ChemRxiv article:
Mason J., Wilders H., Fallon D.J., Thomas R.P., Bush J.T., Tomkinson N.C.O., Rianjongdee, F.; Automated LC-MS Analysis and Data Extraction for High-Throughput Chemistry; ChemRxiv (**2023**), DOI: 10.26434/chemrxiv-2023-1x288
		
License
---------------

Apache 2.0

Copyright 
---------------

2023 GlaxoSmithKline Research & Development Limited


