.. PyParse documentation master file, created by
   sphinx-quickstart on Sun Feb 19 15:26:24 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyParse's documentation!
===================================

Authors: Joseph Mason, Francesco Rianjongdee, Harry Wilders, David Fallon


Description
--------------- 

This script carries out an automated analysis of LCMS data in the Waters .rpt file 
format, using a .csv platemap to inform which compounds are found in which well of
a plate. It was designed specifically for the analysis of data generated from 
high-throughput chemistry, and is suitable for reaction optimisations, parallel synthesis
and library validation experiments. 

Example Usage 
---------------
		
python PyParse.py example_rpt.rpt example_platemap.csv -o new_output_directory -r 8 -c 12 -pt corrP/STD

(Used for an LCMS plate with 8 rows and 12 columns, to generate visualisations based 
on the corrected ratio of product to internal standard. Saves all output tables/data/visualisations
to "new_output_directory".)
		
		
License
---------------

Apache 2.0

Copyright 
---------------

2023 GlaxoSmithKline Research & Development Limited



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Contents
---------

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   installation
   output_desc
   add_param
   pyparse_det
