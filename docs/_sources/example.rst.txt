Using the Example Dataset
===========================

Introduction
--------------

To help you get started, an example dataset is provided in the repository 
under the folder "example_dataset". 

Included are:
	* An "example_rpt" file, which includes the LCMS traces for 24 compounds published in the public domain.
	* An "example_platemap" file, which is the corresponding platemap to enable PyParse analysis of the above data.
	* The output you should expect to see
	
Using the conda environment provided, you can use the sample data to run your first PyParse analysis!

.. code-block::
	:caption: Using the sample dataset provided
	
	python PyParse.py example_dataset/example_rpt.rpt example_dataset/example_platemap.csv -o example_dataset/output -r 2 -c 12 -pt Parea -moa area
	
The code provided uses the following options:

	* "-r 2" specifies that the data is split over two rows on the analysis plate
	* "-c 12" specifies that there are 12 columns in the analysis plate
	* "-pt Parea" specifies that the heatmap should be shaded by the percentage area of the product
	* "-moa area" specifies that where there are two peaks that could correspond to the product, choose the one with the larger percentage peak area. 
	
	
Expected Output
-------------------

Once the analysis is complete, navigate to the folder, where you should find a new sub-folder called output. 
Open this folder, and double-click on the "html_output" file. The HTML report will open, and contain the following
visualisations. 

You can also find these visualisations, and the underlying processed data, in the same subfolder. 

Heatmap
----------

.. figure:: images/example_heatmap.jpg
	:alt: An example heatmap.
	
	Example Heatmap
	
	
Check to make sure that your heatmap looks the same as the one shown above. If so, you are ready to begin your own 
analyses!!
