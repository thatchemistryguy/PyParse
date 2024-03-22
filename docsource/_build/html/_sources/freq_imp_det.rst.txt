.. _freq_impurity_detection_label:

Automatic Impurity Detection and Reporting
============================================

Background
------------------------

For a plate of reactions, there will typically be at least one reagent which is commonly used
across multiple wells of the plate. With this commonality is a higher chance that a particular
by-product or side-product (generically termed "impurity" for simplicity) is synthesised across multiple wells. 

Compared to single reactions that compare One Variable At a Time (OVAT), high-throughput experimentation
generates large quantities of data in a single experiment. What was previously a one-off impurity can become 
an observation backed up with data from multiple reactions. 

PyParse, through the hit-validation algorithm, overlays all peaks that weren't assigned to a compound specified 
by the user and identifies those that frequently occur. 

The Hit Validation Graph 
--------------------------

Below is a hit validation graph for the impurities in a 384-well plate containing reactions on random members of a diverse 
fragment library. The graph does not contain any real common impurities - it is a random assortment of peaks. 

.. figure:: images/384_hitvalidationgraph.jpg
    :alt: Hit Validation graph for 384-well plate

    Hit Validation Graph for 384-Well Plate

Here is the hit validation graph, generated using the same algorithm, for the 48-Well reaction optimisation plate 
included in the GitHub repository. 

.. note:: 
    For further details on the reaction optimisation plate, see: Mason J., Wilders H., Fallon D.J., Thomas R.P., 
    Bush J.T., Tomkinson N.C.O., Rianjongdee, F.; Automated LC-MS Analysis and Data Extraction for High-Throughput Chemistry; 
    Digital Discovery (2023), 2, 1894 - 1899; https://doi.org/10.1039/D3DD00167A

.. figure:: images/48_hitvalidationgraph.jpg
    :alt: Hit Validation graph for reaction optimisation plate

    Hit Validation Graph for Reaction Optimisation Plate

In this graph, the impurities line up at specific retention times; this consistency 
across multiple wells suggest that the impurity observed in each well may be the same.
For example, it can be seen that impurity6 (at 0.79 min) 
For these impurities, it is relevant to investigate the m/z values that are common to 
all peaks at each retention time. 




