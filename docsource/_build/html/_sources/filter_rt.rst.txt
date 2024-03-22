

Filtering by Retention Time
==============================

The option "filter_by_rt" (-f) can be used to filter out regions of all chromatograms provided as an
input. Not only will PyParse remove any peaks in these regions from consideration, but it will also
recalculate all peak percentage areas to reflect the fact that this region was removed. 

Usage 
-------

To use this option, a range (or multiple ranges) should be provided in the form:
"[first_retention_time]-[second_retention_time]"

For example, to filter out all peaks between 0 and 0.4 minutes, you would add the following argument:

.. code-block::
	:caption: Filter out a specific range of the chromatogram
	
	-f 0-0.4

.. figure:: images/filter_single_range.jpg
    :alt: An example chromatogram where a specific range was filtered out

    An example chromatogram where a specific range was filtered out.

To filter out multiple ranges, input each range separated by a space:

.. code-block::
	:caption: Filter out a multiple ranges of the chromatogram
	
	-f 0-0.4 1.0-1.1

.. figure:: images/filter_multiple_range.jpg
    :alt: An example chromatogram where a multiple ranges were filtered out


The result of this filtering is that the purity reported by PyParse is higher:

.. figure:: images/filter_oldpurity.jpg
    :alt: Original purity value reported for Product1


.. figure:: images/filter_newpurity.jpg
    :alt: New purity value reported for Product1

    Top) Original purity value for Product1; Bottom) New purity value for Product1

The purpose of this feature is to enable the user to remove solvent peaks (e.g. DMSO, toluene)
from consideration. However, it should be noted that this feature should only be used with 
extreme caution to avoid misrepresenting the results!
