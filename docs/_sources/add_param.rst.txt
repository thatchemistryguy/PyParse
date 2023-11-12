
Additional Parameters
========================

Introduction
--------------

There are optional parameters in addition to those detailed in :ref:`running_an_analysis_label` which are more rarely 
re-configured. 

Generic Parameters
-------------------

* -verbose (-v): "True" or "False"
	Specify whether logging should be run in debug mode or not

* -generate_csv (-g): "True" or "False"
	Choose whether to generate and save a csv of the outputTable
	
* -generate_zip (-z): "True" or "False"
	Choose whether to generate and save a zip file containing all visualisations and outputs


Parameters Used for the Analysis
----------------------------------

* -validate (-V): "True" or "False"
	Specify whether or not to run the hit validation algorithm. 
* -mass_abs_tol (-mat): Float value 
	How close an observed m/z should be to a calculate mass to be considered a match
* -time_abs_tol (-tat): Float value
	How close two retention times should be to be considered part of the same cluster/relate to the same compound
* -uv_abs_tol (-uat): Integer value
	How close should the wavelengths of the maxima in a UV absorption profile be for two eluting peaks to be considered a match
* -min_peak_area (-mpa): Float value
	Minimum percentage UV peak area required for inclusion in analysis
* -min_massconf_threshold (-mmt): Integer value
	The lowest mass confidence (percentage height of m/z peak divided by the sum of percentage heights of all observed m/z peaks) at which a hit will still be recorded. 
* -min_uv_threshold (-mut): Integer value
	Minimum height of a peak in the UV absoption profile (AU vs wavelength) required to be saved for use in the hit validation algorithm.
* -uv_match_threshold (-umt): Float value
	Parameter to determine how many UV maxima a hit needs to match those found for the bulk values before it is considered a tentative hit and flagged to the user for further analysis.  
* -uv_cluster_threshold (-uct): Float value
	Parameter to determine how many times a UV maximal wavelength should be observed for a compound (across all hits in a cluster) compared to the number of times the most commonly observed wavelength is seen, as a fraction, to be considered a “required value” for that compound. 
* -massconf_theshold (-mt): Float value
	The value for which, when it’s multiplied by the mean mass confidence for a cluster, the mass confidence for a hit in that cluster is deemed too small and is flagged to the user for further analysis.
* -cluster_size_threshold (-cst): Float value
	The value for which, when it’s multiplied by the number of hits in the largest cluster, the number of hits in a cluster is deemed too small and is discarded from further analysis. 
* -min_no_of_wells (-now): Integer value
	Minimum number of hit-containing wells required for a compound to allow the full validation process to run.
* -mass_or_area (-moa): “mass_conf” or “area”
	Parameter to determine how which hit is chosen when more than one hit is found in a well for a compound
* -calc_higherions (-chi): "True" or "False"
	Use [M+2H]2+ and [M+3H]3+ ions when looking for hits for a compound. Specify when a set of compounds 
	are likely to ionise multiple times on the mass spectrometer. 





Parameters Used for the Visualisations
---------------------------------------

* -points_per_trace (-ppt): Integer value
	The maximum number of data points of a UV chromatogram to save. This value should be reduced if the program is found to run exceptionally slowly. 