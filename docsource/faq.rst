Frequently Asked Questions
===============================

Which UPLC/LCMS Vendors are Supported?
-----------------------------------------

	Currently, the .rpt file format from Waters\ |trademark| and the .daml file format from Shimadzu.

Which plot type heatmap should I look at?
--------------------------------------------

	This depends on your plate and the information you want to extract. Typically, the heatmap 
	plotted by product percentage area (Parea) is useful where there is no internal standard present. 
	In cases where an internal standard was added to the plate, the heatmap plotted by the corrected ratio of 
	product to internal standard (corrP/STD) is most useful. 

What does "corrP/STD" actually mean?
-----------------------------------------

	This is the ratio of the absolute peak of the product to internal standard in that well, 
	divided by the maximum observed ratio of the product to internal standard observed in any well 
	for that particular product. This normalises arbitrary ratio values to onto a "zero-to-one" scale. 

My plate design has deliberate gaps/empty wells in it - what do I do?
----------------------------------------------------------------------

	Continue generating the platemap as per the guide on :ref:`preparing_a_platemap_label`, omitting
	any wells which shouldn't be used in the analysis. 

I'd like to view the results in Excel/Spotfire/another program - how do I do this?
------------------------------------------------------------------------------------

	The compoundtable.csv and outputTable.csv output files contain all the information used
	to generate the visualisations detailed in the section on :ref:`pyparse_outputs_label`. 
	Import these files into your preferred data analysis program for custom visualisations.
	
The program didn't find my product!!
--------------------------------------

Check that the compound strongly ionises to give an [M+H]+, [M+2H]2+ or [M+3H]3+ ion, or Boc group degradation fragment ions. 
PyParse uses a generous limit to ensure that background MS noise isn't captured as a compound. If the compound
doesn't ionise, add the retention time to the platemap to force PyParse to accept peaks at that time which don't ionise correctly
:ref:`adding_a_rt_label`


The output says “Peak overlap detected!” for my compound – what does this mean?
----------------------------------------------------------------------------------

PyParse found that in the Best Well for that product, the compound peak overlapped with another 
peak in the LCMS. The calculated purity given in the “Purity of Best Well” column may therefore be 
incorrect. Check the LCMS yourself, and re-run using a different modifier (typically formic acid, ammonium bicarbonate or 
trifluoroacetic acid).


.. |trademark|	unicode:: U+2122 .. TRADEMARK SYMBOL

 