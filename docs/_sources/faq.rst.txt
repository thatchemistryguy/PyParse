Frequently Asked Questions
===============================


Which analysis type should I choose?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	This depends on your plate and the information you want to extract – speak to a lead user or super 
	user for advice, or refer to the “Help” tab in Pipeline Pilot for detail on each analysis type. 


My plate design has deliberate gaps/empty wells in it - what do I do?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	Continue generating the platemap as per the guide on :ref:`preparing_a_platemap_label`, omitting
	any wells which shouldn't be used in the analysis. 

I'd like to view the results in Excel/Spotfire/another program - how do I do this?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	The compoundtable.csv and outputTable.csv output files contain all the information used
	to generate the visualisations detailed in the section on :ref:`pyparse_outputs_label`. 
	Import these files into your preferred data analysis program for custom visualisations.
	
The program didn't find my product!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Check that the compound strongly ionises to give an [M+H]+, [M+2H]2+ or [M+3H]3+ ion, or Boc group degradation fragment ions. 
PyParse uses a generous limit to ensure that background MS noise isn't captured as a compound. If the compound
doesn't ionise, add the retention time to the platemap to force PyParse to accept peaks at that time which don't ionise correctly
:ref:`adding_a_rt_label`


The output says “Peak overlap detected!” for my compound – what does this mean?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
PyParse found that in the Best Well for that product, the compound peak overlapped with another 
peak in the LCMS. The calculated purity given in the “Purity of Best Well” column may therefore be 
incorrect. Check the LCMS yourself, and re-run using a different modifier (typically formic acid, ammonium bicarbonate or 
trifluoroacetic acid).

 