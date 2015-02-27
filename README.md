# TDTR_vH3_pub
Efficient, flexible TDTR data analysis in MATLAB

Hi all,

This is the next version of my MATLAB scripts for TDTR data analysis, after vH2. The TDTR_vH2_pub repository is still up on GitHub, so feel free to check that one out for past README's and documentation. For vH3, the introduction is handled with the TDTR_vH3 slides PDF file I've posted in this repository. Please look at those, and then at the process and analyze MATLAB scripts, to get started with my system.

What's new in vH3 from vH2? Mainly:

* Moved around the inputs to the writeLCTE database function, since some were obsolete or irritatingly placed.
* Cleaned up MAIN script by partitioning into sub-scripts.
* Can now model, fit, and do sensitivities for beam offset TDTR data.
* Added variables to some of the cellparam bundles to handle beam offset data.
* Automatic fitting is more flexible, to keep up with the capabilities of manual fitting.
* New functionality: parametric sensitivity plots, where sensitivities relative to a specified signal component (ratio, Vin, Vout, beam offset, FWHM of beam offset) can be plotted at a specified time delay or spatial offset, as a function of a thermal or system model parameter (transducer thickness, heat capacity, and so on).

This capability is controlled through the sense_setup.m or sense_setup_retro.m script, and I've included several other sense_setup scripts that I used to validate the parametric sensitivity plots against prior publications from our group.

Warning: not every function in these scripts has been re-tested since vH3. In particular, the user should take care with errorbar and parametric sensitivity calculations relative to beam offset data, to ensure the code is running properly.

The two-temperature TDTR scripts, TDTR_TTM_vH2 as a whole, have NOT been updated since vH2.

Again, look through the posted slides for guidance, try out some of the customized process and analyze scripts using the included data and results folders, and examine the template process and analyze scripts to get started. You may contact me if you have questions about the scripts or about how to expand its capabilities to fit your application.

-Greg Hohensee
