--------------------------------------------------------------------
This directory contains programs and data files needed to add a 
NEW X-section to UCI reference J-code(v8) AND then fast-JX (v5.3).
--------------------------------------------------------------------

This is a three-step process:
   (1) generate X-sections at 10 /cm (wavenumber) resolution,
         see example.  This is your responsibility!
   (2) run 'add_xsect.f' on output of (1) to generate flux-weighted
         X-sections in the 77-bins used in UCI std chemistry
   (3) run 'add_FJX.f' on output of (2) to generate fast-JX tables.

--------------------------------------------------------------------
add_xsect.f 
   reads in data sets: sflx_10cm.dat & wbins_77.dat
   reads in user-supplied Xsections at 10/cm resolution
         e.g. given as NO2-IUPAC.in
   writes (stdout) the 77-bin X-sections
         e.g. given as NO2-IUPAC.out

>>> We are no longer supplying our 10 /cm fits to all the
>>> cross-sections, since this is redundant with the JX X-ections.

>>> Further, the choices of interpolation and extrapolation  
>>> needed to generate these are a scientific choice, one that 
>>> the user needs to make.


--------------------------------------------------------------------
add_FJX.f
   reads in data set:  FJX_bins.dat
   uses include file:  cmn_bins.f
   reads in (stdin) the output from add_xsect
         e.g. given as NO2-IUPAC.out
   writes (stdout) the pair of 18-bin fast-JX tables
         e.g. given as NO2_FJX.out

--------------------------------------------------------------------
Another program 'FJX_bin.f' is meant to read the full spectral data 
   file from the UCI ref J-code (fort10.x) and generate the full 
   fast-JX spectral data set (FJX_spec.dat) <<<<unlikely needed>>>>

--------------------------------------------------------------------




