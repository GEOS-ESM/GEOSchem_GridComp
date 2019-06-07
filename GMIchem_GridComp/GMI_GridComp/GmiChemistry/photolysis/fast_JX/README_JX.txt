-----------------------------------------------
Fast-JX Photolysis Code, V 2.  (XZ & MP, 7/04)
-----------------------------------------------

Fast-JX scheme is new:
      it calculates photolysis rates in both troposphere AND stratosphere;
      the new spectral data has been updated to JPL-02 and IUPAC(2004);
      see notes on jx_spec_04a.dat.	

The main features of fast-J (troposphere, wvl > 291 nm) are combined with fast-J2
(includes stratosphere, all J's to 60 km).  The new fast-Jx now calculates short
wavelengths the same as long wavelengths:  
	1) pseudo-spherical attenuation of the solar beam (curved atmos)
	2) multiple scattering in a plane parallel atmosphere
	3) full scattering phase function (+ 8-stream anisotropic)
	4) optimized wavelength bins:  7 bins 291-850 nm (troposphere)
                                     +11 bins 180-291 nm (stratosphere)

Fast-JX has been compared with the full 77-bin (145 bins counting S-Runge bins)
standard UCI photolysis code.  Both for all altitudes (18 JX-bins to do 0-60 km)
and for the Fast-J version for troposphere (use last 7 bins to do 0-18 km).
This was done for SZA=0 and the US Std atmosphere. Agreement is within 3% over
all regions of relevance (e.g., J for CFCs only in stratosphere). 
(The lowest level does not match, since Fast-JX reports the middle of that
level (~0.5 km) and the std code reports the surface (0 km), but all else
looks fine!).

Features:
  >>> the new fast-JX release can be used for troposphere only CTMs with the 
      same code and spectral dataset:  select wavelength bins 12:18 in the 
      jx_spec.dat file, no other changes needed.

    > O3 climatology integrated rather than interpolated, to ensure
      preservation of the supplied column. (An updated climatology needed)

    > Simplified level scheme to allow easy replacement of climatological
      O3 and T with model-supplied values.

    > Supports scattering from more than one type of aerosol at once;
      data file now includes phase functions for 14 types.

    > Finer tau-steps at cloud tops to improve accuracy at higher solar
      zenith angles.

    > Alternate flux attenuation scheme using second-order finite
      difference method - better for calculating heating rates.

    > Additional diagnostic routines and comments to aid in understanding


The current release comes in the form used in the UCI CTM, with a model framework 
supplied so that the code can be run in "stand-alone" mode, with the parameters 
set in 'cmn_h.t' and user controls in 'input' - both files contain additional 
comments. The following files  are required:

 Fast-JX bits:

     newcol.f        Fast-JX CTM code

     jv_cmn.h        Fast-JX include file - main variables and commons
     jv_mie.h        Fast-JX include file - Mie variables

     jx_spec.dat     Fast-JX input - spectral data file
                  **************************************************************
                  ***These are the new 2004 updates (designated 04b) to all the 
                  *** cross-sections, including some IUPAC ones.  see the notes
                  *** at the end of the table to calculate branching rates and
                  *** dissociation J's for NO3, acetone, MeVK.
                  **************************************************************

     jv_atms.dat     Fast-JX input - T and O3 climatology (dated, M&M-I data)
                     This file is not essential except to fill in above CTM levels,
                     you can write your own. 

     ratj.d          Specimen file for supplying photolysis reactions,
                       factors and labels  - keyed to Wild etal's UCI CTM
                  **************************************************************
                  ***this is extremely obscure and set for UCI's chem code
                  ***as enclosed it matches the jx_spec.dat cross-sections
                  *** and will print out all of them 
                  **************************************************************

 Model bits:

     standalone.f    Fortran routine to simulate model environment;
                       sets up the spatial and temporal variables

     cmn_h.f         Simulated header files 'included' in model code,
     cmn_t.f           and including documentation on the parameters
     cmn_w.f           and variables required

     input           Control file for run

 Makefile/Run:

     Just compile standalone.f + newcol.f together and run.

     The standalone code will write the calculated photolysis rates at each
     supplied model level to a file 'output', and write diagnostics and/or
     climatologies to stdout (unit=6).

------------------------------------------------------------
 2004-b version cross-section notes for pratmo and fast-JX
------------------------------------------------------------

The x-section data have been revised using JPL-02, IUPAC (up to 2004)
and Gierczak's acetone tables.  

(1) The NO2 temperature dependence in JPL-02 and before appears to be in
error, producing large x-sections for colder temperatures - opposite to
that in the IUPAC tabulations.  The systematic error is about 2%, smaller
than the JPL vs IUPAC recommendations at 298K (J-IUPAC > J-JPL by 5-7%).
>>>>have decided to use the JPL 298K values for ALL temperatures.

(2) The wavelength-density dependence of the quantum yield for dissociation
of acetone (CH3C(O)CH3 aka C3H6O) is awkward.  Have devised and tested
the new method of computing J at 1000 mbar and at 300 mbar and then
interpolating in log[M #/cm3]. The choice of 300 mbar was chosen to
minimize the error over the range 100-1000 mbar.  See notes at end
of jx_spec.dat

(3) MethylVinylKetone has a pressure dependent quantum yield, the cross
sections used are for the low density limit (M = 0) and the rational
polynomial formula is given in the jx_spec.dat notes.

(4) NO3 cross-sections are combined into a total with fixed branching factors.

(5) ClONO2 dissociation yields (Cl+NO3 vs. ClO+NO2) are very wavelength 
dependent and change with large zenith angle (e.g., polar twilight), hence
the new method give x-sections for both branches separately.

(6) HO2NO2 (aka HNO4) photolysis in the near IR at a 'daytime' clear-sky
rate of 1E-5 /sec has been added as a cross-section at the longest wavelength 
(412-850 nm).  This may overestimate the effects of Rayleigh and other 
scattering:  e.g., a clear-sky, 10% albedo enhancement is 1.3E-5.

(7) Some of the original Harvard x-sections for methacrolein, glyoxal,
etc have been dropped, because I could not find back in the IUPAC or JPL
listings.  Many of these are too slow to be important and should be dropped.
The methacrolein (MACR) is included, but has only a 'guesstimate' for the
quantum yield of about 0.008(??).

(8) NEW cross-sections can be added easily as before.  The key is to use
the fast-J2 bins (these are unchanged) and to generate x-sections by 
averaging with solar-flux weighting.

-----------------------------------------------------------

                                                Oliver Wild   (29/07/99) 
                                                Huisheng Bian (29/08/02)
                                                Michael Prather (Jul'04)