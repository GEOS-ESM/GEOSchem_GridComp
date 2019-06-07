-----------------------------------------------------------------
The UCI fast-JX code has been updated to version 5.3c (MJP 4/2006)
-----------------------------------------------------------------

1) The solar spectrum has been updated ! from the old WMO 1985
  heritage.  The new one is SUSIM average of 29Mar92 & 11Nov94.
  This new spectrum is very similar to SOLAR2000 model average
  over a a couple decades (Roland Uhl, p.c.) but is very high
  resolution. 
     a) the 77-bin ref spectrum and the fast-JX solar bins are updated
     b) the fast-JX X-sections have been redone from the UCI ref code
     c) the old high-res (10/cm) solar spectrum in /add_X is NOT changed
        hence the original 77-bin reference X-sections are unchanged
     d) THUS, the /add_X solar spectum STILL needs to be updated to SUSIM.

2) The std atmospheres have been updated to Gordon Labow's (w/ J Logan) 
   new O3 climatology.

-----------------------------------------------------------------
The UCI fast-JX code has been updated to version 5.3b (MJP 7/2005)
-----------------------------------------------------------------

The new code has many advantages, compares closely (under tests
so far) with the old code.  The codes for adding a new X-section
have been cleaned and corrected (bug found in averaging X-sects
had impact at <1% (max)) - please see the sub-dir /add_X.

Fixes compared with previous versions:
      > sub-b code has vastly cleaned up the interfaces and is
          now ready for parallel implementation.  Also the 
          standalone version is simpler!
      > vastly simplified calling sequences, reduced use of commons
      > simplified input datasets
      > updated cross-sections AND code for generating new X-sections.
      >>> new, improved method for adding layers for thick clouds:
            uses log-spacing factor of 1.18 (error in tests <0.3%)
            reduces the number of layers (linear with comp costs)
            log-spacing eliminates need for upper bound on cloud OD
            eliminated some unnecessary bits of previous fast-J

The common blocks (all local) needed for Mie scattering have been
replaced with local internal arrays passed through calls.

The new code has been implemented in UCI's CTM with openMP.

-----------------------------------------------------------------
The new standalone version of fast-JX (fast-JX53.f) requires:
      parm_CTM.f     parameters for CTM dimensions
      parm_MIE.f     parameters for Mie scattering arrays
      cmn_metdat.f   common blocks w/ general CTM data for fast-JX
      cmn_JVdat.f    common blocks w. specific data for J-values

The needed data sets standalone runs are included:
      FJX_spec.dat   spectral data (X-sections+) ver 2005e
      FJX-scat.dat   scattering data for clouds+aerosols 
            this data set has been separated from the X-sections

  these are for the standalone run and will not be needed in CTM run
      FJX_test.dat   set up for p, T, grid etc
      atmos_std.dat  UCI's old std atmospheres
      ECT42_grid.dat used for sample run
      ratj.dat       UCI's old chemistry set of J's, 
                        shows mapping from fast-JX
      
      std.out        output from this test of fast-JX53 standalone

-----------------------------------------------------------------

