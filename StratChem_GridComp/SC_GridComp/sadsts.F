      SUBROUTINE sadsts(temp,press,pi,h2ogas,hno3gas,hno3cond,
     &                  stssad,rmed,denssts,constantn,nsts,sigsts,rsts,
     &                  sasulf,h2so4gas,rlbs,siglbs,nlbs,constantnlbs,
     &                  h2o_clim,mdens)

c... subroutine:   sadsts
c
c... written by: David B. Considine
c...             Earth System Science Interdisciplinary Center and
c...             Department of Meteorology
c...             University of Maryland
c...             College Park, MD 20742
c
c... email:      dbc@welkin.gsfc.nasa.gov
c... phone:      (301) 405-5282
c... 
c... date:       3/5/99  
c
c... purpose:
c...     
c...    This subroutine is written for the GSFC 3D CTM.  It is a modified
c...    version of a subroutine originally written for the Global Modeling
c...    Initiative of the Atmospheric Effects of Aviation
c...    Program.  It is part of the heterogeneous chemistry
c...    package to be used in the GSFC 3D CTM.
c
c... description:  
c
c...    Subroutine sadsts calculates the surface area density
c...    of supercooled ternary sulfate aerosols when called.
c...    It uses the equilibrium model of Ken Carslaw to calculate
c...    the total volume density of the STS aerosol and 
c...    convertss the volume density to surface area density
c...    based on size distribution parameters passed to the
c...    code from the calling subroutine.  The code takes as
c...    input an observed sulfate aerosol surface area density
c...    from observations.  It uses a sulfate aerosol size 
c...    distribution assumption to converts from surface area
c...    density to volume density.  It then takes this volume density
c...    field and uses binary H2SO4/H2O assumptions to calculate
c...    the condensed phase H2SO4 amount consistent with the
c...    sulfate aerosol surface area density.  This H2SO4 amount
c...    is then used along with HNO3 and H2O to calculate the
c...    surface area consistent with a ternary composition assumption.
c...    NOTE: the field stssad is the DIFFERENCE between the calculated
c...    surface area density assuming sts composition and the background
c...    sulfate surface area density, when the difference produces a
c...    positive result.  This allows the user of the code to simply
c...    sum the background and sts surface area densities to get a
c...    valid total.

c... variable declarations
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

c... declare input variables
      LOGICAL :: constantn,constantnlbs

      REAL(KIND=DBL) :: hno3gas, h2ogas, hno3cond
      
      REAL :: temp, press, mdens, h2so4gas, sasulf, nsts, sigsts
      REAL :: rsts, rlbs, siglbs, nlbs, zt, h2o_clim, pi

c... declare output variables
      REAL :: stssad,rmed

c... declare internal variables
      LOGICAL :: err

      REAL(KIND=DBL) :: h2so4ppbv, ws, wn, v, parthno3
      REAL(KIND=DBL) :: logsigstssq, rstscm, denssts
      REAL(KIND=DBL) :: h2oppmv, hno3ppbv, hno3amb

c... set up call to Carslaw's code (subroutine ternary)
c... get total (gas + condensed phase HNO3)
      hno3amb = hno3gas+hno3cond

c... converts h2o mixing ratio to ppmv:
      h2oppmv = h2ogas*1.e6

c... converts hno3 concentration to ppbv.  Since we are doing an equilibrium
c... calculation here, pass the sum of gas and condensed phase HNO3
c... to Carslaw's code.  Carslaw's code will decide partitioning between
c... gas and condensed phases.
      hno3ppbv = hno3amb*1.e9

c... get h2so4 concentration in ppbv from liquid binary sulfate
c... surface area density

c... Note added 7/7/02: This code uses
c... a background h2o climatology, which is passed to calcsulf
c... here.  The way calcsulf works, passing a strongly dehydrated
c... h2o distribution to calcsulf would be a mistake.  This is
c... because the routine uses a climatological background sulfate
c... surface area density distribution to infer h2so4 mixing ratios.
c... Given a particular sad, the less water at that location, the
c... more inferred h2so4.  Thus, the routine could return unrealistically
c... high h2so4 mixing ratios, leading to substantial overestimates
c... of the amount of ternary calculated below.
c      h2so4gas = 0.0

      call calcsulf (temp,press,h2o_clim,sasulf,h2so4gas,rlbs
     &  ,siglbs,nlbs,constantnlbs,mdens,pi)

      h2so4ppbv = h2so4gas*1.e9

c... now call subroutine ternary:
      ws = 0.0
      wn = 0.0
      v = 0.0
      parthno3 = 0.0
      err = .false.
      denssts = 0.0
     
      call ternary(temp,press,h2oppmv,hno3ppbv,h2so4ppbv,ws,wn,v,
     &     parthno3,err,denssts,pi)

c... parthno3 is the fraction of the hno3 molecules that are in gas phase.
c... Use this to calculate gas and condensed phase HNO3 amounts
c... I have noticed that parthno3 is sometimes slightly greater than 1,
c... which could create negative hno3cond values.  Trap for this 
c... condition 
      if(parthno3.gt.1.0) parthno3 = 1.0
      hno3gas = hno3amb*parthno3
      hno3cond = hno3amb-hno3gas

c... converts sts volume to surface area density
c... use constant STS particle number density assumption if flag set
      logsigstssq = log(sigsts)*log(sigsts)
      rstscm = rsts*1.e-4

c... converts volume units from microns cubed/cm**3 to cm**3/cm**3
      v = v*1.e-12

      if(constantn) then
         stssad = (4.0*pi*nsts)**(1.0/3.0)*
     &      (3.0*v)**(2.0/3.0)*exp(-1.0*logsigstssq)
         rmed = (3.0*v/(4.0*pi*nsts))**(1.0/3.0)
     &      *exp(-3.0/2.0*logsigstssq)
       else
         stssad = 3.0/rstscm*v
     &      *exp(-5.0*logsigstssq/2.0)
       endif

c... return the difference between the calculated sts sad and sasulf,
c... providing stssad > sasulf
      stssad = stssad-sasulf

      if(stssad.lt.0.0) then
        stssad = 0.0
        rmed = 0.0
        hno3cond = 0.0
        hno3gas = hno3amb
      endif

      RETURN
      END SUBROUTINE sadsts
