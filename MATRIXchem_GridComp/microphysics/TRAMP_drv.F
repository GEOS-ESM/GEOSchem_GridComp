#ifndef GEOS5_PORT
#include "rundeck_opts.h"
#endif

      MODULE AMP_AEROSOL
!@sum Driver for Aerosol Microphysics
!@auth Susanne Bauer
#ifndef GEOS5_PORT
      USE TRACER_COM
      USE AERO_CONFIG, ONLY: NMODES
      USE AERO_PARAM,  ONLY: NEMIS_SPCS
      USE MODEL_COM,   ONLY: LM, jhour, jdate
#else
      USE AERO_CONFIG, ONLY: NMODES
      USE AERO_PARAM,  ONLY: NEMIS_SPCS
#endif
      IMPLICIT NONE
      SAVE

C**************  Latitude-Dependant (allocatable) *******************
#ifndef GEOS5_PORT
      ! Mie lookup tables
      REAL*8, DIMENSION(15,17,23,6)      :: AMP_EXT, AMP_ASY, AMP_SCA   !(15,17,23,6) (RE,IM,size,lambda)
      REAL*8, DIMENSION(15,17,23)        :: AMP_Q55
      REAL*8, DIMENSION(23,26,26,26,6)   :: AMP_EXT_CS,AMP_ASY_CS,AMP_SCA_CS !(23,26,26,26,6) (radius,OC,SO4,H2O,lambda)
      REAL*8, DIMENSION(23,26,26,26)     :: AMP_Q55_CS
      ! 1 Dim arrays for Radiation
      REAL*8, DIMENSION(LM,nmodes)       :: Reff_LEV, NUMB_LEV
      REAL*8, DIMENSION(LM,nmodes)       :: MIX_OC, MIX_SU, MIX_AQ
      COMPLEX*8, DIMENSION(LM,nmodes,6)  :: RindexAMP
      REAL*8, DIMENSION(LM,nmodes,7)     :: dry_Vf_LEV
      ! FALSE : one Radiation call
      ! TRUE  : nmodes Radiation calls
      INTEGER                            :: AMP_RAD_KEY = 1 ! 1=Volume Mixing || 2=Core - Shell || 3=Maxwell Garnett
      INTEGER                            :: AMP_DIAG_FC = 1 ! 2=nmode radiation calls  ||  1=one radiation call

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)       :: AQsulfRATE !(i,j,l)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: DIAM       ![m](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: AMP_dens   !density(i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: AMP_TR_MM  !molec. mass(i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: NACTV      != 1.0D-30  ![#/m^3](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: CCNSS      != 1.0D-30  ![#/m^3](i,j,l,nmodes,3)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: VDDEP_AERO != 1.0D-30  ![m/s](i,j,nmodes,2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: NUMB_SS  ! Sea salt number concentration [#/gb]
#else
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: DIAM       ![m](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: NACTV      != 1.0D-30  ![#/m^3](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:)   :: CCNSS      != 1.0D-30  ![#/m^3](i,j,l,nmodes,3)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: VDDEP_AERO != 1.0D-30  ![m/s](i,j,nmodes,2)
#endif

#ifndef NO_HDIURN
      REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: DIURN_LWP  ! lwp hourly output
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)     :: DIURN_LWC  ! lwc hourly output
#endif
!-------------------------------------------------------------------------------------------------------------------------
!     The array VDDEP_AERO(X,Y,Z,I,1) contains current values for the dry deposition velocities 
!     for aerosol number concentrations for mode I. 
!     The array VDDEP_AERO(X,Y,Z,I,2) contains current values for the dry deposition velocities 
!     for aerosol mass   concentrations for mode I. 
!     Values in VDDEP_AERO are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!     The array NACTV(X,Y,Z,I) contains current values of the number of aerosol particles 
!     activated in clouds for each mode I for use outside of the MATRIX microphysical module.
!     Values in NACTV are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------      
!     1 - BC  2-BCmix 3-OC 4-OCmix 5-SS1  6-SS2 7-D1 8-D2
!-------------------------------------------------------------------------------------------------------------------------      
!-------------------------------------------------------------------------------------------------------------------------
!     The array DIAM(x,y,z) contains current values of some measure of average ambient mode diameter for each mode  
!       for use outside of the MATRIX microphysical module where it is calculated. 
!       Values in DIAM are saved at the top of the subr. MATRIX before microphysical evolution 
!       for the current time step is done. 
!
!     The current measure of particle diameter is the diameter of average mass:
! 
!        DIAM(x,y,z) = [ (6/pi) * (Mi/Ni) * (1/D) ]^(1/3)
!
!     with Mi the total mass concentration (including water) in mode i, Ni the number concentration in mode i, and
!     D a constant ambient particle density, currently set to D = DENSP = 1.4 g/cm^3. 
!-------------------------------------------------------------------------------------------------------------------------      
      END MODULE AMP_AEROSOL

#ifndef GEOS5_PORT
      SUBROUTINE MATRIX_DRV
      USE TRACER_COM
      USE TRDIAG_COM, only : taijs=>taijs_loc,taijls=>taijls_loc
     *     ,ijts_AMPp,ijlt_AMPm,ijlt_AMPext,ijts_AMPpdf
     *     ,itcon_AMP,itcon_AMPm
      USE AMP_AEROSOL

      USE MODEL_COM, only : im,jm,lm     ! dimensions
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
     $                     ,dtsrc
      USE GEOM, only: axyp,imaxj,BYAXYP
      USE CONSTANT,   only:  lhe,mair,gasc   
      USE FLUXES, only: tr3Dsource,trsource,trsrfflx,trflux1
      USE DYNAMICS,   only: pmid,pk,byam,gz, am   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           BYAM  1/Air mass (m^2/kg)
      USE AERO_CONFIG
      USE AERO_INIT
      USE AERO_PARAM, only: IXXX, IYYY, ILAY, NEMIS_SPCS
      USE AERO_SETUP 
      USE PBLCOM,     only: EGCM !(LM,IM,JM) 3-D turbulent kinetic energy [m^2/s^2]
      USE DOMAIN_DECOMP_ATM,only: GRID, GET, am_i_root
#ifndef NO_HDIURN
c for the hourly diagnostic
#ifdef CLD_AER_CDNC 
      USE CLOUDS_COM, only: CDN3D  ! CDNC
#endif
      USE DIAG_COM, only: adiurn=>adiurn_loc,ndiuvar,iwrite,jwrite,itwrite,ndiupt,idd_diam
     *                    ,ijdd, idd_ccn, idd_cdnc, idd_lwp, idd_numb, idd_mass, idd_so2,
     *                     idd_lwc, idd_ncL, idd_pres
     *     ,hdiurn=>hdiurn_loc
#endif
#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups,
     &     inc_subdd,find_groups
#endif
      IMPLICIT NONE

      REAL(8):: TK,RH,PRES,TSTEP,AQSO4RATE
      REAL(8):: AERO(NAEROBOX)     ! aerosol conc. [ug/m^3] or [#/m^3]
      REAL(8):: GAS(NGASES)        ! gas-phase conc. [ug/m^3]
      REAL(8):: EMIS_MASS(NEMIS_SPCS) ! mass emission rates [ug/m^3]
      REAL(8):: SPCMASS(NMASS_SPCS+2)
      REAL(8):: DT_AERO(NDIAG_AERO,NAEROBOX) !NDIAG_AERO=15
      REAL(8):: yS, yM, ZHEIGHT1,WUP,AVOL 
      REAL(8) :: PDF1(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF2(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      INTEGER:: j,l,i,n,J_0, J_1, I_0, I_1, m,nAMP
C**** functions
      REAL(8):: QSAT
#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      type(subdd_type), pointer :: subdd

      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &    PN_subdd, APN_subdd
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &    PN_3d_subdd, APN_3d_subdd,p_3d_subdd, rh_3d_subdd,
     &    c1_3d_subdd, c2_3d_subdd, c5_3d_subdd

#endif
#ifndef NO_HDIURN
c for the hourly diagnostic
      REAL*8 :: TMP(NDIUVAR)
      INTEGER, PARAMETER :: NLOC_DIU_VAR  = nmodes + nmodes + 38 +  2
      INTEGER, PARAMETER :: NLOC_DIU_VARL = 5
      INTEGER :: idxd(NLOC_DIU_VAR),idxl(NLOC_DIU_VARL) 
      INTEGER :: ih, ihm, kr,ii
#endif
      REAL*8 :: HD_NUMB(nmodes) 

      CALL GET(grid, J_STRT =J_0, J_STOP =J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      NACTV(I_0:I_1,J_0:J_1,:,:)      = 0.d0 
      CCNSS(I_0:I_1,J_0:J_1,:,:,:)      = 0.d0 
      VDDEP_AERO(I_0:I_1,J_0:J_1,:,:) = 0.d0 
      DIAM(I_0:I_1,J_0:J_1,:,:)       = 0.d0
      AMP_dens(I_0:I_1,J_0:J_1,:,:)   = 0.d0
      AMP_TR_MM(I_0:I_1,J_0:J_1,:,:)  = 0.d0
      NUMB_SS(I_0:I_1,J_0:J_1,:,:)    = 0.d0
#ifdef CACHED_SUBDD
      PN_subdd(I_0:I_1,J_0:J_1)    = 0.d0
      APN_subdd(I_0:I_1,J_0:J_1)    = 0.d0
      PN_3d_subdd(I_0:I_1,J_0:J_1,:)    = 0.d0
      APN_3d_subdd(I_0:I_1,J_0:J_1,:)    = 0.d0
      c1_3d_subdd(I_0:I_1,J_0:J_1,:)    = 0.d0
      c2_3d_subdd(I_0:I_1,J_0:J_1,:)    = 0.d0
      c5_3d_subdd(I_0:I_1,J_0:J_1,:)    = 0.d0
      P_3d_subdd(I_0:I_1,J_0:J_1,:)    = 0.d0
      RH_3d_subdd(I_0:I_1,J_0:J_1,:)    = 0.d0
#endif
      DO L=1,LM                            
      DO J=J_0,J_1                          
      DO I=I_0,I_1                 

      IXXX = I
      IYYY = J
      ILAY = L
      DT_AERO(:,:) = 0.d0
      EMIS_MASS(:) = 0.d0
      AERO(:)      = 0.d0
      HD_NUMB(:)   = 0.d0
! meteo
      TK = pk(l,i,j)*t(i,j,l)           !should be in [K]
      RH = MIN(1.,q(i,j,l)/QSAT(TK,lhe,pmid(l,i,j))) ! rH [0-1]
      PRES= pmid(l,i,j)*100.                  ! pmid in [hPa]
      TSTEP=dtsrc
      ZHEIGHT1 = GZ(i,j,l) /1000./9.81
      WUP = SQRT(.6666667*EGCM(l,i,j))  ! updraft velocity

c avol [m3/gb] mass of air pro m3      
      AVOL = am(l,i,j)*axyp(i,j)/mair*1000.d0*gasc*tk/pres 
! in-cloud SO4 production rate [ug/m^3/s] ::: AQsulfRATE [kg] 
      AQSO4RATE = AQsulfRATE (i,j,l)* 1.d9  / AVOL /dtsrc
c conversion trm [kg/gb] -> [ug /m^3]
      GAS(1) = trm(i,j,l,n_H2SO4)* 1.d9 / AVOL! [ug H2SO4/m^3]
c conversion trm [kg/kg] -> [ug /m^3]
      GAS(2) = trm(i,j,l,n_HNO3)*1.d9 / AVOL!   [ug HNO3/m^3]
c conversion trm [kg/gb] -> [ug /m^3]
      GAS(3) = trm(i,j,l,n_NH3)* 1.d9 / AVOL!   [ug NH3 /m^3]
!  [kg/s] -> [ug/m3/s]

       DO n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
c conversion trm [kg/gb] -> AERO [ug/m3]
         if(AMP_NUMB_MAP(nAMP).eq. 0) then
       AERO(AMP_AERO_MAP(nAMP)) =trm(i,j,l,n)*1.d9 / AVOL ! ug/m3
          else

       AERO(AMP_AERO_MAP(nAMP)) =trm(i,j,l,n)/ AVOL       !  #/m3
          endif
       ENDDO

       if (L.eq.1) then     
!      Emis Mass [ug/m3/s] <-- trflux1[kg/s]
#ifdef TRACERS_AMP_M4
      EMIS_MASS(2) =  MAX(trflux1(i,j,n_M_ACC_SU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(3) =  MAX(trflux1(i,j,n_M_BC1_BC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(4) =  MAX(trflux1(i,j,n_M_OCC_OC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(5) =  MAX(trflux1(i,j,n_M_DD1_DU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(6) =  MAX(trflux1(i,j,n_M_SSS_SS)*1.d9 / AVOL,0.d0)
      EMIS_MASS(10)=  MAX(trflux1(i,j,n_M_DD2_DU)*1.d9 / AVOL,0.d0)
#else
      EMIS_MASS(1) =  MAX(trflux1(i,j,n_M_AKK_SU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(2) =  MAX(trflux1(i,j,n_M_ACC_SU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(3) =  MAX(trflux1(i,j,n_M_BC1_BC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(4) =  MAX(trflux1(i,j,n_M_OCC_OC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(5) =  MAX(trflux1(i,j,n_M_DD1_DU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(6) =  MAX(trflux1(i,j,n_M_SSA_SS)*1.d9 / AVOL,0.d0)
      EMIS_MASS(7) =  MAX(trflux1(i,j,n_M_SSC_SS)*1.d9 / AVOL,0.d0)
      EMIS_MASS(10)=  MAX(trflux1(i,j,n_M_DD2_DU)*1.d9 / AVOL,0.d0)
#endif
        endif
!      Emis Mass [ug/m3/s] <-- trflux1[kg/s]
#ifdef TRACERS_AMP_M4
      EMIS_MASS(2) =  EMIS_MASS(2) + ((tr3Dsource(i,j,l,nVolcanic,n_M_ACC_SU)+
     *                                 tr3Dsource(i,j,l,nBiomass,n_M_ACC_SU))*1.d9 / AVOL)
      EMIS_MASS(3) =  EMIS_MASS(3) + (tr3Dsource(i,j,l,nBiomass,n_M_BC1_BC)*1.d9 / AVOL)
      EMIS_MASS(9) =  EMIS_MASS(9) + (tr3Dsource(i,j,l,nBiomass,n_M_OCC_OC)*1.d9 / AVOL)
#else
      EMIS_MASS(1) =  EMIS_MASS(1) + ((tr3Dsource(i,j,l,nVolcanic,n_M_AKK_SU)+
     *                                 tr3Dsource(i,j,l,nBiomass,n_M_AKK_SU))*1.d9 / AVOL)
      EMIS_MASS(2) =  EMIS_MASS(2) + ((tr3Dsource(i,j,l,nVolcanic,n_M_ACC_SU)+
     *                                 tr3Dsource(i,j,l,nBiomass,n_M_ACC_SU))*1.d9 / AVOL)
      EMIS_MASS(3) =  EMIS_MASS(3) + (tr3Dsource(i,j,l,nBiomass,n_M_BC1_BC)*1.d9 / AVOL)
c     Biomass BC OC is Mixed
c      EMIS_MASS(8) =  EMIS_MASS(8) + (tr3Dsource(i,j,l,nBiomass,n_M_BOC_BC)*1.d9 / AVOL)
c      EMIS_MASS(9) =  EMIS_MASS(9) + (tr3Dsource(i,j,l,nBiomass,n_M_BOC_OC)*1.d9 / AVOL)
c     Biomass BC OC is NOT mixed
      EMIS_MASS(3) =  EMIS_MASS(3) + (tr3Dsource(i,j,l,nBiomass,n_M_BC1_BC)*1.d9 / AVOL)
      EMIS_MASS(4) =  EMIS_MASS(4) + (tr3Dsource(i,j,l,nBiomass,n_M_OCC_OC)*1.d9 / AVOL)
#endif
       CALL SPCMASSES(AERO,GAS,SPCMASS)

       CALL MATRIX(AERO,GAS,EMIS_MASS,TSTEP,TK,RH,PRES,AQSO4RATE,WUP,DT_AERO) 
c       CALL SIZE_PDFS(AERO,PDF1,PDF2)
 
       DO n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
          if(AMP_NUMB_MAP(nAMP).eq. 0) then
! Mass loop
      tr3Dsource(i,j,l,nChemistry,n) =((AERO(AMP_AERO_MAP(nAMP)) *AVOL *1.d-9) -trm(i,j,l,n)) /dtsrc 
          else
! Number loop
      tr3Dsource(i,j,l,nChemistry,n) =((AERO(AMP_AERO_MAP(nAMP)) *AVOL)  -trm(i,j,l,n)) /dtsrc
#ifndef NO_HDIURN
      if (DIAM(i,j,l,AMP_MODES_MAP(nAMP))*1.e6.ge.0.1 ) HD_NUMB(AMP_MODES_MAP(nAMP)) = AERO(AMP_AERO_MAP(nAMP))
#endif
#ifdef CACHED_SUBDD
C****
C**** Collect some high-frequency outputs
C****
       if (DIAM(i,j,l,AMP_MODES_MAP(nAMP))*1.e6.ge.0.1 ) HD_NUMB(AMP_MODES_MAP(nAMP)) = AERO(AMP_AERO_MAP(nAMP))
        PN_subdd(i,j) = PN_subdd(i,j) + (HD_NUMB(AMP_MODES_MAP(nAMP))) * ZHEIGHT1
        PN_3d_subdd(i,j,l) = PN_3d_subdd(i,j,l) + (HD_NUMB(AMP_MODES_MAP(nAMP))) 
#endif /* CACHED_SUBDD */
          endif   
       ENDDO
#ifdef CACHED_SUBDD
       APN_subdd(i,j) = APN_subdd(i,j) + sum(NACTV(i,j,l,:)) * ZHEIGHT1
       APN_3d_subdd(i,j,l) = APN_3d_subdd(i,j,l) + sum(NACTV(i,j,l,:))
       c1_3d_subdd(i,j,l) = c1_3d_subdd(i,j,l) + sum(CCNSS(i,j,l,:,1))
       c2_3d_subdd(i,j,l) = c2_3d_subdd(i,j,l) + sum(CCNSS(i,j,l,:,2))
       c5_3d_subdd(i,j,l) = c5_3d_subdd(i,j,l) + sum(CCNSS(i,j,l,:,3))
       p_3d_subdd(i,j,l) = PRES
       rh_3d_subdd(i,j,l) = RH
#endif
       NUMB_SS(i,j,l,1) = AERO(22) *AVOL
       NUMB_SS(i,j,l,2) = AERO(25) *AVOL ! but has only tiny number in it
    

      tr3Dsource(i,j,l,nChemistry,n_H2SO4) =((GAS(1)*AVOL *1.d-9)
     *        -trm(i,j,l,n_H2SO4)) /dtsrc 
      tr3Dsource(i,j,l,nChemistry,n_NH3)   =((GAS(3)*AVOL *1.d-9)
     *        -trm(i,j,l,n_NH3)) /dtsrc
      tr3Dsource(i,j,l,3,n_HNO3)           =((GAS(2)*AVOL *1.d-9)
     *        -trm(i,j,l,n_HNO3))/dtsrc

c       DT_AERO(:,:) = DT_AERO(:,:) * dtsrc !DT_AERO [# or ug/m3/s] , taijs [kg m2/kg(air)], byam [kg/m2]

c Update physical properties per mode
       do n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
c Diagnostic of Processes - Sources and Sincs - timestep included
          if(AMP_NUMB_MAP(nAMP).eq. 0) then  !taijs [kg/s] -> in acc [kg/m2*s]
            do m=1,7
              taijs(i,j,ijts_AMPp(m,n)) =taijs(i,j,ijts_AMPp(m,n)) +(DT_AERO(m+8,AMP_AERO_MAP(nAMP))* AVOL * 1.d-9)
              if (itcon_amp(m,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(m+8,AMP_AERO_MAP(nAMP))* AVOL * 1.d-9),itcon_amp(m,n),n)
            end do
             
          else
            taijs(i,j,ijts_AMPp(1,n)) =taijs(i,j,ijts_AMPp(1,n))+(DT_AERO(2,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(1,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(2,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(1,n),n)
            taijs(i,j,ijts_AMPp(2,n)) =taijs(i,j,ijts_AMPp(2,n))+(DT_AERO(3,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(2,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(3,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(2,n),n)
            taijs(i,j,ijts_AMPp(3,n)) =taijs(i,j,ijts_AMPp(3,n))+(DT_AERO(1,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(3,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(1,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(3,n),n)
            do m=4,7
              taijs(i,j,ijts_AMPp(m,n)) =taijs(i,j,ijts_AMPp(m,n))+(DT_AERO(m,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(m,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(m,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(m,n),n)
            end do

          endif
       select case (trname(n)) !taijs [kg * m2/kg air] -> in acc [kg/kg air]
      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ','N_DS2_1 ','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ','N_MXX_1 ','N_OCS_1 ')
c - 3d acc output
        taijls(i,j,l,ijlt_AMPm(1,n))=taijls(i,j,l,ijlt_AMPm(1,n)) + DIAM(i,j,l,AMP_MODES_MAP(nAMP))
        taijls(i,j,l,ijlt_AMPm(2,n))=taijls(i,j,l,ijlt_AMPm(2,n)) + (NACTV(i,j,l,AMP_MODES_MAP(nAMP))*AVOL*byam(l,i,j)/axyp(i,j))

c - 2d PRT Diagnostic
        if (itcon_AMPm(1,n) .gt.0) call inc_diagtcb(i,j,(DIAM(i,j,l,AMP_MODES_MAP(nAMP))*1d6),itcon_AMPm(1,n),n) 
        if (itcon_AMPm(2,n) .gt.0) call inc_diagtcb(i,j,NACTV(i,j,l,AMP_MODES_MAP(nAMP))*AVOL ,itcon_AMPm(2,n),n) 
       end select

      enddo !n
c - special diag: Size distribution pdfs
c       if (l.eq.1) taijs(i,j,ijts_AMPpdf(l,:))=taijs(i,j,ijts_AMPpdf(l,:)) + (PDF1(:)*AVOL*byam(l,i,j))
c - N_SSA, N_SSC, M_SSA_SU
        taijls(i,j,l,ijlt_AMPext(1))=taijls(i,j,l,ijlt_AMPext(1)) + (NACTV(i,j,l,SEAS_MODE_MAP(1))*AVOL*byam(l,i,j)/axyp(i,j))
        taijls(i,j,l,ijlt_AMPext(2))=taijls(i,j,l,ijlt_AMPext(2)) + (NACTV(i,j,l,SEAS_MODE_MAP(2))*AVOL*byam(l,i,j)/axyp(i,j))
        taijls(i,j,l,ijlt_AMPext(3))=taijls(i,j,l,ijlt_AMPext(3)) +  DIAM(i,j,l,SEAS_MODE_MAP(1)) 
        taijls(i,j,l,ijlt_AMPext(4))=taijls(i,j,l,ijlt_AMPext(4)) +  DIAM(i,j,l,SEAS_MODE_MAP(2)) 
        taijls(i,j,l,ijlt_AMPext(5))=taijls(i,j,l,ijlt_AMPext(5)) + (AERO(22) *AVOL*byam(l,i,j)/axyp(i,j)) 
        taijls(i,j,l,ijlt_AMPext(6))=taijls(i,j,l,ijlt_AMPext(6)) + (AERO(25) *AVOL*byam(l,i,j)/axyp(i,j)) 

        taijls(i,j,l,ijlt_AMPext(7))=taijls(i,j,l,ijlt_AMPext(7)) + (SUM(CCNSS(i,j,l,:,1))*AVOL*byam(l,i,j)/axyp(i,j))
        taijls(i,j,l,ijlt_AMPext(8))=taijls(i,j,l,ijlt_AMPext(8)) + (SUM(CCNSS(i,j,l,:,2))*AVOL*byam(l,i,j)/axyp(i,j))
        taijls(i,j,l,ijlt_AMPext(9))=taijls(i,j,l,ijlt_AMPext(9)) + (SUM(CCNSS(i,j,l,:,3))*AVOL*byam(l,i,j)/axyp(i,j))

#ifndef NO_HDIURN
c     Hourly Station Diagnostic -------------------------------------------------------------------------------
      idxd=(/(idd_diam+ii-1,ii=1,nmodes),(idd_numb+ii-1,ii=1,nmodes),(idd_mass+ii-1,ii=1,38),
     *        idd_so2, idd_lwp /)
      idxl=(/ (idd_ncL+l-1),  (idd_ccn+l-1), (idd_cdnc+l-1),(idd_lwc+l-1),(idd_pres+l-1)/)
      ih=jhour+1
      ihm=ih+(jdate-1)*24

      do kr=1,ndiupt
        if(i.eq.ijdd(1,kr).and.j.eq.ijdd(2,kr)) then

        if (l.eq.1) then   ! only surface diagnostic
c surface
            tmp(idd_diam:idd_diam+nmodes-1)=diam(i,j,1,1:nmodes)
            tmp(idd_numb:idd_numb+nmodes-1)= HD_NUMB(1:nmodes)    !  #/m3
            tmp(idd_mass)  =  AERO(1)               ! ug/m3;
            tmp(idd_mass+1)=  AERO(2)               ! ug/m3;
            tmp(idd_mass+2)=  AERO(3)               ! ug/m3;
            tmp(idd_mass+3)=  AERO(4)               ! ug/m3;
            tmp(idd_mass+4)=  AERO(6)               ! ug/m3;
            tmp(idd_mass+5)=  AERO(8)               ! ug/m3;
            tmp(idd_mass+6)=  AERO(9)               ! ug/m3;
            tmp(idd_mass+7)=  AERO(11)               ! ug/m3;
            tmp(idd_mass+8)=  AERO(12)               ! ug/m3;
            tmp(idd_mass+9)=  AERO(14)               ! ug/m3;
            tmp(idd_mass+10)= AERO(15)               ! ug/m3;
            tmp(idd_mass+11)= AERO(17)               ! ug/m3;
            tmp(idd_mass+12)= AERO(18)               ! ug/m3;
            tmp(idd_mass+13)= AERO(20)               ! ug/m3;
            tmp(idd_mass+14)= AERO(21)               ! ug/m3;
            tmp(idd_mass+15)= AERO(23)               ! ug/m3;
            tmp(idd_mass+16)= AERO(24)               ! ug/m3;
            tmp(idd_mass+17)= AERO(26)               ! ug/m3;
            tmp(idd_mass+18)= AERO(27)               ! ug/m3;
            tmp(idd_mass+19)= AERO(29)               ! ug/m3;
            tmp(idd_mass+20)= AERO(30)               ! ug/m3;
            tmp(idd_mass+21)= AERO(32)               ! ug/m3;
            tmp(idd_mass+22)= AERO(33)               ! ug/m3;
            tmp(idd_mass+23)= AERO(35)               ! ug/m3;
            tmp(idd_mass+24)= AERO(36)               ! ug/m3;
            tmp(idd_mass+25)= AERO(38)               ! ug/m3;
            tmp(idd_mass+26)= AERO(39)               ! ug/m3;
            tmp(idd_mass+27)= AERO(40)               ! ug/m3;
            tmp(idd_mass+28)= AERO(42)               ! ug/m3;
            tmp(idd_mass+29)= AERO(43)               ! ug/m3;
            tmp(idd_mass+30)= AERO(44)               ! ug/m3;
            tmp(idd_mass+31)= AERO(46)               ! ug/m3;
            tmp(idd_mass+32)= AERO(47)               ! ug/m3;
            tmp(idd_mass+33)= AERO(49)               ! ug/m3;
            tmp(idd_mass+34)= AERO(50)               ! ug/m3;
            tmp(idd_mass+35)= AERO(51)               ! ug/m3;
            tmp(idd_mass+36)= AERO(52)               ! ug/m3;
            tmp(idd_mass+37)= AERO(53)               ! ug/m3;
            tmp(idd_so2)  = trm(i,j,l,n_SO2)* 1.d9 / AVOL!ug/m3
            tmp(idd_lwp)  = DIURN_LWP(i,j)           ! kg/m2

            ADIURN(idxd(:),kr,ih) =ADIURN(idxd(:),kr,ih)+tmp(idxd(:))
            HDIURN(idxd(:),kr,ihm)=HDIURN(idxd(:),kr,ihm)+tmp(idxd(:))

      end if  ! level 1 only

c profile
            tmp(idd_ccn+L-1)   =  sum(nactv(i,j,l,:)) *1.e-6 ! #/cm3
            tmp(idd_ncL+L-1)   =  sum(HD_NUMB(:))     *1.e-6 ! #/cm3    
#ifdef CLD_AER_CDNC      
            tmp(idd_cdnc+L-1)  =  CDN3d(l,i,j)             
#else
            tmp(idd_cdnc+L-1)  =  0.d0
#endif    
            tmp(idd_lwc+L-1)   =  DIURN_LWC(i,j,l)         ! kg/kg(air)
            tmp(idd_pres+L-1)  =  pres         ! hPa

         
            ADIURN(idxl(:),kr,ih) =ADIURN(idxl(:),kr,ih)+tmp(idxl(:))
            HDIURN(idxl(:),kr,ihm)=HDIURN(idxl(:),kr,ihm)+tmp(idxl(:))

      end if
      end do
#endif


      ENDDO !i
      ENDDO !j
      ENDDO !l


#ifdef CACHED_SUBDD
C****
C**** Collect some high-frequency outputs
C****
      call find_groups('aijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do kr=1,subdd%ndiags
      select case (subdd%name(kr))

      case ('pn')
        call inc_subdd(subdd,kr,PN_subdd)
        PN_subdd(:,:) = 0.d0
      case ('apn')
        call inc_subdd(subdd,kr,APN_subdd)
        APN_subdd(:,:) = 0.d0
      end select
      enddo
      enddo


      call find_groups('aijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do kr=1,subdd%ndiags
      select case (subdd%name(kr))
      case ('apn_3d')
        call inc_subdd(subdd,kr,apn_3d_subdd)
      case ('pn_3d')
        call inc_subdd(subdd,kr,pn_3d_subdd)
      case ('p_3d')
        call inc_subdd(subdd,kr,p_3d_subdd)
      case ('rh_3d')
        call inc_subdd(subdd,kr,rh_3d_subdd)
      case ('ccn01_3d')
        call inc_subdd(subdd,kr,c1_3d_subdd)
      case ('ccn02_3d')
        call inc_subdd(subdd,kr,c2_3d_subdd)
      case ('ccn05_3d')
        call inc_subdd(subdd,kr,c5_3d_subdd)
      end select
      enddo
      enddo

#endif /* CACHED_SUBDD */



      RETURN
      END SUBROUTINE MATRIX_DRV
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE AMPtrdens(i,j,l,n)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the actual density per mode
!----------------------------------------------------------------------------------------------------------------------
      USE TRACER_COM
      USE AMP_AEROSOL, only : AMP_dens
      USE AERO_CONFIG, ONLY: NMODES

      IMPLICIT NONE
      Integer :: i,j,l,n,x,nAMP
      real*8, dimension(:), allocatable :: trpdens_local

      allocate(trpdens_local(ntm))
      do x=1,ntm
        trpdens_local(x)=trpdens(x)
      enddo
 
      nAMP=n-ntmAMPi+1
         if(AMP_MODES_MAP(nAMP).gt.0)
     &   AMP_dens(i,j,l,AMP_MODES_MAP(nAMP)) = 
     &   sum(trpdens_local(AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP)) * trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))) 
     & / (sum(trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))) + 1.0D-30)
         if (AMP_dens(i,j,l,AMP_MODES_MAP(nAMP)).le.0) AMP_dens(i,j,l,AMP_MODES_MAP(nAMP)) = trpdens_local(AMP_MODES_MAP(nAMP))

      deallocate(trpdens_local)

      RETURN
      END SUBROUTINE AMPtrdens
c -----------------------------------------------------------------
c -----------------------------------------------------------------
      SUBROUTINE AMPtrmass(i,j,l,n)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the actual molecular mass per mode
!----------------------------------------------------------------------------------------------------------------------
      USE TRACER_COM
      USE AMP_AEROSOL, only : AMP_TR_MM
      USE AERO_CONFIG, ONLY: NMODES

      IMPLICIT NONE
      Integer :: i,j,l,n,x,nAMP
      real*8, dimension(:), allocatable :: tr_mm_local

      allocate(tr_mm_local(ntm))
      do x=1,ntm
        tr_mm_local(x)=tr_mm(x)
      enddo

      nAMP=n-ntmAMPi+1
         if(AMP_MODES_MAP(nAMP).gt.0.and.sum(trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))).gt.0. )
     &   AMP_TR_MM(i,j,l,AMP_MODES_MAP(nAMP)) = 
     &   sum(tr_mm_local(AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP)) * trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))) 
     & / sum(trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP)))
         if (AMP_TR_MM(i,j,l,AMP_MODES_MAP(nAMP)).le.0) AMP_TR_MM(i,j,l,AMP_MODES_MAP(nAMP)) = tr_mm_local(AMP_MODES_MAP(nAMP))

      deallocate(tr_mm_local)

      RETURN
      END SUBROUTINE AMPtrmass
c -----------------------------------------------------------------
#endif GEOS5_PORT

      SUBROUTINE SPCMASSES(AERO,GAS,SPCMASS)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the total mass concentration of each model species:
!     SULF, BCAR, OCAR, DUST, SEAS, NO3, NH4. Aerosol water is not treated. 
!----------------------------------------------------------------------------------------------------------------------
      USE AERO_SETUP, ONLY: SULF_MAP, BCAR_MAP, OCAR_MAP, DUST_MAP, SEAS_MAP
      USE AERO_PARAM
      USE AERO_CONFIG
      IMPLICIT NONE
      REAL(8) :: AERO(NAEROBOX)
      REAL(8) :: GAS(NGASES)    
      REAL(8) :: SPCMASS(NMASS_SPCS+2)
      SPCMASS(1) = SUM( AERO( SULF_MAP(:) ) ) + GAS( GAS_H2SO4 )
      SPCMASS(2) = SUM( AERO( BCAR_MAP(:) ) )
      SPCMASS(3) = SUM( AERO( OCAR_MAP(:) ) )
      SPCMASS(4) = SUM( AERO( DUST_MAP(:) ) )
      SPCMASS(5) = SUM( AERO( SEAS_MAP(:) ) )
      SPCMASS(6) = AERO( MASS_NO3 )           + GAS( GAS_HNO3 )
      SPCMASS(7) = AERO( MASS_NH4 )           + GAS( GAS_NH3  )

 
      RETURN
      END SUBROUTINE SPCMASSES
     
      SUBROUTINE SIZE_PDFS(AERO,PDF1,PDF2)
      USE AERO_PARAM, ONLY: PI6, DENSP, IXXX, IYYY, ILAY
      USE AERO_CONFIG, ONLY: NMODES, NAEROBOX,NBINS
      USE AERO_SETUP, ONLY: SIG0, CONV_DPAM_TO_DGN, NUMB_MAP, MODE_NAME
      USE AMP_AEROSOL, only: DIAM
      IMPLICIT NONE

      ! Arguments.
       REAL(8), INTENT(IN) :: AERO(NAEROBOX)! aerosol conc. [ug/m^3] or [#/m^3]

      ! Local variables. 

      INTEGER :: I, N
!      INTEGER, PARAMETER :: NBINS = 30! 200    ! number of bins [1]     defined in config 
      REAL(8) :: DGRID(NBINS)              ! fixed diameter        grid [um]
      REAL(8) :: MGRID(NBINS)              ! fixed mass/particle   grid [ug/particle]
      REAL(8) :: DLOWER(NBINS)             ! lower boundary fixed diameter grid [um]
      REAL(8) :: DUPPER(NBINS)             ! upper boundary fixed diameter grid [um]
      REAL(8) :: NTOT(NMODES)              ! number concentration for each mode [#/m^3]
      REAL(8) :: PDF(NBINS,2,NMODES)       ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF1(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF2(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: DNDLOGD(NMODES)           ! dN/dlog10(Dp) [ #/m^3]
      REAL(8) :: DMDLOGD(NMODES)           ! dM/dlog10(Dp) [ug/m^3]
      REAL(8) :: RDMIN                     ! reciprocal of DMIN to optimize coagulation [1/um]
      REAL(8) :: RDLOGDSC                  ! reciprocal of log10 of the grid spacing [1]
      REAL(8) :: SCALE, F, SUM1, SUM2      ! scratch variables 
      REAL(8) :: DMINL, DMAXL, DG          ! diameters [um]  
      REAL(8) :: FLN                       ! function for lognormal distribution [1]  
      REAL(8), PARAMETER :: DMIN =  0.001D+00   ! smallest particle diameter of the discrete grid [um]
      REAL(8), PARAMETER :: DMAX = 20.000D+00   ! largest  particle diameter of the discrete grid [um]


        DMAXL = DMAX
        DMINL = DMIN

      SCALE    = ( DMAXL / DMINL )**(1.0D+00/REAL(NBINS-1))
      RDLOGDSC = 1.0D+00 / LOG10( SCALE )
      RDMIN    = 1.0D+00 / DMINL
      DO I=1, NBINS
        DGRID(I)  = DMINL * SCALE**(I-1)                  ! [um]
        DLOWER(I) = DGRID(I) / SCALE**0.5D+00             ! [um]
        DUPPER(I) = DGRID(I) * SCALE**0.5D+00             ! [um]
        MGRID(I)  = 1.0D-06 * DENSP * PI6 * DGRID(I)**3   ! [ug/particle]
        DO N=1, NMODES
          DG = 1.0D+06 * DIAM(IXXX,IYYY,ILAY,N) * CONV_DPAM_TO_DGN(N)   ! convert [m] to [um] and Dbar to Dg
          NTOT(N) = AERO( NUMB_MAP(N) )
          F = NTOT(N) * FLN( DGRID(I), DG, SIG0(N) )
          PDF(I,1,N) = F * ( DUPPER(I) - DLOWER(I) )
          PDF(I,2,N) = PDF(I,1,N) * MGRID(I)
          DNDLOGD(N) = PDF(I,1,N) * RDLOGDSC * 1.0D-06           ! convert from [#/m^3] to [#/cm^3]
          DNDLOGD(N) = MAX( DNDLOGD(N), 1.0D-30 )
          DMDLOGD(N) = PDF(I,2,N) * RDLOGDSC                     ! [ug/m^3]
          DMDLOGD(N) = MAX( DMDLOGD(N), 1.0D-30 )
        ENDDO
c        WRITE(IUNIT,91) I, DGRID(I), DNDLOGD(:)
c        WRITE(JUNIT,91) I, DGRID(I), DMDLOGD(:)
      ENDDO

        PDF1(:) = 0.0D+00
        PDF2(:) = 0.0D+00
      DO N=1, NMODES
        DO I=1, NBINS
          PDF1(I) = PDF1(I)  + PDF(I,1,N)
          PDF2(I) = PDF2(I)  + PDF(I,2,N)
          SUM1 = SUM1 + PDF(I,1,N)
          SUM2 = SUM2 + PDF(I,2,N)
        ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE SIZE_PDFS


      REAL(8) FUNCTION FLN(X,XG,SIGMAG)
      REAL(8) :: X      ! particle radius or diameter [any units]
      REAL(8) :: XG     ! geometric mean radius or diameter [any units]
      REAL(8) :: SIGMAG ! geometric standard deviation [monodisperse = 1.0]
      REAL(8), PARAMETER :: SQRTTWOPI = 2.506628275D+00
      FLN = EXP(-0.5D+00*(LOG(X/XG)/LOG(SIGMAG))**2) / (X*LOG(SIGMAG)*SQRTTWOPI)
      RETURN
      END FUNCTION FLN

#ifndef GEOS5_PORT
      subroutine alloc_tracer_amp_com(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth Susanne Bauer
!@ver  1.0
      use domain_decomp_atm, only : dist_grid, get
      use model_com, only     : im,lm
      use amp_aerosol
      use aero_config, only   : nmodes

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

! I,J,L
      allocate(  AQsulfRATE(I_0H:I_1H,J_0H:J_1H,LM)   )
! other dimensions
      allocate(  DIAM(I_0H:I_1H,J_0H:J_1H,LM,nmodes)  )
      allocate(  AMP_TR_MM(I_0H:I_1H,J_0H:J_1H,LM,nmodes)  )
      allocate(  AMP_dens(I_0H:I_1H,J_0H:J_1H,LM,nmodes)  )
      allocate(  NACTV(I_0H:I_1H,J_0H:J_1H,LM,nmodes) )
      allocate(  CCNSS(I_0H:I_1H,J_0H:J_1H,LM,nmodes,3) )
      allocate(  VDDEP_AERO(I_0H:I_1H,J_0H:J_1H,nmodes,2))
      allocate(  NUMB_SS(I_0H:I_1H,J_0H:J_1H,LM,2))
#ifndef NO_HDIURN
      allocate(  DIURN_LWP(I_0H:I_1H,J_0H:J_1H))     ! lwp hourly output
      allocate(  DIURN_LWC(I_0H:I_1H,J_0H:J_1H,LM))  ! lwc hourly output
#endif

      NACTV   = 1.0D-30
      DIAM    = 1.0D-30
      NUMB_SS = 1.0D-30       
      return
      end subroutine alloc_tracer_amp_com
#endif    
