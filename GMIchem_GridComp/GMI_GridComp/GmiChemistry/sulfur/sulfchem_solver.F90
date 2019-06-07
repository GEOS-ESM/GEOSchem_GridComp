!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   Xiaohong Liu
!
! FILE
!   sulfchem_solver.F
!
! ROUTINES
!   Do_Sulf_Solver
!
! HISTORY
!  - December 8, 2005 * Bigyani Das
!     Added two variables do_aerocom, do_dust_emiss to the argument list
!     of Do_Sulf_Solver and Do_Sulfchem
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Sulf_Solver
!
! DESCRIPTION
!   This is the main control routine for the chemistry solver, "Sulfchem".
!
! ARGUMENTS
!   time       : chemistry time step (s)
!   nymd       : current  year/month/day (YYYYMMDD)
!   tempc_nb   : temperature (K)
!   presc_nb   : atmospheric pressure at the center of each grid box (mb)
!   const_nb   : species concentration, known at zone centers
!                (molecules/cm^3 for gas-species and g/g for aerosols)
!   humidity   : specific humidity (g/kg)
!   sulf_emiss : array of sulfur emissions (molecules/cm^3/s)
!   qjgmi      : photolysis rate constants (s^-1)
!   qkgmi      : thermal    rate constants (units vary)
!   clfrac     : total cloud fraction [0 - 1]
!   lwccol     : in-cloud liquid water in each grid box (g/g)
!   relhum     : relative humidity [0 - 1]
!
!-----------------------------------------------------------------------------

      subroutine Do_Sulf_Solver  &
     &  (time, nymd, nhms, tempc_nb, presc_nb, const_nb, humidity,  &
     &   sulf_emiss, qjgmi, qkgmi, phot_opt, clfrac, lwccol, relhum,  &
     &   aqua_infile_name, num_time_steps, londeg, latdeg, mass,  &
     &   pr_sulf_src, do_aerocom, do_dust_emiss,  &
     &   dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3, &
     &   pr_diag, loc_proc, ilong, itloop, &
     &   i1, i2, ju1, j2, k1, k2, i1_gl, i2_gl, ju1_gl, j2_gl, &
     &   num_qjs, num_qks, num_species)

      use GmiTimeControl_mod, only : GmiSplitDateTime
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

#     include "setkin_par.h"
#     include "sulfchem.h"

!     ----------------------
!     GMIMOD-related values.
!     ----------------------

      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, phot_opt
      integer, intent(in) :: ilong, itloop
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: num_qjs, num_qks, num_species
      integer :: nymd
      integer :: nhms
      integer :: num_time_steps

      character*(*) :: aqua_infile_name

      real*8  :: time
      real*8  :: tempc_nb  (i1:i2, ju1:j2, k1:k2)
      real*8  :: presc_nb  (i1:i2, ju1:j2, k1:k2)
      real*8  :: const_nb  (i1:i2, ju1:j2, k1:k2, num_species)
      real*8  :: humidity  (i1:i2, ju1:j2, k1:k2)
      real*8  :: sulf_emiss(i1:i2, ju1:j2, k1:k2, num_species)
      real*8  :: clfrac    (i1:i2, ju1:j2, k1:k2)
      real*8  :: lwccol    (i1:i2, ju1:j2, k1:k2)
      real*8  :: relhum    (i1:i2, ju1:j2, k1:k2)
      real*8  :: londeg    (i1_gl:i2_gl)
      real*8  :: latdeg    (ju1_gl:j2_gl)
      real*8  :: mass      (i1:i2, ju1:j2, k1:k2)
      logical :: pr_sulf_src
      logical :: do_aerocom
      logical :: do_dust_emiss
      real*8  :: dms_oh    (i1:i2, ju1:j2, k1:k2)
      real*8  :: dms_no3   (i1:i2, ju1:j2, k1:k2)
      real*8  :: so2_oh    (i1:i2, ju1:j2, k1:k2)
      real*8  :: so2_h2o2  (i1:i2, ju1:j2, k1:k2)
      real*8  :: so2_o3    (i1:i2, ju1:j2, k1:k2)

                             ! photolysis rate constants (s^-1)
      type (t_GmiArrayBundle), intent(inOut) :: qjgmi(num_qjs)
                             ! thermal rate constants (units vary)
      type (t_GmiArrayBundle), intent(in) :: qkgmi(num_qks)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, ik
      integer :: idumday, idumyear
      integer :: im, it

      real*8 :: oh_saved (i1:i2, ju1:j2, k1:k2)
      real*8 :: ho2_saved(i1:i2, ju1:j2, k1:k2)
      real*8 :: no3_saved(i1:i2, ju1:j2, k1:k2)

      real*8  :: hox_qj_scaler(i1:i2, ju1:j2)
      real*8  :: no3_scaler   (i1:i2, ju1:j2)

      real*8, allocatable :: loc_qjgmi(:,:,:,:)
      real*8, allocatable :: loc_qkgmi(:,:,:,:)


!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_Sulf_Solver called by ', loc_proc
      end if


      hox_qj_scaler = 0.0d0
      no3_scaler    = 0.0d0

     if (phot_opt == 7) then

!        ====================
         call GmiSplitDateTime  &
!        ====================
     &     (nymd, idumyear, im, idumday)


!        Interpolation of qj (91x5) to the GMIMOD pressure levels.

         do ik = k1, k2
           do ij = ju1, j2
             do il = i1, i2

               k = Int (2.0d0  &
     &           + (presc_nb(il,ij,ik) * 1.0d-3 - 0.1d0) / 0.3d0)
               k  = Min (k, 4)

               qjgmi(1)%pArray3D(il,ij,ik) =  &
     &               qjh2o2_2d(ij,k,im)  &
     &               + (qjh2o2_2d(ij,k+1,im) - qjh2o2_2d(ij,k,im))  &
     &               * (presc_nb(il,ij,ik) * 1.0d-3 - ppr(k))  &
     &               / (ppr(k+1) - ppr(k))

             end do
           end do
         end do
      endif


!     Calculate the scaling factor for OH, HO2, NO3, Qj diurnal variations

!     =====================
      call Scale_HOx_NO3_Qj  &
!     =====================
     &  (nymd, nhms, londeg, latdeg, hox_qj_scaler, no3_scaler, &
     &   i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl, ilong)

      oh_saved  = const_nb(i1:i2,ju1:j2,k1:k2,IOH)
      ho2_saved = const_nb(i1:i2,ju1:j2,k1:k2,IHO2)
      no3_saved = const_nb(i1:i2,ju1:j2,k1:k2,INO3)

      do ik = k1, k2
        do ij = ju1, j2
          do il = i1, i2

            const_nb(il,ij,ik,IOH) =  &
     &        const_nb(il,ij,ik,IOH) * hox_qj_scaler(il,ij)
            const_nb(il,ij,ik,IHO2) =  &
     &        const_nb(il,ij,ik,IHO2) * hox_qj_scaler(il,ij)
            qjgmi(1)%pArray3D(il,ij,ik) =  &
     &        qjgmi(1)%pArray3D(il,ij,ik) * hox_qj_scaler(il,ij)
            const_nb(il,ij,ik,INO3) =  &
     &        const_nb(il,ij,ik,INO3) * no3_scaler(il,ij)

          end do
        end do
      end do

       allocate(loc_qjgmi(i1:i2, ju1:j2, k1:k2, num_qjs))
       do it = 1, num_qjs
          loc_qjgmi(:,:,:,it) = qjgmi(it)%pArray3D(:,:,:)
       end do

       allocate(loc_qkgmi(i1:i2, ju1:j2, k1:k2, num_qks))
       do ic = 1, num_qks
          loc_qkgmi(:,:,:,it) = qkgmi(it)%pArray3D(:,:,:)
       end do

!     ================
      call Do_Sulfchem  &
!     ================
     &  (itloop, time, tempc_nb, presc_nb, const_nb, humidity,  &
     &   sulf_emiss, loc_qjgmi, loc_qkgmi, clfrac, lwccol, relhum,  &
     &   aqua_infile_name, pr_diag, loc_proc, num_time_steps,  &
     &   mass, pr_sulf_src, do_aerocom, do_dust_emiss,  &
     &   dms_oh, dms_no3, so2_oh, so2_h2o2, so2_o3)

       do it = 1, num_qjs
          qjgmi(it)%pArray3D(:,:,:) = loc_qjgmi(:,:,:,it)
       end do
      deallocate(loc_qjgmi, loc_qkgmi)

      const_nb(i1:i2,ju1:j2,k1:k2,IOH)  = oh_saved
      const_nb(i1:i2,ju1:j2,k1:k2,IHO2) = ho2_saved
      const_nb(i1:i2,ju1:j2,k1:k2,INO3) = no3_saved

      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Scale_HOx_NO3_Qj
!
! CODE DEVELOPER
!   Xiaohong Liu
!
! DESCRIPTION
!   This routine computes the scaling factor for HOx-Qj based on
!
!   OH(instantaneous) = OH(max) * cos(theta)
!
!   where cos(theta) is cosine of solar zenith angle.
!
!   thus, OH(monthly_ave) * 24 (hr) = OH(max) * Sum (cos(theta) * dt)
!                                     ( t (hr) = triz -> tset)
!
!   OH(max) = OH(monthly_ave) * 24 (hr) / (Sum (cos(theta) * dt) )
!
!   so, hox_qj_scaler = 24 / (Sum (cos(theta) * dt) ) * cos(theta)
!
!   and, computes the scaling factor for NO3 assuming NO3 even distributed
!        over night.
!   so, no3_scaler = 24 / (24-(tset-triz))
!
!-----------------------------------------------------------------------------

      subroutine Scale_HOx_NO3_Qj  &
     &  (nymd, nhms, londeg, latdeg, hox_qj_scaler, no3_scaler, &
     &   i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, j2_gl, ilong)

      use GmiTimeControl_mod, only : GmiSplitDateTime
      use GmiTimeControl_mod, only : GetDaysFromJanuary1, GetSecondsFromJanuary1

      implicit none

#     include "gmi_time_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, ilong

      integer :: nymd
      integer :: nhms

      real*8  :: londeg (i1_gl:i2_gl)
      real*8  :: latdeg (ju1_gl:j2_gl)
      real*8  :: hox_qj_scaler(i1:i2, ju1:j2)
      real*8  :: no3_scaler   (i1:i2, ju1:j2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: il, ij, it
      integer :: idumyear, idummin, idumsec
      integer :: im, id, ih
      integer :: nday_jan1
      integer :: nsec_jan1
      integer :: nt

      real*8 :: tloc
      real*8 :: stime_loc
      real*8 :: days
      real*8 :: decl     ! solar declination
      real*8 :: rdistsq  ! solar distance correction
      real*8 :: time
      real*8 :: triz(ju1:j2)
      real*8 :: tset(ju1:j2)
      real*8 :: tot_cossza(ju1:j2)
      real*8, save :: lon_dummy(1) = 0.0d0
      real*8 :: cossza(1)
      real*8 :: csza(i1:i2)

!     ----------------
!     Begin execution.
!     ----------------

      csza(:) = 0.0d0

!     ====================
      call GmiSplitDateTime  &
!     ====================
     &  (nymd, idumyear, im, id)

!     ====================
      call GmiSplitDateTime  &
!     ====================
     &  (nhms, ih, idummin, idumsec)

!     =================
      call GetDaysFromJanuary1  &
!     =================
     &  (nday_jan1, nymd)


!     Calculate the sun rize and sun set time (local time: hr)

      do ij = ju1, j2

!       ==============
        call SunRizSet  &
!       ==============
     &    (im, id, latdeg(ij), triz(ij), tset(ij))

      end do


      nt = 100       ! divide one-day (24 hr) into 100 time intervals
      tot_cossza(:) = 0.0d0

!     =============
      do it = 0, nt
!     =============

        tloc = it * 24.0d0/nt      ! hr
        stime_loc = nday_jan1 * SECPDY + tloc * SECPHR

        days = stime_loc / SECPDY

!       ============
        call Solrpos  &
!       ============
     &    (days, decl, rdistsq)


!       ===============
        do ij = ju1, j2
!       ===============

!         ---------------------------------------------
!         Calculate solar zenith angle at this latitude
!         ---------------------------------------------

!         ===========
          call Solrza  &
!         ===========
     &      (days, decl, latdeg(ij), lon_dummy, 1, cossza)


          if(tloc > triz(ij) .and. tloc < tset(ij)) then

            tot_cossza(ij) = tot_cossza(ij) +  &
     &                       cossza(1) * 24.0d0/nt

          end if

!       ======
        end do
!       ======

!     ======
      end do
!     ======


!     =================
      call GetSecondsFromJanuary1  &
!     =================
     &  (nsec_jan1, nymd, nhms)

      time = nsec_jan1
      days = time / SECPDY

!     ============
      call Solrpos  &
!     ============
     &  (days, decl, rdistsq)


!     ===============
      do ij = ju1, j2
!     ===============

!       -----------------------------------------------------
!       Calculate solar zenith angle at this latitude for all
!       longitudes.
!       -----------------------------------------------------
!       ===========
        call Solrza  &
!       ===========
     &    (days, decl, latdeg(ij), londeg(i1), ilong, csza)


        do il = i1, i2

          tloc = ih + londeg(il) / 360.d0 * 24.0d0
          if(tloc > 24.0d0) tloc = tloc - 24.0d0

          if(tloc > triz(ij) .and. tloc < tset(ij)) then
            hox_qj_scaler(il,ij) = 24.0d0 /  &
     &                             tot_cossza(ij) * csza(il)
          else
            hox_qj_scaler(il,ij) = 1.0d-5
          end if

          if(tloc < triz(ij) .or. tloc > tset(ij)) then
            no3_scaler(il,ij) = 24.0d0 /  &
     &                          (24.0d0 - (tset(ij) - triz(ij)))
          else
            no3_scaler(il,ij) = 1.0d-5
          end if

        end do

!     ======
      end do
!     ======


      return

      end


!LL Subroutine: SunRizSet.f--------------------------------------------
!LL Purpose: Compute the local time of sunrise and sunset given
!LL          latitude and date.
!LL
!LL Author: Wenwei Pan
!LL         Center for Meteorology and Oceanography
!LL         Center for Global Change Science
!LL         Department of Earth, Atmospheric and Planetary Sciences
!LL
!LL Version: 1.0
!LL Date: April 3, 1997
!LL
!LL Tested under compiler: f77 V4.0
!LL Tested under OS version: SunOS 5.5.1 Generic sun4u sparc SUNW,Ultra-1
!LL
!LL Code modification history:
!LL Version Date
!LL---------------------------------------------------------------------

      subroutine SunRizSet(month,iday,rlat,triz,tset)

      implicit none

! Input variables.
      integer  &
     &     month,  &
     &     iday
      real*8 ::  &
     &     rlat           !Latitude in degrees.

! Output variables.
      real*8 :: triz,tset      !local time of sunrise and sunset(hr)

! Local variables.
      real*8 ::  &
     &     fac,  &
     &     Pi,  &
     &     x,  &
     &     xlat,            & !Latitude in radians.
     &     phi,cosphi,sinphi,  &
     &     tphi,costph,sintph,  &
     &     ttphi,costtp,sinttp,  &
     &     delta,           & !Solar inclination angle in radians.
     &     h              !Hour angle at sunrise and sunset in radians.

!*----------------------------------------------------------------------

      FAC = .01745d0
      PI = 3.1415926d0
      X = FLOAT((MONTH-1)*30+IDAY+MONTH/2)
      XLAT = RLAT*FAC
      PHI = 2.0d0*3.14159d0*X/365.0d0
      TPHI = PHI*2.0d0
      TTPHI = PHI*3.0d0
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
      COSTPH = COS(TPHI)
      SINTPH = SIN(TPHI)
      COSTTP = COS(TTPHI)
      SINTTP = SIN(TTPHI)
      DELTA = .006918d0 - 0.399912d0*COSPHI + 0.070257d0*SINPHI -  &
     &  0.006758d0*COSTPH + 0.000907d0*SINTPH - 0.002697d0*COSTTP +  &
     &  0.00148d0*SINTTP
      IF(DELTA.GE.0.0d0)THEN                   !SUMMER
         IF(XLAT.GE.(PI/2.0d0-DELTA))THEN      !NO SUNSET
            H=PI
         ELSEIF(XLAT.LE.(DELTA-PI/2.0d0))THEN  !NO SUNRISE
            H=0.0d0
         ELSE
            H=ACOS(-TAN(XLAT)*TAN(DELTA))
         ENDIF
      ELSE
         IF(XLAT.GE.(PI/2.0d0+DELTA))THEN      !NO SUNRISE
            H=0.0d0
         ELSEIF(XLAT.LE.(-DELTA-PI/2.0d0))THEN !NO SUNSET
            H=PI
         ELSE
            H=ACOS(-TAN(XLAT)*TAN(DELTA))
         ENDIF
      ENDIF
      Triz=12.0d0*(-H+PI)/PI
      Tset=12.0d0*(H+PI)/PI

      return

      end

