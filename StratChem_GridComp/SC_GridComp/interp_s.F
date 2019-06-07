      SUBROUTINE interp_s(k,sza,o3column,s,jo2,gcSC)

! ----------------------------------------------------------------------------
! NAME:
!   interp_s
!
! PURPOSE:
!   Interpolate S values for each wavelength in table to specified O3
!   column and zenith angle
!
! INPUTS:
!   k         Current layer number
!   szaRad    Solar zenith angle [radians]
!   o3column  Overhead o3 column value [cm^{-2}]
!   gcSC      The StratChem grid component, which contains
!     sza_tab Solar zenith angle table
!     o3_tab  Overhead O3 values table
!     sdat    Radiative source function 
!     o2jdat  Table of J(O2) values
!
! OUTPUTS:
!   s         S value for each wavelength at current k, interpolated to
!               the given o3column and sza
!   jo2       J(O2) values interpolated as above
!
! 
! PROCEDURE:
!   Bi-linear interpolation, for sza > 94 s=0, for O3 out of range use min/max
!
! MODIFICATION HISTORY: 
!   25 Aug 1993  Kawa
!   10 Jul 1996  Kawa For 28 levels and to handle J(O2) separately
!   11 May 2012  Nielsen Accomodation for GEOS-5 FV cubed release
! ----------------------------------------------------------------------------
       
      USE SC_GridCompMod
      IMPLICIT NONE

      TYPE(SC_GridComp), INTENT(IN) :: gcSC   ! Grid Component

      INTEGER, INTENT(IN) :: k
      REAL, INTENT(IN) :: sza, o3column 
      REAL, INTENT(OUT) :: s(gcSC%nlam), jo2

      INTEGER :: ijj, ik, ikk, ikkm, il, is
      REAL :: omt, omu, t, u
      REAL, PARAMETER :: PI = 3.14159265

! For each input solar zenith angle, find the first element of gcSC%sza_tab that 
! is greater.  Use this element and previous one to determine the interpolated value.
! -----------------------------------------------------------------------------------
      DO is = 1,gcSC%nsza
         ijj = is 
         IF(gcSC%sza_tab(is) > sza) EXIT 
      END DO
      
! Zenith angle test       
! -----------------
      IF(sza > gcSC%sza_tab(gcSC%nsza)) THEN

! Cell is dark, set s and jo2=0        
! -----------------------------
         s(1:gcSC%nlam) = 0.
         jo2 = 0.

      ELSE  

! Cell is illuminated     
! -------------------
         t = (sza-gcSC%sza_tab(ijj-1))/(gcSC%sza_tab(ijj)-gcSC%sza_tab(ijj-1))
         omt = 1.-t
         
! For each overhead O3 column, find the first element in gcSC%o3_tab that is
! greater. Use this element and previous one to determine the interpolated value.
! -------------------------------------------------------------------------------
         DO is = 1,gcSC%numo3
            ikk = is 
            IF(gcSC%o3_tab(is,k) > o3column) EXIT
         END DO

         ikkm = ikk-1 
         IF(ikk > 1 .AND. o3column <= gcSC%o3_tab(gcSC%numo3,k)) THEN
            u = (o3column-gcSC%o3_tab(ikkm,k))/(gcSC%o3_tab(ikk,k)-gcSC%o3_tab(ikkm,k))
            omu = 1.-u

! Do bilinear interpolation for each wavelength.
! ----------------------------------------------
            DO il = 1,gcSC%nlam       
               s(il) = omt*omu*gcSC%sdat(ijj-1,ikkm,k,il)+t*omu*gcSC%sdat(ijj,ikkm,k,il)+
     &                 t*u*gcSC%sdat(ijj,ikk,k,il)+omt*u*gcSC%sdat(ijj-1,ikk,k,il)
            END DO
            jo2 = omt*omu*gcSC%o2jdat(ijj-1,ikkm,k)+t*omu*gcSC%o2jdat(ijj,ikkm,k)+
     &           +t*u*gcSC%o2jdat(ijj,ikk,k)+omt*u*gcSC%o2jdat(ijj-1,ikk,k)
    
! Extrapolate ahead of table
! --------------------------
         ELSE IF (ikk == 1) THEN
            DO il = 1,gcSC%nlam
               s(il) = omt*gcSC%sdat(ijj-1,1,k,il)+t*gcSC%sdat(ijj,1,k,il)
            END DO
            jo2 = omt*gcSC%o2jdat(ijj-1,1,k)+t*gcSC%o2jdat(ijj,1,k)

! Extrapolate beyond table
! ------------------------
         ELSE
            DO il = 1,gcSC%nlam
               s(il) = omt*gcSC%sdat(ijj-1,gcSC%numo3,k,il)+t*gcSC%sdat(ijj,gcSC%numo3,k,il)
            END DO 
            jo2 = omt*gcSC%o2jdat(ijj-1,gcSC%numo3,k)+t*gcSC%o2jdat(ijj,gcSC%numo3,k)
         END IF  
      END IF
      
      RETURN
      END SUBROUTINE interp_s
