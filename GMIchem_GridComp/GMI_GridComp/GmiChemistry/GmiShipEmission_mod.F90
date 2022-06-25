!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiShipEmission_mod
!
! !INTERFACE:
!
      module GmiShipEmission_mod
!
      implicit none
!
        INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public   :: calcShipEmission
!EOP
!=============================================================================
      CONTAINS
!=============================================================================
!BOP
!
! !IROUTINE: calcShipEmission
!
! !INTERFACE:
!
      subroutine calcShipEmission (emiss_o3, emiss_hno3, &
                     latdeg, jno2val,       &
                     emissionArray, emiss_ship_no, ship_o3_index, ship_hno3_index, &
                     cellArea, i1, i2, ju1, j2, num_emiss)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none
!
! !INPUT PARAMETERS:
      integer                , intent(in) :: i1, i2, ju1, j2
      integer                , intent(in) :: ship_o3_index
      integer                , intent(in) :: ship_hno3_index
      REAL(KIND=DBL), POINTER, intent(in) :: cellArea(:,:)
      integer                , intent(in) :: num_emiss
      real*8                 , intent(in) :: latdeg (i1:i2, ju1:j2)
      real                   , intent(in) :: jno2val(i1:i2, ju1:j2)
      real                   , intent(in) :: emiss_ship_no(i1:i2, ju1:j2)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8                 , intent(inOut) :: emiss_o3(i1:i2, ju1:j2)
      real*8                 , intent(inOut) :: emiss_hno3(i1:i2, ju1:j2)
      type (t_GmiArrayBundle), intent(inOut) :: emissionArray(num_emiss)
!
! !DESCRIPTION: 
! This routine calculates the ship emissions.
! Only able to preprocess ship emissions if doing surface emissions inside 
! the chemistry solver (smvgear). Otherwise, NO is directly emitted to code, 
! which causes excessive O3 production.
! Preprocess ship NO emissions as:
!     \begin{itemize}
!     \item E(hno3) = E(no)
!     \item E(o3) = $E(no) \times 10 \times (jno2/0.0095)^2$
!     \end{itemize}
! Following loosely 
!     Chen et al., An investigation of the chemistry of ship emission plumes
!     during ITCT 2002, J. Geophys. Res., Vol. 110, No. D10, 
!     D10S90 10.1029/2004JD005236, 20 May 2005.
!
! !LOCAL VARIABLES:
      integer        :: i,j
      real*8         :: tempVar
!
! !REVISION HISTORY:
!  Initial code. Bryan Duncan 2-28-06
!
!EOP
!-------------------------------------------------------------------------
!BOC

      emiss_o3(:,:) = 0.0d0

      do j=ju1,j2
         do i=i1,i2

            if ((latdeg(i,j).lt.60.0d0) .and. (latdeg(i,j).gt.-60.0d0)) then
               if (jno2val(i,j).le.0.0095d0) then
                  tempVar       = jno2val(i,j)/0.0095d0
                  emiss_o3(i,j) = 10.0d0 * emiss_ship_no(i,j) * tempVar*tempVar
               else
                  emiss_o3(i,j) = 10.0d0 * emiss_ship_no(i,j)
               endif
            endif

         enddo 
      enddo

      ! Ship NO values are kg/m2/s 
      ! Both the Emissions diagnostic, and the Emissions Array are in units kg/s

      emiss_o3(:,:) = emiss_o3(:,:) * cellArea(:,:)

      emissionArray(ship_o3_index)%pArray3D(:,:,1 ) = emiss_o3(:,:)
      emissionArray(ship_o3_index)%pArray3D(:,:,2:) = 0.0


      ! Fill HNO3 array entry
      emissionArray(ship_hno3_index)%pArray3D(:,:,1 ) = emiss_ship_no(:,:) * 2.1 * cellArea
      emissionArray(ship_hno3_index)%pArray3D(:,:,2:) = 0.0

      ! Fill the emissions diagnostic
      emiss_hno3(:,:) =                                 emiss_ship_no(:,:) * 2.1 * cellArea

      return

      end subroutine calcShipEmission
!EOC
!-------------------------------------------------------------------------
      end module GmiShipEmission_mod

