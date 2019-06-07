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
      subroutine calcShipEmission (emiss_o3, thisRecord, curRecord, &
                     latdeg, jno2val,       &
     &               emissionArray, emiss_ozone, o3_index, &
                     i1, i2, ju1, j2, k1, k2, ju1_gl, j2_gl, num_emiss)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none
!
! !INPUT PARAMETERS:
      integer                , intent(in) :: i1, i2, ju1, j2, k1, k2
      integer                , intent(in) :: ju1_gl, j2_gl
      integer                , intent(in) :: o3_index, curRecord
      integer                , intent(in) :: num_emiss
      real*8                 , intent(in) :: latdeg (i1:i2, ju1:j2)
      real*8                 , intent(in) :: jno2val(i1:i2, ju1:j2)
!
! !INPUT/OUTPUT PARAMETERS:
      integer                , intent(inOut) :: thisRecord
      real*8                 , intent(inOut) :: emiss_o3(i1:i2, ju1:j2)
      real*8                 , intent(inOut) :: emiss_ozone(i1:i2, ju1:j2)
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

      ! Check if new record (month or day).
      if (thisRecord /= curRecord) then
          emiss_ozone(:,:) = emissionArray(o3_index)%pArray3D(:,:,1)
          thisRecord       = curRecord
      end if
        
      emiss_o3(:,:) = 0.0d0

      do j=ju1,j2
         do i=i1,i2

            if ((latdeg(i,j).lt.60.0d0) .and. (latdeg(i,j).gt.-60.0d0)) then
               if (jno2val(i,j).le.0.0095d0) then
                  tempVar       = jno2val(i,j)/0.0095d0
                  emiss_o3(i,j) = 10.0d0 * emiss_ozone(i,j) * tempVar*tempVar
               else
                  emiss_o3(i,j) = emiss_ozone(i,j) * 10.0d0
               endif
            endif

         enddo 
      enddo

      emissionArray(o3_index)%pArray3D(:,:,1) = emiss_o3(:,:)

      return

      end subroutine calcShipEmission
!EOC
!-------------------------------------------------------------------------
      end module GmiShipEmission_mod

