#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1              !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE
module GmiSwapSpeciesBundlesMod
!
! !USES:
   USE ESMF
   USE MAPL

      USE GmiStringManipulation_mod, ONLY : stringLowerCase
      USE GmiArrayBundlePointer_mod, ONLY : t_GmiArrayBundle
      USE Species_ArrayMod
!
      implicit none
!
      private
      public  :: SwapSpeciesBundles
      public  :: speciesReg_for_CCM

#include "GmiParameters.h"
#include "gmi_phys_constants.h"
!
! !DESCRIPTION:
! Provides tools to transfer species concentrations between GEOS-5
! and GMI.
!
! !AUTHOR:
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  SwapSpeciesBundles
!
! !INTERFACE:

      SUBROUTINE SwapSpeciesBundles(direction, gmiBundle, geos5BundleGG, geos5BundleXX, &
                     specHum, mapSpecies, lchemvar, do_synoz, numSpecies, rc)

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER, INTENT(IN) :: numSpecies
      INTEGER, INTENT(IN) :: direction       ! 1=To GMI  -1=FromGMI
#define GEOS_TO_GMI 1
#define GMI_TO_GEOS -1
      INTEGER, INTENT(IN) :: mapSpecies(numSpecies) ! map GMI species indices to CCM indices
      LOGICAL, INTENT(IN) :: do_synoz
      character (len=*), intent(in) :: lchemvar(numSpecies) ! GMI list of species
!
! !OUTPUT PARAMETERS:
      INTEGER, INTENT(OUT) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
      REAL                    :: specHum(:,:,:)
      type (t_GmiArrayBundle) :: gmiBundle(:)
      type(Species_Array)     :: geos5BundleGG(:)
      type(Species_Array)     :: geos5BundleXX(:)
!
! !DESCRIPTION:
!  Copy chemistry bundle to GMI's species concentration bundle and 
!  vice--versa. Use species names for determining the indicies, because the
!  ordering in the AGCM is not the same as in GMI. The vertical ordering
!  in GMI is surface-to-lid.
!
!  When synthetic ozone is not active, set it's concentration to zero, and
!  assume that it does not exist in the ACGM's bundle. [If it disappears
!  from a future GMI strat/trop mechanism it can be eliminated here, too.]
!
! !LOCAL VARIABLES:
      INTEGER :: STATUS
      INTEGER :: i, ic
      INTEGER :: k, km
      INTEGER :: nqGG, nqXX
      LOGICAL :: found
      LOGICAL :: FeedBack_QV    ! TRUE for CCM, FALSE for CTM
      CHARACTER(LEN=255) :: speciesName
      CHARACTER(LEN=255) :: MyName
      CHARACTER(LEN=ESMF_MAXSTR) :: Iam
!
!EOP
!---------------------------------------------------------------------------
!BOC
      rc=0
      MyName="GMI_SwapSpeciesBundles"
      Iam="SwapSpeciesBundles"

      CALL CheckConfig(FeedBack_QV, __RC__)

      IF ( direction /= GEOS_TO_GMI .AND. direction /= GMI_TO_GEOS ) THEN
         rc=1
         PRINT *,TRIM(MyName),": Invalid direction specified for swapping bundles."
         RETURN
      END IF

      km = size(specHum, 3)

      nqGG = size(geos5BundleGG)
      nqXX = size(geos5BundleXX)

      DO ic=1,numSpecies
         found=.FALSE.
         speciesName=TRIM(lchemvar(ic))
    
         IF (TRIM(speciesName) == "H2O") THEN
    
            SELECT CASE (direction)
            CASE ( GEOS_TO_GMI )
               gmiBundle(ic)%pArray3D(:,:,km:1:-1) = specHum(:,:,1:km)*(MWTAIR/MWTH2O)
            CASE ( GMI_TO_GEOS )
              IF ( FeedBack_QV ) THEN
               specHum(:,:,1:km) = gmiBundle(ic)%pArray3D(:,:,km:1:-1)*(MWTH2O/MWTAIR)
              END IF
            END SELECT
    
            found=.TRUE.
    
         ELSE IF(TRIM(speciesName) == "SYNOZ" .OR. do_synoz .eqv. .TRUE.) THEN
    
            found=.FALSE.
    
         ELSE
            i = mapSpecies(ic)
    
            SELECT CASE (direction)
            CASE ( GEOS_TO_GMI )
               IF ( i .LE. nqGG ) THEN
                 gmiBundle(ic)%pArray3D(:,:,km:1:-1) = geos5BundleGG(i     )%data3d(:,:,1:km)
               ELSE
                 gmiBundle(ic)%pArray3D(:,:,km:1:-1) = geos5BundleXX(i-nqGG)%data3d(:,:,1:km)
               END IF
! ########
!                IF( MAPL_AM_I_ROOT() ) THEN
!                  print*,'+++SPEC,iGMI,iGEOS '//TRIM(speciesName), ic, i
!                END IF
! ########

            CASE ( GMI_TO_GEOS )
               IF ( i .LE. nqGG ) THEN
                 geos5BundleGG(i     )%data3d(:,:,1:km) = gmiBundle(ic)%pArray3D(:,:,km:1:-1)
               ELSE
                 geos5BundleXX(i-nqGG)%data3d(:,:,1:km) = gmiBundle(ic)%pArray3D(:,:,km:1:-1)
               END IF
            END SELECT
    
            found=.TRUE.
         END IF
      END DO
    
      RETURN

      END SUBROUTINE SwapSpeciesBundles
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: speciesReg_for_CCM
!
! !INTERFACE:
      FUNCTION speciesReg_for_CCM(nameSpecies, numSpecies, vnameGG, vnameXX ) RESULT(regCCM)
!
! !INPUT PARAMETERS:

    CHARACTER (LEN=*), INTENT(IN) :: nameSpecies(numSpecies) ! GMI mechanism list of species  - setkins names
    INTEGER,           INTENT(IN) :: numSpecies              ! number of setkins species

    CHARACTER (LEN=*), INTENT(IN) :: vnameGG(:)		     ! GEOS list of species - part A (transported)
    CHARACTER (LEN=*), INTENT(IN) :: vnameXX(:)		     ! GEOS list of species - part B (non-transported)
!
! !RETURNED VALUE:
    INTEGER :: regCCM(numSpecies)
!
! !DESCRIPTION:
! For each GMI species, this function determines its corresponding index
! in the CCM species list.
!
! !LOCAL VARIABLES:
    INTEGER :: ii, ic
    LOGICAL :: found
    CHARACTER(LEN=255) :: speciesName
    INTEGER, PARAMETER :: wrongIndex = -999
    INTEGER  nqGG, nqXX
!
!EOP
!------------------------------------------------------------------------------
!BOC
      regCCM(:) = wrongIndex

      nqGG = size(vnameGG)
      nqXX = size(vnameXX)

! Scan the GMI mechanism's species list, with numSpecies = NSP
! ------------------------------------------------------------
      GMIList: DO ic = 1, numSpecies

         speciesName = nameSpecies(ic)

! Substitute the GEOS-5 where it differs
! --------------------------------------
         IF(TRIM(speciesName) ==     "O3") speciesName = "OX"

         found = .FALSE.

! Look for a match in the GEOS set
! --------------------------------
         GEOS5ListGG: DO ii = 1, nqGG

            IF(TRIM(stringLowerCase(speciesName)) == TRIM(stringLowerCase(vnameGG(ii)))) THEN
               regCCM(ic) = ii
               found = .TRUE.
               EXIT GEOS5ListGG
            END IF

         END DO GEOS5ListGG

         GEOS5ListXX: DO ii = 1, nqXX

            IF(TRIM(stringLowerCase(speciesName)) == TRIM(stringLowerCase(vnameXX(ii)))) THEN
               regCCM(ic) = ii + nqGG
               found = .TRUE.
               EXIT GEOS5ListXX
            END IF

         END DO GEOS5ListXX

         IF (.NOT. found .AND. TRIM(nameSpecies(ic)) /= "H2O" ) THEN
            PRINT *,"speciesReg_for_CCM: "//TRIM(nameSpecies(ic))//" not found in GMIChem internal state"
         END IF

      END DO GMIList

      END FUNCTION speciesReg_for_CCM
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  CheckConfig
!
! !INTERFACE:

      SUBROUTINE CheckConfig(FeedBack_QV, rc)

      IMPLICIT NONE
!
! !OUTPUT PARAMETERS:
      LOGICAL, INTENT(OUT) :: FeedBack_QV
      INTEGER, INTENT(OUT) :: rc
!
! !DESCRIPTION:
!  Report the configuration settings
!    FeedBack_QV - default TRUE - whether to copy H2O field from GMI back into SPEC HUM field of GEOS
!
! !LOCAL VARIABLES:
      LOGICAL, SAVE              :: finished_reading = .FALSE.
      LOGICAL, SAVE              :: saved_feedback_qv

      INTEGER                    :: STATUS
      CHARACTER(LEN=255)         :: rcfilen = 'GMI_GridComp.rc'
      type (ESMF_Config)         :: gmiConfig
      character(len=ESMF_MAXSTR) :: Iam
!
!EOP
!---------------------------------------------------------------------------
!BOC
      rc=0
      Iam="GmiSwapSpeciesBundlesMod:CheckConfig"

      !!!
      !!!  NOTE: Consider making FeedBack_QV depend on the SimType
      !!!        i.e. if CTM then set to FALSE
      !!!

      IF ( .NOT. finished_reading ) THEN

        gmiConfig = ESMF_ConfigCreate(__RC__)

        call ESMF_ConfigLoadFile(gmiConfig, TRIM(rcfilen), __RC__)

        call ESMF_ConfigGetAttribute(gmiConfig, value=saved_feedback_qv, &
                                     label="FeedBack_QV:", default=.TRUE., __RC__)

        call ESMF_ConfigDestroy(gmiConfig, __RC__)

        finished_reading = .TRUE.

      END IF

      FeedBack_QV = saved_feedback_qv

      RETURN

      END SUBROUTINE CheckConfig
!EOC
!------------------------------------------------------------------------------

end module GmiSwapSpeciesBundlesMod
