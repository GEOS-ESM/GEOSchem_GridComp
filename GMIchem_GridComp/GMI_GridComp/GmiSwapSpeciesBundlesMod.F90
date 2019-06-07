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
   USE MAPL_Mod

      use GmiStringManipulation_mod, only : stringLowerCase
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      Use Chem_ArrayMod
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

      SUBROUTINE SwapSpeciesBundles(direction, gmiBundle, geos5Bundle, &
                     specHum, mapSpecies, lchemvar, do_synoz, numSpecies, rc)

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      INTEGER, INTENT(IN) :: numSpecies
      INTEGER, INTENT(IN) :: direction       ! 1=To GMI  -1=FromGMI
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
      type(chem_array)        :: geos5Bundle(:)
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
      INTEGER :: i, ic
      INTEGER :: k, km
      LOGICAL :: found
      CHARACTER(LEN=255) :: speciesName
      CHARACTER(LEN=255) :: MyName
!
!EOP
!---------------------------------------------------------------------------
!BOC
      rc=0
      MyName="GMI_SwapSpeciesBundles"

      IF (ABS(direction) /= 1) THEN
         rc=1
         PRINT *,TRIM(MyName),": Invalid direction specified for swapping bundles."
         RETURN
      END IF

      km = size(specHum, 3)

      DO ic=1,numSpecies
         found=.FALSE.
         speciesName=TRIM(lchemvar(ic))
    
         IF (TRIM(speciesName) == "H2O") THEN
    
            SELECT CASE (direction)
            CASE ( 1)
               gmiBundle(ic)%pArray3D(:,:,km:1:-1) = specHum(:,:,1:km)*(MWTAIR/MWTH2O)
            CASE (-1)
               specHum(:,:,1:km) = gmiBundle(ic)%pArray3D(:,:,km:1:-1)*(MWTH2O/MWTAIR)
            END SELECT
    
            found=.TRUE.
    
         ELSE IF(TRIM(speciesName) == "SYNOZ" .OR. do_synoz .eqv. .TRUE.) THEN
    
            found=.FALSE.
    
         ELSE
            i = mapSpecies(ic)
    
            SELECT CASE (direction)
            CASE ( 1)
               gmiBundle(ic)%pArray3D(:,:,km:1:-1) = geos5Bundle(i)%data3d(:,:,1:km)
! ########
!                IF( MAPL_AM_I_ROOT() ) THEN
!                  print*,'+++SPEC,iGMI,iGEOS '//TRIM(speciesName), ic, i
!                END IF
! ########

            CASE (-1)
               geos5Bundle(i)%data3d(:,:,1:km) = gmiBundle(ic)%pArray3D(:,:,km:1:-1)
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
      FUNCTION speciesReg_for_CCM(nameSpecies, vname, numSpecies, iStart, iEnd) RESULT(regCCM)
!
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN) :: numSpecies ! number of GMI species
    INTEGER, INTENT(IN) :: iStart     ! Starting index for GMI transported species in w_c%qa
    INTEGER, INTENT(IN) :: iEnd       ! Last index for GMI non-transported species in w_c%qa
    CHARACTER (LEN=*), INTENT(IN) :: nameSpecies(numSpecies) ! GMI mechanism list of species
    CHARACTER (LEN=*), INTENT(IN) :: vname(*)		     ! GMIChem internal list of species
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
!
!EOP
!------------------------------------------------------------------------------
!BOC
      regCCM(:) = wrongIndex

! Scan the GMI mechanism's specie list, with numSpecies = NSP
! -----------------------------------------------------------
      GMIList: DO ic = 1, numSpecies

         speciesName = nameSpecies(ic)

! Substitute the GEOS-5 where it differs
! --------------------------------------
         IF(TRIM(speciesName) ==     "O3") speciesName = "OX"
         IF(TRIM(speciesName) ==  "CFCl3") speciesName = "CFC11"
         IF(TRIM(speciesName) == "CF2Cl2") speciesName = "CFC12"

         found = .FALSE.

! Look for a match in the Chem_Registry.rc
! ----------------------------------------
         GEOS5List: DO ii = iStart, iEnd

            IF(TRIM(stringLowerCase(speciesName)) == TRIM(stringLowerCase(vname(ii)))) THEN
               regCCM(ic) = ii
               found = .TRUE.
               EXIT GEOS5List
            END IF

         END DO GEOS5List

         IF (.NOT. found .AND. TRIM(nameSpecies(ic)) /= "H2O" ) THEN
            PRINT *,"speciesReg_for_CCM: "//TRIM(nameSpecies(ic))//" not found in GMIChem internal state"
         END IF

      END DO GMIList

      END FUNCTION speciesReg_for_CCM
!EOC
!------------------------------------------------------------------------------
end module GmiSwapSpeciesBundlesMod
