 module GmiSpeciesRegistry_mod

 use GmiPrintError_mod        , only : GmiPrintError
 use GmiStringManipulation_mod, only : stringLowerCase

 implicit none

 private
 public  :: getSpeciesIndex
 public  :: set_numSpecies, set_labelsSpecies
 public  :: set_molWeightSpecies, get_molWeightSpecies
 public  :: UNKNOWN_SPECIES

# include "GmiParameters.h"
!# include "setkin_par.h"
!# include "setkin_lchem.h"

 integer, parameter :: UNKNOWN_SPECIES         = -2

 integer :: numSpecies
 real*8                                  , pointer :: molWeightSpecies(:) => null()
 character (len=MAX_LENGTH_SPECIES_NAME) , pointer :: labelsSpecies   (:) => null()

 contains

!----------------------------------------------------

 subroutine set_numSpecies(num_species)

 integer, intent(in) :: num_species

 numSpecies = num_species

 return

 end subroutine set_numSpecies

!----------------------------------------------------

 subroutine set_labelsSpecies(const_labels)

 character (len=MAX_LENGTH_SPECIES_NAME), intent(in) :: const_labels(numSpecies)

 allocate(labelsSpecies(numSpecies))

 labelsSpecies(1:numSpecies) = const_labels(1:numSpecies)

 return

 end subroutine set_labelsSpecies

!----------------------------------------------------

 subroutine set_molWeightSpecies(mw)

 real*8, intent(in) :: mw(numSpecies)

 allocate(molWeightSpecies(numSpecies))
 
 molWeightSpecies(1:numSpecies) = mw(1:numSpecies)

 return

 end subroutine set_molWeightSpecies

!----------------------------------------------------

 subroutine get_molWeightSpecies(mw)

 real*8, intent(out) :: mw(numSpecies)

 allocate(molWeightSpecies(numSpecies))
 
 mw(1:numSpecies) = molWeightSpecies(1:numSpecies)

 return

 end subroutine get_molWeightSpecies

!----------------------------------------------------

 function getSpeciesIndex(name) result(index)
 character (len=*), intent(in) :: name
 integer                      :: index
 character (len=128)  :: err_msg

 integer :: ii

 index = UNKNOWN_SPECIES

 if (trim(stringLowerCase(name)) == 'xxx') then
    index = -1  ! this is mainly used for emiss_map
 else
    do ii = 1, numSpecies
       if (trim(stringLowerCase(labelsSpecies(ii))) == trim(stringLowerCase(name))) then
          index = ii
          exit
       end if
    end do
 end if

 if (index == UNKNOWN_SPECIES) then
    err_msg = 'The species does not exist: '// name
    call GmiPrintError(err_msg, .true., 1, index, 0, 0, 0.0d0, 0.0d0)
 end if

 end function getSpeciesIndex

!----------------------------------------------------
!
! function getSpeciesIndex(name) result(index)
! character (len=*), intent(in) :: name
! integer                      :: index
!
! integer :: ii
!
! index = UNKNOWN_SPECIES
!
! do ii = 1, NSP
!    if (trim(lchemvar(ii)) == trim(name)) then
!       index = ii
!       exit
!    end if
! end do
!
! end function getSpeciesIndex
!
!----------------------------------------------------

 end module GmiSpeciesRegistry_mod
