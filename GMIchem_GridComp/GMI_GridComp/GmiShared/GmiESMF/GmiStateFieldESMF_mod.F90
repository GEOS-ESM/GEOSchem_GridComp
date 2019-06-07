#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiStateFieldESMF_mod
!
! !INTERFACE:
!
      module GmiStateFieldESMF_mod
!
! !USES:
      use ESMF
      use MAPL_Mod

      implicit none

! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: setDataToStateField
      public  :: getDataFromStateField
      public  :: initDataInStateField

      interface initDataInStateField
         module procedure initDataInStateField3D_r4
         module procedure initDataInStateField3D_r8
         module procedure initDataInStateField2D_r4
         module procedure initDataInStateField2D_r8
      end interface
!
      interface setDataToStateField
         module procedure setDataToStateField3D_r4
         module procedure setDataToStateField3D_r8
         module procedure setDataToStateField2D_r4
         module procedure setDataToStateField2D_r8
      end interface
!
      interface getDataFromStateField
         module procedure getDataFromStateField3D_r4
         module procedure getDataFromStateField3D_r8
         module procedure getDataFromStateField2D_r4
         module procedure getDataFromStateField2D_r8
      end interface
!
#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Basic routines for manipulating fields in a ESMF state.
! Fortran arrays can be 2D or 3D and be in single or double precision.
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initDataInStateField3D_r4
!
! !INTERFACE:
!
      subroutine initDataInStateField3D_r4(state, grid, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
      type (ESMF_grid),  intent(in) :: grid
      character(len=*),  intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Initializes a ESMF state field from a predefined pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:initDataInStateField3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

!ALT      call ESMF_FieldSetCommit  (field, grid, farrayPtr=PTR, RC=STATUS)
          call ESMF_FieldEmptyComplete  (field, farrayPtr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine initDataInStateField3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initDataInStateField3D_r8
!
! !INTERFACE:
!
      subroutine initDataInStateField3D_r8(state, grid, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
      type (ESMF_grid),  intent(in) :: grid
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Initializes a ESMF state field from a predefined pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:initDataInStateField3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

!ALT      call ESMF_FieldSetCommit  (field, grid, farrayPtr=PTR, RC=STATUS)
          call ESMF_FieldEmptyComplete  (field, farrayPtr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine initDataInStateField3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initDataInStateField2D_r4
!
! !INTERFACE:
!
      subroutine initDataInStateField2D_r4(state, grid, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
      type (ESMF_grid),  intent(in) :: grid
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Initializes a ESMF state field from a predefined pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:initDataInStateField2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

!ALT      call ESMF_FieldSetCommit  (field, grid, farrayPtr=PTR, RC=STATUS)
          call ESMF_FieldEmptyComplete  (field, farrayPtr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine initDataInStateField2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initDataInStateField2D_r8
!
! !INTERFACE:
!
      subroutine initDataInStateField2D_r8(state, grid, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
      type (ESMF_grid),  intent(in) :: grid
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Initializes a ESMF state field from a predefined pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:initDataInStateField2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

!ALT      call ESMF_FieldSetCommit  (field, grid, farrayPtr=PTR, RC=STATUS)
          call ESMF_FieldEmptyComplete  (field, farrayPtr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine initDataInStateField2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setDataToStateField3D_r4
!
! !INTERFACE:
!
      subroutine setDataToStateField3D_r4(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Updates a ESMF state field from a provided pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r4), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:setDataToStateField3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=STATUS)
      VERIFY_(STATUS)

      ptr3D(:,:,:) = PTR(:,:,:)
      
      return

      end subroutine setDataToStateField3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setDataToStateField3D_r8
!
! !INTERFACE:
!
      subroutine setDataToStateField3D_r8(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Updates a ESMF state field from a provided pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r8), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:setDataToStateField3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=STATUS)
      VERIFY_(STATUS)

      ptr3D(:,:,:) = PTR(:,:,:)

      return

      end subroutine setDataToStateField3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setDataToStateField2D_r4
!
! !INTERFACE:
!
      subroutine setDataToStateField2D_r4(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Updates a ESMF state field from a provided pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r4), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:setDataToStateField2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=STATUS)
      VERIFY_(STATUS)

      ptr2D(:,:) = PTR(:,:)

      return

      end subroutine setDataToStateField2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setDataToStateField2D_r8
!
! !INTERFACE:
!
      subroutine setDataToStateField2D_r8(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
      character(len=*), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Updates a ESMF state field from a provided pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r8), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:setDataToStateField2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=STATUS)
      VERIFY_(STATUS)

      ptr2D(:,:) = PTR(:,:)

      return

      end subroutine setDataToStateField2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getDataFromStateField3D_r4
!
! !INTERFACE:
!
      subroutine getDataFromStateField3D_r4(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Gets data from a ESMF state field by passing it to an empty pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:getDataFromStateField3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)
      
      return

      end subroutine getDataFromStateField3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getDataFromStateField3D_r8
!
! !INTERFACE:
!
      subroutine getDataFromStateField3D_r8(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Gets data from a ESMF state field by passing it to an empty pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:getDataFromStateField3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine getDataFromStateField3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getDataFromStateField2D_r4
!
! !INTERFACE:
!
      subroutine getDataFromStateField2D_r4(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Gets data from a ESMF state field by passing it to an empty pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:getDataFromStateField2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine getDataFromStateField2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getDataFromStateField2D_r8
!
! !INTERFACE:
!
      subroutine getDataFromStateField2D_r8(state, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_State), intent(inOut) :: state
!
! !DESCRIPTION:
! Gets data from a ESMF state field by passing it to an empty pointer.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:getDataFromStateField2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_StateGet (state, TRIM(fieldName), field=field, rc=status)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine getDataFromStateField2D_r8
!EOC
!------------------------------------------------------------------------------

      end module GmiStateFieldESMF_mod
