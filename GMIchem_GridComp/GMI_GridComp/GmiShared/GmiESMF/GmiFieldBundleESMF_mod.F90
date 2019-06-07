#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiFieldBundleESMF_mod
!
! !INTERFACE:
!
      module GmiFieldBundleESMF_mod
!
! !USES:
      use ESMF
      use MAPL_Mod

      implicit none

! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: addTracerToBundle
      public  :: updateTracerToBundle
      public  :: obtainTracerFromBundle

      interface addTracerToBundle
         module procedure addTracerToBundle3D_r4
         module procedure addTracerToBundle3D_r8
         module procedure addTracerToBundle2D_r4
         module procedure addTracerToBundle2D_r8
      end interface

      interface updateTracerToBundle
         module procedure updateTracerToBundle_ByIndex3D_r4
         module procedure updateTracerToBundle_ByIndex3D_r8
         module procedure updateTracerToBundle_ByName3D_r4
         module procedure updateTracerToBundle_ByName3D_r8
         module procedure updateTracerToBundle_ByIndex2D_r4
         module procedure updateTracerToBundle_ByIndex2D_r8
         module procedure updateTracerToBundle_ByName2D_r4
         module procedure updateTracerToBundle_ByName2D_r8
      end interface

      interface obtainTracerFromBundle
         module procedure obtainTracerFromBundle_ByIndex3D_r4
         module procedure obtainTracerFromBundle_ByIndex3D_r8
         module procedure obtainTracerFromBundle_ByName3D_r4
         module procedure obtainTracerFromBundle_ByName3D_r8
         module procedure obtainTracerFromBundle_ByIndex2D_r4
         module procedure obtainTracerFromBundle_ByIndex2D_r8
         module procedure obtainTracerFromBundle_ByName2D_r4
         module procedure obtainTracerFromBundle_ByName2D_r8
      end interface

#     include "GmiParameters.h"

! !DESCRIPTION:
! Basic routines for manipulating ESMF bundles.
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine addTracerToBundle3D_r8(bundle, PTR, grid, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
      type (ESMF_Grid),           intent(in) :: grid
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Adds tracers to a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type (ESMF_Field)            :: field
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:addTracerToBundle3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      field = ESMF_FieldCreate(grid=grid, fArrayptr=PTR, &
              name = fieldName, RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_FieldBundleAdd ( bundle, field, rc=STATUS )
      VERIFY_(STATUS)

      return

      end subroutine addTracerToBundle3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine addTracerToBundle3D_r4(bundle, PTR, grid, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
      type (ESMF_Grid),           intent(in) :: grid
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Adds tracers to a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type (ESMF_Field)            :: field
      character(len=ESMF_MAXSTR), parameter   :: IAm='GMI:addTracerToBundle3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      field = ESMF_FieldCreate(grid=grid, fArrayptr=PTR, &
              name = fieldName, RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_FieldBundleAdd ( bundle, field, rc=STATUS )
      VERIFY_(STATUS)

      return

      end subroutine addTracerToBundle3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName3D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r8), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByName3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=STATUS)
      VERIFY_(STATUS)

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByName3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName3D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r4), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByName3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=STATUS)
      VERIFY_(STATUS)

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByName3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex3D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r8), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByIndex3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=STATUS)
      VERIFY_(STATUS)

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex3D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r4), pointer :: ptr3D(:,:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByIndex3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr3D, RC=STATUS)
      VERIFY_(STATUS)

      ptr3D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByName3D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByName3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
     call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
     VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByName3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByName3D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByName3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
     call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
     VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByName3D_r4
!!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByIndex3D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByIndex3D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByIndex3D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByIndex3D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByIndex3D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByIndex3D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine addTracerToBundle2D_r8(bundle, PTR, grid, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
      type (ESMF_Grid),           intent(in) :: grid
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Adds tracers to a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type (ESMF_Field)            :: field
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:addTracerToBundle2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      field = ESMF_FieldCreate(grid=grid, fArrayptr=PTR, &
              name = fieldName, RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_FieldBundleAdd ( bundle, field, rc=STATUS )
      VERIFY_(STATUS)

      return

      end subroutine addTracerToBundle2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine addTracerToBundle2D_r4(bundle, PTR, grid, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
      type (ESMF_Grid),           intent(in) :: grid
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Adds tracers to a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type (ESMF_Field)            :: field
      character(len=ESMF_MAXSTR), parameter   :: IAm='GMI:addTracerToBundle2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      field = ESMF_FieldCreate(grid=grid, fArrayptr=PTR, &
              name = fieldName, RC=STATUS )
      VERIFY_(STATUS)

      call MAPL_FieldBundleAdd ( bundle, field, rc=STATUS )
      VERIFY_(STATUS)

      return

      end subroutine addTracerToBundle2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName2D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r8), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByName2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=STATUS)
      VERIFY_(STATUS)

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByName2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByName2D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r4), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByName2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=STATUS)
      VERIFY_(STATUS)

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByName2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex2D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r8), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByIndex2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=STATUS)
      VERIFY_(STATUS)

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine updateTracerToBundle_ByIndex2D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
      integer,                    intent(in) :: index
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Updates a tracer in a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      real(r4), pointer :: ptr2D(:,:)
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:updateTracerToBundle_ByIndex2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=ptr2D, RC=STATUS)
      VERIFY_(STATUS)

      ptr2D = PTR

      return

      end subroutine updateTracerToBundle_ByIndex2D_r4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByName2D_r8(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByName2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
     call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
     VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByName2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByName2D_r4(bundle, PTR, fieldName)
!
! !INPUT PARAMETERS:
      character(len=ESMF_MAXSTR), intent(in) :: fieldName
!
! !OUTPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByName2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
     call ESMF_FieldBundleGet (bundle, fieldName, field=field, RC=STATUS)
     VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByName2D_r4
!!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByIndex2D_r8(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(r8), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByIndex2D_r8'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByIndex2D_r8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
      subroutine obtainTracerFromBundle_ByIndex2D_r4(bundle, PTR, index)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: index
!
! !OUTPUT PARAMETERS:
      real(r4), pointer :: PTR(:,:)
!
! !INPUT/OUTPUT PARAMETERS:
      type (ESMF_FieldBundle), intent(inOut) :: bundle
!
! !DESCRIPTION:
! Obtains a tracer from a bundle.
!
! !LOCAL VARIABLES:
      integer :: rc, status
      type(ESMF_FIELD)           :: field
      type(ESMF_ARRAY)           :: array
      character(len=ESMF_MAXSTR), parameter :: IAm='GMI:obtainTracerFromBundle_ByIndex2D_r4'
!
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_FieldBundleGet (bundle, fieldIndex=index, field=field, RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_FieldGet  (field, array=array,                 RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ArrayGet  (array, farrayptr=PTR, RC=STATUS)
      VERIFY_(STATUS)

      return

      end subroutine obtainTracerFromBundle_ByIndex2D_r4
!EOC
!------------------------------------------------------------------------------

      end module GmiFieldBundleESMF_mod
