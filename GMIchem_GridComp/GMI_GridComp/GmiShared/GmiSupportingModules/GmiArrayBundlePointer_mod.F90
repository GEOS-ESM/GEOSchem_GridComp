
!-------------------------------------------------------------------------
! Module implements a derived type to support an "array of arrays" in F90.
! The derived type is not opaque - a compromise to minimize the impact on
! legacy code that uses the original data structure (a 4D array).
!-------------------------------------------------------------------------

module GmiArrayBundlePointer_mod

   implicit none
   private

   public :: t_GmiArrayBundle     ! derived type
   public :: ConstructArrayBundle ! constructor
   public :: CleanArrayPointer    ! destructor
   public :: setArrayPointer      ! accessor
   public :: getArrayPointer      ! accessor

   ! Derived type has public components for ease of use.
   ! Strong encapsulation would require intrusive changes
   ! in bothe GEOS and GMI.
   type t_GmiArrayBundle
      real*8, pointer  :: pArray3D(:,:,:) => null()
   end type t_GmiArrayBundle

   interface CleanArrayPointer
      module procedure CleanScalar
      module procedure CleanVector
   end interface

contains

!-------------------------------------------------------------------------

   ! Initialize a bundle from an already allocated 4D array.
   function ConstructArrayBundle(oldArray, &
                i1, i2, ju1, j2, k1, k2, numSpecies) result (bundle)
      integer :: i1, i2, ju1, j2, k1, k2
      integer :: numSpecies
      real*8, target :: oldArray(i1:i2,ju1:j2,k1:k2,numSpecies) ! I,J,K,species
      type (t_GmiArrayBundle), pointer :: bundle(:)

      integer :: i

      allocate(bundle(numSpecies))
      do i = 1, numSpecies
         call setArrayPointer(bundle(i), oldArray(i1,ju1,k1,i),i1 ,i2, ju1, j2, k1, k2)
      end do

   end function ConstructArrayBundle

!-------------------------------------------------------------------------

   ! Establish a link to a specified 3D target
   ! (pretend that encapsulation is enforced)
   subroutine setArrayPointer(this, newTarget, i1 ,i2, ju1, j2, k1, k2)
      integer :: i1, i2, ju1, j2, k1, k2
      type (t_GmiArrayBundle) :: this
      real*8, target :: newTarget(i1:i2,ju1:j2,k1:k2)

      ! The following could fail on some compilers if a temporary copy of
      ! oldArray is generated.  Unlikely in this context, though.
      this%pArray3D => newTarget

   end subroutine setArrayPointer

!-------------------------------------------------------------------------

   ! Retrieve pointer (pretend that encapsulation is enforced)
   function getArrayPointer(this) result(ptr)
      type (t_GmiArrayBundle) :: this
      real*8, pointer :: ptr(:,:,:)

      ptr => this%pArray3D

   end function getArrayPointer

!-------------------------------------------------------------------------

   ! Nullify the internal pointer (pretend encapsulation is enforced)
   subroutine CleanScalar(this)
      type (t_GmiArrayBundle) :: this
      nullify(this%pArray3D)
   end subroutine CleanScalar

!-------------------------------------------------------------------------

   ! Clean an array of derived types. First clean each element, and
   ! then deallocate the array.
   subroutine CleanVector(this, status)
      integer, intent(out) :: status
      type (t_GmiArrayBundle), pointer :: this(:)
      integer :: i

      do i = 1, size(this)
          deallocate(this(i)%pArray3D)
      end do
      deallocate(this, STAT=status)

   end subroutine CleanVector
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

end module GmiArrayBundlePointer_mod

