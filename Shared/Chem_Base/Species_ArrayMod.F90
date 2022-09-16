Module Species_ArrayMod

   private
   public Species_Array

!  The Species_Array, a light weight ESMF-like array
!  -------------------------------------------------
   type Species_Array
!       integer :: rank
!       logical :: did_alloc = .false. ! useful to keep track of allocations
!       real, pointer :: data2d(:,:)   => null()
        real, pointer :: data3d(:,:,:) => null()

!       A per-tracer scavenging efficiency in convective updrafts [km-1]
!       real :: fscav = 0.0

!       A per-tracer large-scale wet removal efficiency [fraction]
!       real :: fwet = 0.0

   end type Species_Array

 end Module Species_ArrayMod
