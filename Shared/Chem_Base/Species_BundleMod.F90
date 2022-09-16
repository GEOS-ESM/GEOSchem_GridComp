#include "unused_dummy.H"
!-------------------------------------------------------------------------
!      NASA/GSFC, Global Modeling & Assimilation Office, Code 900.3      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Species_BundleMod.F90 --- Implements Species Bundle Class
! 
! !INTERFACE:
!
      MODULE  Species_BundleMod
            
! !USES:

      Use ESMF
      Use Runtime_RegistryMod
      Use Species_ArrayMod
      Use m_chars, only: uppercase
      Use MAPL, only: MAPL_UNDEF

      Implicit NONE

! !PUBLIC TYPES:
!
      PRIVATE 
      PUBLIC  Species_Bundle           ! species bundle type
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  Species_Grid             ! grid definition
      PUBLIC  Species_BundleCreate     ! initializes chemical bundle
      PUBLIC  Species_BundleDestroy    ! "cleans" chemical bundle
      PUBLIC  Species_BundleSetPtr     ! Fix internal pointers

!
! !DESCRIPTION: This module implements the species bundle class.
!
! !REVISION HISTORY: 
!
! 2022-09-13 manyin  First crack, based on Chem_Bundle
!EOP
!-------------------------------------------------------------------------

   real, parameter ::  missing_val = MAPL_UNDEF ! hardwire this for now
   integer, parameter :: nch = 256

!   Grid
!   ----
    type Species_Grid

!     Zonal grid
!     ----------
      integer       :: i1, i2, iml               ! local indices
      integer       :: ig                        ! ghosting
      integer       :: im                        ! global dimension
      integer       :: iLeft                     ! i1's index on global grid
      real, pointer :: lon(:,:) => null()        ! longitudes (deg)

!     Meridional grid
!     ---------------
      integer       :: j1, j2, jml               ! local indices
      integer       :: jg                        ! ghosting
      integer       :: jm                        ! global dimension
      real, pointer :: lat(:,:) => null()        ! latitudes (deg)

!     Vertical grid
!     -------------
      integer       :: km
      real, pointer :: lev(:) => null()
      character(len=nch) :: levUnits 
      real          :: ptop          ! Top pressure [Pa]

!     Horizontal gridbox area
!     -----------------------
      real, pointer :: cell_area(:,:) => null()

!     Cubed sphere or not
!     -------------------
      logical :: Cubed_Sphere = .FALSE.

    end type Species_Grid

!   Species vector
!   --------------
    type Species_Bundle

!     Registry
!     --------
      type(Runtime_Registry) :: reg

!     Grid
!     ----
      type(Species_Grid) :: grid

      real, pointer   :: cosz(:,:) => null()  ! cosine solar zenith angle
      real, pointer   :: sinz(:,:) => null()  !   sine solar zenith angle

      type(ESMF_Grid) :: grid_esmf

!     Whether this class allocated the memory for q, delp
!     ---------------------------------------------------
      logical :: did_allocation = .false.
      logical :: has_rh = .false.            ! for backward compatibility
      logical :: diurnal_bb = .false.        ! whether using diurnal biomass burning
      logical :: do_concentration = .false.  ! using concentration: MR * airdens instead of MR alone
      real    :: missing_value = MAPL_UNDEF

!     Tracer array
!     ------------
      real, pointer :: delp(:,:,:) => null()! Layer thickness [Pa] (not ghosted)
      real, pointer :: rh(:,:,:) => null()  ! Layer thickness [Pa] (not ghosted)
      real, pointer :: airdens(:,:,:) => null() ! Air density

      type(Species_array), pointer :: qa(:) => null()
                                            ! access 4D array in q as a 
                                            ! collection of 3D arrays; used
                                            ! for gradually removing the 4D
                                            ! arrays     
      integer :: nq

!     Two calendar elements (from ESMF)
!     ---------------------------------
      LOGICAL :: isLeapYear
      REAL :: dayOfYear

    end type Species_Bundle

   CONTAINS

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Species_BundleCreate --- Creates Species Bundle
! 
! !INTERFACE:
!
  subroutine  Species_BundleCreate  ( reg,                 &
                                   i1, i2, ig, im,      &
                                   j1, j2, jg, jm, km,  &
                                   w_s, rc,             &
                                   skipAlloc, lat, lon, &
                                   cell_area,           &
                                   lev, levUnits, ptop, &
                                   do_Conc )  ! Optional
!
! !USES:
!
  implicit NONE
!
! !INPUT PARAMETERS: 
!

  type(Runtime_registry)        :: reg             ! Species Registry

                                                ! Distributed grid info:
  integer,      intent(in)   :: i1, i2          !   local  zonal indices
  integer,      intent(in)   :: ig              !   zonal ghosting
  integer,      intent(in)   :: im              !   global zonal dimension
  integer,      intent(in)   :: j1, j2          !   local meridional indices
  integer,      intent(in)   :: jg              !   meridional ghosting
  integer,      intent(in)   :: jm              !   global zonal dimension
  integer,      intent(in)   :: km              !   vertical dimension

  logical, OPTIONAL, intent(in) :: skipAlloc    ! Do not allocate arrays
  real,    OPTIONAL, intent(in) :: lon(i1:i2,j1:j2) ! longitude in degrees
  real,    OPTIONAL, intent(in) :: lat(i1:i2,j1:j2) ! latitude in degrees
  real,    OPTIONAL, pointer    :: cell_area(:,:) ! grid box area
  real,    OPTIONAL, intent(in) :: lev(1:km)    ! levels
  character(len=*), OPTIONAL, intent(in) :: levUnits ! level units
  real,    OPTIONAL, intent(in) :: ptop         ! top pressure in Pa
  logical, OPTIONAL, intent(in) :: do_Conc ! Do_concentration: In case you need MR*airdens instead of MR alone
!
! !OUTPUT PARAMETERS:
!
  type(Species_Bundle), intent (out) :: w_s   ! Chemical bundle

  integer, intent(out)            :: rc    ! error return code:
                                           !  0 - all is well
                                           !  1 - already allocated
                                           !  2 - could not allocate memory
                                           !  3 - invalid dimensions

!
! !DESCRIPTION: Creates a Species Bundle, allocating the necessary memory or
!               optionally associating internal pointers with model declared 
!  flat arrays.
!
! !REVISION HISTORY: 
!
!  2022-09-13  manyin  First crack, based on Chem_Bundle
!EOP
!-------------------------------------------------------------------------

     character(len=*), parameter ::  myname = 'Species_BundleCreate'

     integer err, i, j, n, nq, ios, ios1, ios2, ios3, ios4
     logical :: do_allocation
     logical :: do_concentration
     real*8 :: delta

!    Sanity check
!    ------------
     rc = 0
     nq = reg%nq
     if ( im<1 .or. jm<1 .or. km<1 .or. nq<1) then
          rc = -3
          return
     endif
              
     w_s%reg = reg
     w_s%missing_value = MAPL_UNDEF

!    Whether or not we allocate memory for arrays
!    --------------------------------------------
     if ( present(skipAlloc) ) then
          do_allocation = .not. skipAlloc
     else
          do_allocation = .true.
     end if

!   Whether or not read airdens to be able to work with concentration instead  of MR
!   -----------------------------------------------
     if ( present(do_Conc) ) then
          do_concentration = do_Conc
     else
          do_concentration = .false.
     end if



!    Initialize dimensional attributes
!    ---------------------------------
     w_s%grid%i1 = i1; w_s%grid%i2 = i2; w_s%grid%ig = ig; w_s%grid%im = im
     w_s%grid%iml = i2 - i1 + 1 
     w_s%grid%j1 = j1; w_s%grid%j2 = j2; w_s%grid%jg = jg; w_s%grid%jm = jm
     w_s%grid%jml = j2 - j1 + 1 
     w_s%grid%km = km

!    Detect cubed sphere for sanity checks latter
!    --------------------------------------------
     if ( jm == im * 6 ) then
          w_s%grid%Cubed_Sphere = .TRUE.
     else
          w_s%grid%Cubed_Sphere = .FALSE.
     end if

!    Horizontal grid (hardwire A-grid for now)
!    -----------------------------------------
     if ( present(ptop) ) then
          w_s%grid%ptop =  ptop
     else
          w_s%grid%ptop =  1.0 ! Pa: reasonable default 
     endif

!    Save lat/lons
!    -------------
     allocate ( w_s%grid%lon(i1:i2,j1:j2), w_s%grid%lat(i1:i2,j1:j2), &
                stat = ios ) ! 
     if ( ios /= 0 ) then
        rc = 2
        return
     end if
     if ( present(lon) ) then
          w_s%grid%lon = lon
     else
          !ALT w_s%grid%lon = MAPL_UNDEF
          delta = 360.0d0/im
          do i = 1, im
            w_s%grid%lon(i,:) = -180.0d0 + (i-1)*delta
          end do
     end if

     if ( present(lat) ) then
          w_s%grid%lat = lat
     else
          !ALT w_s%grid%lat = MAPL_UNDEF
          if(jm==1) then
            delta = 0.0d0
          else
            delta = 180.0d0/(jm-1)
          endif
          do j = 1, jm
            w_s%grid%lat(:,j) = -90.0d0 + (j-1)*delta
          end do
     end if

     if ( present(cell_area) ) then ! will be left unallocated otherwise
          w_s%grid%cell_area => cell_area
     end if

     if ( present(lev) ) then
        allocate ( w_s%grid%lev(1:km), stat = ios )
        if ( ios /= 0 ) then
           rc = 3
           return
        end if
        w_s%grid%lev = lev 
     end if
     if ( present(levUnits) ) then
          w_s%grid%levUnits = levUnits
     else
          w_s%grid%levUnits = 'none'
     end if

     w_s%did_allocation = .false.
     w_s%has_rh = .false.       ! will be set to TRUE when set
     w_s%do_concentration = .false.

     w_s%nq = nq

     allocate(w_s%qa(nq),stat=ios2)

     if ( do_allocation ) then

        allocate(w_s%delp(i1:i2,j1:j2,km),stat=ios )
        do n = 1, nq
           allocate(w_s%qa(n)%data3d(i1-ig:i2+ig,j1-jg:j2+jg,km),stat=ios1 )
           ios2 = ios2 + ios1
        end do
        allocate(w_s%rh(i1:i2,j1:j2,km),stat=ios3 )
            
        if ( ios + ios2 + ios3 == 0 ) then 
             w_s%did_allocation = .true.
             w_s%delp = 0.0
             do n = 1, nq
                w_s%qa(n)%data3d = 0.0
             end do
             w_s%rh = 0.0             
        else
             w_s%did_allocation = .false.
             rc = 4
        end if

        if (do_concentration ) then
           w_s%do_concentration = .true.
           allocate(w_s%airdens(i1:i2,j1:j2,km),stat=ios4 )
           if ( ios4 == 0 ) w_s%airdens = 0.0
        endif

         
     end if

!    Set array of pointers: may be null() if no allocation took place
!    ----------------------------------------------------------------
     call Species_BundleSetPtr ( w_s, rc ) 

  
   end subroutine Species_BundleCreate

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Species_BundleDestroy --- Deallocates memory used by chemical state
! 
! !INTERFACE:
!
  subroutine  Species_BundleDestroy ( w_s, rc )
!
! !USES:
!
  implicit NONE
!
! !INPUT/OUTPUT PARAMETERS: 
!
  type(Species_Bundle), intent (inout) :: w_s   ! chemical bundle

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - 
! !DESCRIPTION: 
!
!  Deallocates memory used by chemical bundle.
!
! !REVISION HISTORY: 
!
!  20Jul1999 da Silva  Initial code.
!  26oct1999 da Silva  Added hs_stdv, ts, lwi, a, b
!  16Jun2016 Buchard   Added do_concentration option
!EOP
!-------------------------------------------------------------------------

   integer n, ier

   rc = 0

   if ( w_s%did_allocation ) then

      if ( associated(w_s%delp) ) deallocate(w_s%delp, stat=ier)
      if ( associated(w_s%rh)  )  deallocate(w_s%rh, stat=ier)

       
     
    
      do n = 1, w_s%reg%nq 
         if ( associated(w_s%qa(n)%data3d)  ) &
              deallocate(w_s%qa(n)%data3d, stat=ier)
      end do
      if(associated(w_s%grid%lon)) deallocate( w_s%grid%lon )
      if(associated(w_s%grid%lat)) deallocate( w_s%grid%lat )
      if(associated(w_s%grid%lev)) deallocate( w_s%grid%lev )
      if(associated(w_s%qa))       deallocate( w_s%qa )


      if ( w_s%do_concentration ) then
          if ( associated(w_s%airdens)  )  deallocate(w_s%airdens, stat=ier)
      endif

      
   else

      if ( associated(w_s%delp) ) nullify(w_s%delp)
      if ( associated(w_s%rh) )   nullify(w_s%rh)
            
      do n = 1, w_s%reg%nq 
         if ( associated(w_s%qa(n)%data3d)  ) nullify(w_s%qa(n)%data3d)
      end do
      deallocate( w_s%grid%lon, w_s%grid%lat, w_s%grid%lev,w_s%qa, stat=ier)  

      if ( w_s%do_concentration ) then
          if ( associated(w_s%airdens)  )  nullify(w_s%airdens)
      endif
 
   end if


   w_s%reg%nq = -1

  end subroutine Species_BundleDestroy

!-------------------------------------------------------------------------
!      NASA/GSFC Global Modeling & Assimilation Office, Code 900.3       !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Species_BundleSetPtr - Make sure internal pointers are set
! 
! !INTERFACE:
!
  subroutine  Species_BundleSetPtr ( w_s, rc )
!
! !USES:
!
  implicit NONE
!
! !INPUT/OUTPUT PARAMETERS: 
!
  type(Species_Bundle), intent (inout) :: w_s   ! chemical bundle

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - 
! !DESCRIPTION: 
!
!  Make sure the internal array of pointers points to the current
!  memory location of the 4D q array
!
! !REVISION HISTORY: 
!
!  22Jul2005 da Silva  Initial code.
!
!EOP
!-------------------------------------------------------------------------

   _UNUSED_DUMMY(w_s)
   rc = 0

   return

   end subroutine Species_BundleSetPtr

 end MODULE Species_BundleMod
