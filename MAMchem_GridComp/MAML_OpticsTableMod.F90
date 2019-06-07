#include "MAPL_Exceptions.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1    !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  MAML_OpticsTableMod --- Reader for aerosol optical properties 
!           lookup tables.
!
! !INTERFACE:
!

   module  MAML_OpticsTableMod

! !USES:

   use ESMF
   use MAPL_Mod
   use m_die, only: die, warn

   use netcdf

   implicit none

! !PUBLIC TYPES:

   private

   public  MAML_OpticsTable           ! Holds Optics Lookup Tables
                           
!
! !PUBLIC MEMBER FUNCTIONS:
!
   public  MAML_OpticsTableCreate     ! Constructor - create a new optics table
   public  MAML_OpticsTableDestroy    ! Destructor - release resources associated with an optics table
   public  MAML_OpticsTableRead       ! Read optics table from a file

!
! !DESCRIPTION:
!
!  This module reads optics/mie aerosol tables.
!
! !REVISION HISTORY:
!
!  01Apr2013 Darmenov - Initial code based on Chem_MieTableMod.
!
!EOP
!-------------------------------------------------------------------------

! Optics/Mie LUT table
! Will be reduced from input files to the desired channels
! --------
  type MAML_OpticsTable
     character(len=1024):: file                                ! lookup table file 

     ! mode properties
     logical            :: monochromatic                       ! flag that indicates monochromatic or band-averaged quantities

     real               :: Dgs_min                             ! diameter of surface distribution - lower bound
     real               :: Dgs_max                             ! diameter of surface distribution - upper bound
     real               :: sigma                               ! geometric standard deviation of lognormal size distribution
     real               :: rh_deliq                            ! RH deliquescence point
     real               :: rh_cryst                            ! RH crystallization point
     integer            :: n_aerosol_components                ! number of aerosol components
     character(len=80), pointer &
                        :: component(:) => null()              ! aerosol components
     real, pointer      :: component_refractive_index_re(:,:)  ! real      part of aerosol component refractive index
     real, pointer      :: component_refractive_index_im(:,:)  ! imaginary part of aerosol component refractive index

     ! mie 
     integer            :: n_cheb                              ! number of terms used in truncated Chebyshev series expansion
     integer            :: n_bands                             ! number of monochromatic (narrow-band) or heterochromatic (wideband) bands
     integer            :: n_refractive_index_re               ! number of refractive indexes - real part
     integer            :: n_refractive_index_im               ! number of refractive indexes - imaginary part

     real, pointer      :: wavelength(:,:)          => null()  ! band range, lower- and upper-bound wavelengths
     real, pointer      :: refractive_index_re(:,:) => null()  ! real part of refractive index
     real, pointer      :: refractive_index_im(:,:) => null()  ! imaginary part of refractive index
     real, pointer      :: c_ext(:,:,:,:)           => null()  ! chebyshev polynomial coefficients of specific extinction
     real, pointer      :: c_sca(:,:,:,:)           => null()  ! chebyshev polynomial coefficients of specific scattering
     real, pointer      :: c_asy(:,:,:,:)           => null()  ! chebyshev polynomial coefficients of asymmetry parameter

!    TODO...
!    integer            :: n_mom                               ! number of moments of phase function
!    integer            :: n_pol                               ! number of elements of scattering phase matrix
!    real, pointer      :: pback(:,:,:,:)  => null()           ! backscatter phase function
!    real, pointer      :: pmom(:,:,:,:,:) => null()           ! moments of phase function
  end type MAML_OpticsTable


contains

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAML_OpticsTableCreate --- Construct Chemistry Registry
!
! !INTERFACE:
!

  function MAML_OpticsTableCreate(file, rc) result (this)

   implicit none

   type(MAML_OpticsTable) :: this

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
   character(len=*), intent(in) :: file

! !OUTPUT PARAMETERS:
   integer, intent(out) :: rc  ! Error return code:
                               !  0 - all is well
                               !  1 - 

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
!  01Apr2013 Darmenov    API.
!
!EOP
!-------------------------------------------------------------------------

                           __Iam__("Chem_MieTableCreate")

   type(MAML_OpticsTable) :: lut

   lut%file = trim(file)

   lut%monochromatic = .false.

   lut%Dgs_min       = -1.0
   lut%Dgs_max       = -1.0
   lut%sigma         = -1.0
   lut%rh_deliq      = -1.0
   lut%rh_cryst      = -1.0

   lut%n_aerosol_components = -1

   lut%n_cheb        = -1
   lut%n_bands       = -1

   lut%n_refractive_index_re  = -1
   lut%n_refractive_index_im  = -1

!  lut%n_mom         = -1
!  lut%n_pol         = -1

   this = lut

   RETURN_(ESMF_SUCCESS)

  end function MAML_OpticsTableCreate


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MAML_OpticsTableDestroy --- Release resources associated with the table
!
! !INTERFACE:
!
  subroutine MAML_OpticsTableDestroy(this, rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(MAML_OpticsTable), intent(inout) :: this

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   integer, intent(out) :: rc    ! Error return code:
                                 !  0 - all is well
                                 !  1 - 

! !DESCRIPTION: Destructor for MAML_OpticsTable object.
!
! !REVISION HISTORY:
!
!  01Apr2013 Darmenov    API.
!
!EOP
!-------------------------------------------------------------------------

                           __Iam__("MAML_OpticsTableDestroy")


!  Set these to invalid values
!  ---------------------------
   this%file = ''

   this%monochromatic = .false.
   this%Dgs_min  = -1.0
   this%Dgs_max  = -1.0
   this%sigma    = -1.0
   this%rh_deliq = -1.0
   this%rh_cryst = -1.0
   this%n_aerosol_components = -1
   this%n_cheb   = -1
   this%n_bands  = -1
   this%n_refractive_index_re  = -1
   this%n_refractive_index_im  = -1

!  this%n_mom = -1
!  this%n_pol = -1

!  Deallocate whatever has been allocated
!  --------------------------------------
   if (associated(this%component)) &
       deallocate(this%component, __STAT__)

   if (associated(this%component_refractive_index_re)) &
       deallocate(this%component_refractive_index_re, __STAT__)

   if (associated(this%component_refractive_index_im)) &
       deallocate(this%component_refractive_index_im, __STAT__) 

   if (associated(this%refractive_index_re)) &
       deallocate(this%refractive_index_re, __STAT__)

   if (associated(this%refractive_index_im)) &
       deallocate(this%refractive_index_im, __STAT__)

   if (associated(this%wavelength)) &
       deallocate(this%wavelength, __STAT__)

   if (associated(this%c_ext)) &
       deallocate(this%c_ext, __STAT__)

   if (associated(this%c_sca)) &
       deallocate(this%c_sca, __STAT__)

   if (associated(this%c_asy)) &
       deallocate(this%c_asy, __STAT__)


!  if (associated(this%bbck))       deallocate(this%bbck,   __STAT__)
!  if (associated(this%pmom))       deallocate(this%pmom,   __STAT__)

!  All done
!  --------
   
   RETURN_(ESMF_SUCCESS)
   
  end subroutine MAML_OpticsTableDestroy 


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAML_OpticsTableRead --- Read and fill in the Optics/Mie table
!
! !INTERFACE:
!
  subroutine MAML_OpticsTableRead(this, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   type(MAML_OpticsTable), intent(inout) :: this

! !INPUT PARAMETERS:   

! !OUTPUT PARAMETERS:
   integer, intent(out)                  :: rc         ! return code


! !DESCRIPTION: Fills in the Optics/Mie table
!
! !REVISION HISTORY:
!
! 01Mar2013 Darmenov    API.
!
!EOP
!-------------------------------------------------------------------------

   integer :: id_nc, id_dim, id_var                 ! IDs for NetCDF objects

   integer :: id_dim_band, len_dim_band             ! length of the wavelength dimension
   integer :: id_dim_range, len_dim_range           ! length of the range dimension
   integer :: id_dim_component, len_dim_component   ! length of the aerosol components dimension
   integer :: id_dim_chars, len_dim_chars           ! length of the characters dimension
   integer :: id_dim_n_re, len_dim_n_re             ! length of the real part of refractive index dimension
   integer :: id_dim_n_im, len_dim_n_im             ! length of the imaginary part of refractive index dimension
   integer :: id_dim_k, len_dim_k                   ! length of the truncated Chebyshev series dimension

   integer, dimension(nf90_max_var_dims) :: var_dim_ids
   integer, dimension(nf90_max_var_dims) :: var_dim_len
   integer :: dim_len
   integer :: n_dims

   character(len=80) :: buff_str

                           __Iam__("MAML_OpticsTableRead")

!  Open the lookup table file
!  --------------------------
   rc = nf90_open(trim(this%file), NF90_NOWRITE, id_nc); VERIFY_(rc)

!  Read ID and size of dimensions
!  -------------------------------
   call get_dim_info_(id_nc, 'band',      id_dim_band,      len_dim_band,      __RC__)
   call get_dim_info_(id_nc, 'range',     id_dim_range,     len_dim_range,     __RC__)
   call get_dim_info_(id_nc, 'nchars',    id_dim_chars,     len_dim_chars,     __RC__)
   call get_dim_info_(id_nc, 'component', id_dim_component, len_dim_component, __RC__)
   call get_dim_info_(id_nc, 'n_re',      id_dim_n_re,      len_dim_n_re,      __RC__)
   call get_dim_info_(id_nc, 'n_im',      id_dim_n_im,      len_dim_n_im,      __RC__)
   call get_dim_info_(id_nc, 'k',         id_dim_k,         len_dim_k,         __RC__)

   ASSERT_(len_dim_range == 2)

   this%n_bands               = len_dim_band
   this%n_aerosol_components  = len_dim_component
   this%n_cheb                = len_dim_k
   this%n_refractive_index_re = len_dim_n_re
   this%n_refractive_index_im = len_dim_n_im 

!  Read lookup table attributes
!  ----------------------------
   rc = nf90_get_att(id_nc, NF90_GLOBAL, 'aerosol_method', buff_str); VERIFY_(rc)
   ASSERT_(buff_str == 'modal')
 
   rc = nf90_get_att(id_nc, NF90_GLOBAL, 'optics', buff_str); VERIFY_(rc)
   ASSERT_(buff_str == 'monochromatic' .or. buff_str == 'band averaged')
   if (buff_str == 'monochromatic') then
       this%monochromatic = .true.
   else    
       this%monochromatic = .false.
   end if    

   rc = nf90_get_att(id_nc, NF90_GLOBAL, 'Dgs_min',    this%Dgs_min)
   VERIFY_(rc)

   rc = nf90_get_att(id_nc, NF90_GLOBAL, 'Dgs_max',    this%Dgs_max)
   VERIFY_(rc)

   rc = nf90_get_att(id_nc, NF90_GLOBAL, 'mode_width', this%sigma)
   VERIFY_(rc)

   rc = nf90_get_att(id_nc, NF90_GLOBAL, 'mode_deliq', this%rh_deliq)
   VERIFY_(rc)

   rc = nf90_get_att(id_nc, NF90_GLOBAL, 'mode_cryst', this%rh_cryst)
   VERIFY_(rc)


!  Read data - aerosol components
!  -----------------------------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'component', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 2)
   ASSERT_(var_dim_ids(1) == id_dim_chars .and. var_dim_ids(2) == id_dim_component)

   allocate(this%component(len_dim_component), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%component); VERIFY_(rc)


!  Read data - real part of aerosol component refractive index
!  -----------------------------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'component_n_re', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 2)
   ASSERT_(var_dim_ids(1) == id_dim_band .and. var_dim_ids(2) == id_dim_component)

   allocate(this%component_refractive_index_re(len_dim_band, len_dim_component), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%component_refractive_index_re); VERIFY_(rc)


!  Read data - imaginary part of aerosol component refractive index
!  -----------------------------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'component_n_im', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 2)
   ASSERT_(var_dim_ids(1) == id_dim_band .and. var_dim_ids(2) == id_dim_component)

   allocate(this%component_refractive_index_im(len_dim_band, len_dim_component), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%component_refractive_index_im); VERIFY_(rc)
 
!  Read data - band range
!  ----------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'wavelength', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 2)
   ASSERT_(var_dim_ids(1) == id_dim_band .and. var_dim_ids(2) == id_dim_range)

   allocate(this%wavelength(len_dim_band, len_dim_range), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%wavelength); VERIFY_(rc)

   if (this%monochromatic) then
       ASSERT_(maxval(abs(this%wavelength(:,1) - this%wavelength(:,2))) < 1e2*tiny(0.0))
   end if

!  Read data - real part of refractive index
!  -----------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'n_re', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 2)
   ASSERT_(var_dim_ids(1) == id_dim_band .and. var_dim_ids(2) == id_dim_n_re)

   allocate(this%refractive_index_re(len_dim_band, len_dim_n_re), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%refractive_index_re); VERIFY_(rc)

!  Read data - imaginary part of refractive index
!  ----------------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'n_im', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 2)
   ASSERT_(var_dim_ids(1) == id_dim_band .and. var_dim_ids(2) == id_dim_n_im)

   allocate(this%refractive_index_im(len_dim_band, len_dim_n_im), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%refractive_index_im); VERIFY_(rc)

!  Read data - chebyshev polynomial coefficients of specific extinction
!  --------------------------------------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'c_ext', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 4)
   ASSERT_(all(var_dim_ids(:4) .eq. (/id_dim_band, id_dim_n_im, id_dim_n_re, id_dim_k/)))

   allocate(this%c_ext(len_dim_band, len_dim_n_im, len_dim_n_re, len_dim_k), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%c_ext); VERIFY_(rc)

!  Read data - chebyshev polynomial coefficients of specific scattering
!  --------------------------------------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'c_sca', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 4)
   ASSERT_(all(var_dim_ids(:4) .eq. (/id_dim_band, id_dim_n_im, id_dim_n_re, id_dim_k/)))

   allocate(this%c_sca(len_dim_band, len_dim_n_im, len_dim_n_re, len_dim_k), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%c_sca); VERIFY_(rc)

!  Read data - chebyshev polynomial coefficients of asymmetry parameter
!  --------------------------------------------------------------------
   var_dim_ids = 0
   var_dim_len = 0

   rc = nf90_inq_varid(id_nc, 'c_asy', id_var); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, ndims=n_dims); VERIFY_(rc)
   rc = nf90_inquire_variable(id_nc, id_var, dimids=var_dim_ids(:n_dims)); VERIFY_(rc)

   ASSERT_(n_dims == 4)
   ASSERT_(all(var_dim_ids(:4) .eq. (/id_dim_band, id_dim_n_im, id_dim_n_re, id_dim_k/)))

   allocate(this%c_asy(len_dim_band, len_dim_n_im, len_dim_n_re, len_dim_k), __STAT__)
   rc = nf90_get_var(id_nc, id_var, this%c_asy); VERIFY_(rc)


!  Get the backscatter phase function values
!  TODO...


!  Close the table file
!  -------------------------------------
   rc = nf90_close(id_nc); VERIFY_(rc)


   RETURN_(ESMF_SUCCESS)

   contains

       subroutine get_dim_info_(id_nc, dim_name, dim_id, dim_len, rc)
           implicit none

           integer, intent(in)          :: id_nc
           character(len=*), intent(in) :: dim_name
           integer, intent(out)         :: dim_id
           integer, intent(out)         :: dim_len
           integer, intent(out)         :: rc

           rc = nf90_inq_dimid(id_nc, trim(dim_name), dim_id); VERIFY_(rc)
           rc = nf90_inquire_dimension(id_nc, dim_id, len=dim_len); VERIFY_(rc)

           RETURN_(ESMF_SUCCESS)
       end subroutine get_dim_info_

   end subroutine MAML_OpticsTableRead


   subroutine polint(x, y, n, xWant, yWant, yErr, myname)

       implicit none

       integer          :: n
       real(kind=8)     :: x(n), y(n)
       real             :: xWant, yWant, yErr
       character(len=*) :: myname

       ! given array x(n) of independent variables and array y(n) of dependent
       ! variables, compute the linear interpolated result yWant at xWant and return
       ! with a dummy error estimate yErr.  Hacked up from Numerical Recipes Chapter 3

       integer :: i, j
       real    :: dx, slope
       character(len=255) :: msg

       ! on out of bounds, set i to lower or upper limit
       i = 0
       if(xWant < x(1)) then 
           write(msg,*) "in polint, wanted: ", xWant, ", got lower bound: ", x(1)
           call warn(myname, msg)
           i = 1
       endif

       if(xWant > x(n)) then 
           write(msg,*) "in polint, wanted: ", xWant, ", got upper bound: ", x(n)
           call warn(myname, msg)
           i = n
       endif

       ! if i is still zero find i less than xWant
       if(i == 0) then
           do j = 1, n
               if(xWant >= x(j)) i = j
           enddo
       endif

       ! slope
       if(i == n) then 
           slope = 0.0
       else
           slope = (y(i+1)-y(i)) / (x(i+1)-x(i))
       endif

       dx = xWant - x(i)
       yWant = y(i) + slope*dx

       yErr = 0.0

       return
   end subroutine polint

 end module MAML_OpticsTableMod

