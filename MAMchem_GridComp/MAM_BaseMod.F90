#include "MAPL_Exceptions.h"
#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: MAM_BaseMod - basic MAM constants and types
!
! !INTERFACE:
!
   module MAM_BaseMod
!
! !USES:
!

   use MAPL_Mod
   use MAPL_SimpleBundleMod

   use MAM_ComponentsDataMod
   use MAM_ConstituentsDataMod
 
   use MAM3_DataMod
!  use MAM4_DataMod
   use MAM7_DataMod

   implicit NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:

!
! !DESCRIPTION: 
!
!  {\tt MAM\_Base} provides a collection of basic constants and 
!  types used in the MAM code.
!
!
! !REVISION HISTORY:
!
!  14Sep2011  A. Darmenov  Initial version
!
!EOP
!-------------------------------------------------------------------------



!  MAM public types
!  ------------------
   public MAM_Range
   public MAM_AerosolSpecies
   public MAM_AerosolMode

   public MAM_Scheme


!  MAM public methods
!  ------------------------
   public MAM_RangeSet
   public MAM_RangeGet

   public MAM_AerosolComponentSet
   public MAM_AerosolComponentGet

   public MAM_AerosolSpeciesGet

   public MAM_AerosolModeSet
   public MAM_AerosolModeGet

   public MAM_SchemeGetModeIndex

   public MAM_SchemeInit

   public MAM_SchemeValidate


!  Status codes
!  ------------
   integer, public, parameter :: MAM_SUCCESS                           = 0
   integer, public, parameter :: MAM_GENERAL_ERROR                     = 9000 + 1
   integer, public, parameter :: MAM_NOT_IMPLEMENTED_ERROR             = 9000 + 2
   integer, public, parameter :: MAM_UNKNOWN_SCHEME_ERROR              = 9000 + 3
   integer, public, parameter :: MAM_UNKNOWN_AEROSOL_MODE_ERROR        = 9000 + 4  
   integer, public, parameter :: MAM_UNKNOWN_AEROSOL_COMPONENT_ERROR   = 9000 + 5
   integer, public, parameter :: MAM_UNKNOWN_AEROSOL_CONSTITUENT_ERROR = 9000 + 6


!  MAM defs
!  ------------------
   integer, public, parameter :: MAM3_SCHEME               = 1000 + 3
   integer, public, parameter :: MAM4_SCHEME               = 1000 + 4
   integer, public, parameter :: MAM7_SCHEME               = 1000 + 7

   integer, public, parameter :: MAM_INTERSTITIAL_AEROSOL  = 5000 + 1
   integer, public, parameter :: MAM_CLOUDBORNE_AEROSOL    = 5000 + 2

   integer, public, parameter :: MAM_NUMBER_AEROSOL        = 6000 + 1
   integer, public, parameter :: MAM_MASS_AEROSOL          = 6000 + 2



   integer, public, parameter :: MAM_MAX_NUMBER_MODES = max(MAM3_MODES, MAM7_MODES)

   integer, public, parameter :: MAM_MAXSTR = 256                  ! maximum string length

   



!  A range/interval 
!  ----------------
   type MAM_Range
       real                                :: low = 0.0           ! lower bound of the range
       real                                :: up  = 0.0           ! upper bound of the range
   end type


!  Aerosol Component
!  -----------------
   type MAM_AerosolComponent
       character(len=MAM_MAXSTR)           :: name                ! name of the aerosol component
       character(len=MAM_MAXSTR)           :: long_name           ! long name

       real                                :: density             ! bulk density 
       real                                :: hygroscopicity      ! hygroscopicity of the aerosol component
       real                                :: solubility          ! fractional solubility
   end type


!  Aerosol Species
!  ---------------
   type MAM_AerosolSpecies
       character(len=MAM_MAXSTR)           :: name                ! name of the aerosol specie
       character(len=MAM_MAXSTR)           :: long_name           ! long name

       type(MAM_AerosolComponent), pointer :: component => null() ! aerosol component

       type(MAM_Range), pointer            :: emission_size_range => null()  ! cutoff size range for emissions
   end type



!  Aerosol Mode
!  ------------
   type MAM_AerosolMode
       character(len=MAM_MAXSTR)           :: name                 ! name of the aerosol mode
       character(len=MAM_MAXSTR)           :: long_name            ! long name 

       real                                :: sigma                ! geometric standard deviation

       real                                :: size                 ! default geometric mean diameter of number size distribution
       type(MAM_Range)                     :: size_range           ! size range limits

       real                                :: rh_deliquescence     ! deliquescence   RH point
       real                                :: rh_crystallization   ! crystallization RH point

       real                                :: f_conv_scav          ! convective scavenging parameter
       real                                :: f_wet                ! wet removal efficiency

       integer                             :: n_species            ! number of species
       type(MAM_AerosolSpecies), pointer   :: species(:) => null() ! species
   end type


!  Gas Species
!  -----------
   type MAM_GasSpecies
       character(len=MAM_MAXSTR)           :: name                ! name of the gas specie
       character(len=MAM_MAXSTR)           :: long_name           ! long name
       character(len=MAM_MAXSTR)           :: units               ! units

       real                                :: mw                  ! molecular weight
   end type


!  MAM Tracer
!  ----------
   type MAM_Tracer
       integer                             :: id                  ! tracer ID
       character(len=MAM_MAXSTR)           :: name                ! name of the tracer
       character(len=MAM_MAXSTR)           :: long_name           ! long name

       character(len=MAM_MAXSTR)           :: units               ! units
       integer                             :: type                ! interstital aerosol, cloud-borne aerosol or gas

       type(MAM_AerosolSpecies), pointer   :: aerosol_species => null() 
       type(MAM_GasSpecies), pointer       :: gas_species => null() 

!      real, pointer, dimension(:,:,:)     :: q => null()
   end type


!  Modal Aerosol Model Scheme/Configuration 
!  ----------------------------------------
   type MAM_Scheme
       integer                             :: id                   ! scheme ID
       character(len=MAM_MAXSTR)           :: name                 ! name of the scheme
       character(len=MAM_MAXSTR)           :: long_name            ! long name of the scheme


       ! Aerosol microphysics
       ! --------------------
       integer                             :: n_aerosol_components ! number of aerosol components
       type(MAM_AerosolComponent), pointer :: aerosol_component(:) => null()

       integer                             :: n_modes              ! number of (aerosol) modes
       type(MAM_AerosolMode), pointer      :: mode(:) => null()
   end type


   interface MAM_SchemeGetAerosolComponentIndex
       module procedure MAM_SchemeGetAerosolComponentIndexFromName
   end interface MAM_SchemeGetAerosolComponentIndex


   interface MAM_SchemeGetModeIndex
       module procedure MAM_SchemeGetModeIndexFromName
   end interface MAM_SchemeGetModeIndex


   contains


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_RangeSet --- Sets range bounds
!
! !INTERFACE:

   subroutine MAM_RangeSet(self, lower_bound, upper_bound, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_Range), intent(inout) :: self

! !INPUT PARAMETERS:   
       real, optional, intent(in)     :: lower_bound
       real, optional, intent(in)     :: upper_bound

! !OUTPUT PARAMETERS:
       integer, optional, intent(out) :: rc 


! !DESCRIPTION: Sets upper and lower bounds of a range. 
!
! !REVISION HISTORY:
!
!  15Sep2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                        __Iam__('MAM_RangeGet')

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(lower_bound)) &
           self%low = lower_bound

       if (present(upper_bound)) &    
           self%up  = upper_bound

       if (present(rc)) rc = MAM_SUCCESS

   end subroutine MAM_RangeSet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_RangeGet --- Gets bounds of a range
!
! !INTERFACE:

   subroutine MAM_RangeGet(self, lower_bound, upper_bound, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       real, optional, intent(inout)  :: lower_bound
       real, optional, intent(inout)  :: upper_bound

! !INPUT PARAMETERS:
       type(MAM_Range), intent(in)    :: self

! !OUTPUT PARAMETERS:
       integer, optional, intent(out) :: rc

! !DESCRIPTION: Gets the upper and lower bounds of a size range/bin. 
!
! !REVISION HISTORY:
!
!  15Sep2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('MAM_RangeGet')


       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(lower_bound)) lower_bound = self%low
 
       if (present(upper_bound)) upper_bound = self%up

       if (present(rc)) rc = MAM_SUCCESS

   end subroutine MAM_RangeGet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_RangeGet --- Gets bounds of a range
!
! !INTERFACE:

   subroutine MAM_RangePrint(self, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
       type(MAM_Range), intent(in)    :: self

! !OUTPUT PARAMETERS:
       integer, optional, intent(out) :: rc


! !DESCRIPTION: Gets the upper and lower bounds of a size range/bin. 
!
! !REVISION HISTORY:
!
!  15Sep2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('MAM_RangePrint')


       real :: lower_bound
       real :: upper_bound

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       call MAM_RangeGet(self, lower_bound = lower_bound, &
                               upper_bound = upper_bound, __RC__)

       write (*, *) lower_bound, upper_bound

       if (present(rc)) rc = STATUS

   end subroutine MAM_RangePrint


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolComponentSet ---
!
! !INTERFACE:

   subroutine MAM_AerosolComponentSet(self, name,           &
                                            long_name,      &
                                            density,        &
                                            hygroscopicity, &
                                            solubility,     & 
                                            rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_AerosolComponent), intent(inout) :: self

! !INPUT PARAMETERS:   
       character(len=*), optional, intent(in)    :: name
       character(len=*), optional, intent(in)    :: long_name

       real,             optional, intent(in)    :: density
       real,             optional, intent(in)    :: hygroscopicity
       real,             optional, intent(in)    :: solubility

! !OUTPUT PARAMETERS:
       integer,          optional, intent(out)   :: rc

! !DESCRIPTION: Sets aerosol component properties
!
! !REVISION HISTORY:
!
!  15Sep2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('MAM_AerosolComponentSet')

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(name)) &
           self%name = trim(name)

       if (present(long_name)) &
           self%long_name = trim(long_name)
   
       if (present(density)) &
           self%density = density
   
       if (present(hygroscopicity)) &
           self%hygroscopicity = hygroscopicity

       if (present(solubility)) &
           self%solubility = solubility

       if (present(rc)) rc = MAM_SUCCESS

   end subroutine MAM_AerosolComponentSet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolComponentGet ---
!
! !INTERFACE:

   subroutine MAM_AerosolComponentGet(self, name,           &
                                            long_name,      &
                                            density,        &
                                            hygroscopicity, &
                                            solubility,     &
                                            rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       character(len=*), optional, intent(inout) :: name
       character(len=*), optional, intent(inout) :: long_name

       real,             optional, intent(inout) :: density
       real,             optional, intent(inout) :: hygroscopicity
       real,             optional, intent(inout) :: solubility

! !INPUT PARAMETERS:
       type(MAM_AerosolComponent), intent(in)    :: self

! !OUTPUT PARAMETERS:
       integer,          optional, intent(out)   :: rc 

! !DESCRIPTION: Get aerosol component properties
!
! !REVISION HISTORY:
!
!  15Sep2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('MAM_AerosolComponentGet')

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(name))&
           name = trim(self%name)

       if (present(long_name)) &
           long_name = trim(self%long_name)
   
       if (present(density)) &
           density = self%density
   
       if (present(hygroscopicity)) &
           hygroscopicity = self%hygroscopicity

       if (present(solubility)) &
           solubility = self%solubility

       if (present(rc)) rc = MAM_SUCCESS

   end subroutine MAM_AerosolComponentGet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolComponentPrint ---
!
! !INTERFACE:

   subroutine MAM_AerosolComponentPrint(self, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
       type(MAM_AerosolComponent), intent(in)  :: self

! !OUTPUT PARAMETERS:
       integer,          optional, intent(out) :: rc 

! !DESCRIPTION: Get aerosol component properties
!
! !REVISION HISTORY:
!
!  15Sep2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('MAM_AerosolComponentPrint')

       character(len=MAM_MAXSTR) :: name
       character(len=MAM_MAXSTR) :: long_name

       real                      :: density
       real                      :: hygroscopicity
       real                      :: solubility


       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       call MAM_AerosolComponentGet(self, name           = name,           &
                                          long_name      = long_name,      &
                                          density        = density,        &
                                          hygroscopicity = hygroscopicity, &
                                          solubility     = solubility, __RC__)

       write (*, *) 'Aerosol Component: '
       write (*, *) '    name           = ', trim(name)
       write (*, *) '    long_name      = ', trim(long_name) 
       write (*, *) '    density        = ', density
       write (*, *) '    hygroscopicity = ', hygroscopicity
       write (*, *) '    solubility     = ', solubility
       write (*, *)

       if (present(rc)) rc = STATUS

   end subroutine MAM_AerosolComponentPrint
   

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolSpeciesSet ---
!
! !INTERFACE:

   subroutine MAM_AerosolSpeciesSet(self, name,                &
                                          long_name,           &
                                          component,           &
                                          emission_size_range, &
                                          rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_AerosolSpecies), intent(inout) :: self

! !INPUT PARAMETERS:   
       character(len=*), optional, intent(in)  :: name
       character(len=*), optional, intent(in)  :: long_name

       type(MAM_Range),  optional, intent(in)  :: emission_size_range

       type(MAM_AerosolComponent), pointer, optional, intent(in) :: component

! !OUTPUT PARAMETERS:
       integer, optional, intent(out)         :: rc

! !DESCRIPTION: Sets aerosol species properties
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                  __Iam__('MAM_AerosolSpeciesSet')

    
       type(MAM_Range), pointer :: buff

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(name))&
           self%name = trim(name)

       if (present(long_name)) &
           self%long_name = trim(long_name)
   
       if (present(component)) &
           self%component => component
   
       if (present(emission_size_range)) then

           if (.not. associated(self%emission_size_range)) then
               allocate(buff, __STAT__)
           end if

           buff = emission_size_range
           self%emission_size_range => buff

           self%emission_size_range = emission_size_range
       end if

       if (present(rc)) rc = MAM_SUCCESS

   end subroutine MAM_AerosolSpeciesSet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolSpeciesGet ---
!
! !INTERFACE:

   subroutine MAM_AerosolSpeciesGet(self, name,                &
                                          long_name,           &
                                          component,           &
                                          emission_size_range, &
                                          density,             &
                                          hygroscopicity,      &
                                          solubility,          &
                                          component_name,      &
                                          component_long_name, &
                                          rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       character(len=*), optional, intent(inout) :: name
       character(len=*), optional, intent(inout) :: long_name

       type(MAM_Range),  optional, intent(inout) :: emission_size_range

       type(MAM_AerosolComponent), pointer, optional, intent(inout) :: component

       real,             optional, intent(inout) :: density
       real,             optional, intent(inout) :: hygroscopicity
       real,             optional, intent(inout) :: solubility
       character(len=*), optional, intent(inout) :: component_name
       character(len=*), optional, intent(inout) :: component_long_name

! !INPUT PARAMETERS:
       type(MAM_AerosolSpecies), intent(in)      :: self

! !OUTPUT PARAMETERS:
       integer,          optional, intent(out)   :: rc

! !DESCRIPTION: Sets aerosol species properties
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                  __Iam__('MAM_AerosolSpeciesGet')

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR
    
       if (present(name))&
           name = trim(self%name)

       if (present(long_name)) &
           long_name = trim(self%long_name)
   
       if (present(component)) &
           component => self%component
   
       if (present(emission_size_range)) then
           if (associated(self%emission_size_range)) then
               emission_size_range = self%emission_size_range
           else
               call MAM_RangeSet(emission_size_range, -MAPL_UNDEF, -MAPL_UNDEF, __RC__)
           end if
       end if

       if(present(density)) then
           if (associated(self%component)) then
               call MAM_AerosolComponentGet(self%component, density=density, __RC__)
           else
               density = MAPL_UNDEF
           end if
       end if

       if (present(hygroscopicity)) then
           if (associated(self%component)) then
               call MAM_AerosolComponentGet(self%component, hygroscopicity=hygroscopicity, __RC__)
           else
               hygroscopicity = MAPL_UNDEF
           end if
       end if

       if (present(solubility)) then
           if (associated(self%component)) then
               call MAM_AerosolComponentGet(self%component, solubility=solubility, __RC__)
           else
               solubility = MAPL_UNDEF
           end if
       end if

       if (present(component_name)) then
           if (associated(self%component)) then
               call MAM_AerosolComponentGet(self%component, name=component_name, __RC__)
           else
               component_name = ''
           end if
       end if

       if (present(component_long_name)) then
           if (associated(self%component)) then
               call MAM_AerosolComponentGet(self%component, long_name=component_long_name, __RC__)
           else
               component_long_name = ''
           end if
       end if

       if (present(rc)) rc = STATUS

   end subroutine MAM_AerosolSpeciesGet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolSpeciesPrint ---
!
! !INTERFACE:

   subroutine MAM_AerosolSpeciesPrint(self, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
       type(MAM_AerosolSpecies), intent(in) :: self

! !OUTPUT PARAMETERS:
       integer, optional, intent(out)       :: rc

! !DESCRIPTION: Sets aerosol species properties
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                  __Iam__('MAM_AerosolSpeciesPrint')

       character(len=MAM_MAXSTR) :: name
       character(len=MAM_MAXSTR) :: long_name

       type(MAM_Range)           :: emission_size_range

       character(len=MAM_MAXSTR) :: component_name


       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       call MAM_AerosolSpeciesGet(self, name = name, &
                                        long_name = long_name, &
                                        emission_size_range = emission_size_range, &
                                        component_name = component_name, __RC__)
       
       write (*, *) 'Aerosol Species: '
       write (*, *) '    name           = ', trim(name)
       write (*, *) '    long_name      = ', trim(long_name) 
       write (*, *) '    emission_range = ', emission_size_range%low, emission_size_range%up
       write (*, *) '    component      = ', component_name
       write (*, *)

       if (present(rc)) rc = STATUS

   end subroutine MAM_AerosolSpeciesPrint


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolModeSet ---
!
! !INTERFACE:

   subroutine MAM_AerosolModeSet(self, name,               &
                                       long_name,          &
                                       sigma,              &
                                       size_default,       &
                                       size_min,           &
                                       size_max,           &
                                       rh_deliquescence,   &
                                       rh_crystallization, &
                                       f_conv_scavenging,  &
                                       f_wet,              &
                                       n_species,          &
                                       rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_AerosolMode), intent(inout)   :: self

! !INPUT PARAMETERS:   
       character(len=*), optional, intent(in)  :: name
       character(len=*), optional, intent(in)  :: long_name

       real,             optional, intent(in)  :: sigma

       real,             optional, intent(in)  :: size_default
       real,             optional, intent(in)  :: size_min
       real,             optional, intent(in)  :: size_max

       real,             optional, intent(in)  :: rh_deliquescence
       real,             optional, intent(in)  :: rh_crystallization

       real,             optional, intent(in)  :: f_conv_scavenging
       real,             optional, intent(in)  :: f_wet 

       integer,          optional, intent(in)  :: n_species

! !OUTPUT PARAMETERS:
       integer,          optional, intent(out) :: rc

! !DESCRIPTION: Sets aerosol species properties
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('MAM_AerosolModeSet')


       type(MAM_AerosolSpecies), pointer, dimension(:) :: buff

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(name))&
           self%name = trim(name)

       if (present(long_name)) &
           self%long_name = trim(long_name)

       if (present(sigma)) &
           self%sigma = sigma

       if (present(size_default)) &
           self%size = size_default
   
       if (present(size_min)) &
           self%size_range%low = size_min

       if (present(size_max)) &
           self%size_range%up = size_max

        if (present(rh_deliquescence)) &
           self%rh_deliquescence = rh_deliquescence

       if (present(rh_crystallization)) &
           self%rh_crystallization = rh_crystallization

       if (present(f_conv_scavenging)) &
           self%f_conv_scav = f_conv_scavenging

       if (present(f_wet)) &
           self%f_wet = f_wet

       if (present(n_species)) then
           if (n_species > 0) then
               allocate(buff(n_species), __STAT__)

               self%n_species = n_species
               self%species => buff
           else 
               __raise__(MAM_GENERAL_ERROR, 'Number of aerosol species must be positive.')      
           end if
       end if

       if (present(rc)) then
           rc = MAM_SUCCESS
       end if

   end subroutine MAM_AerosolModeSet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolModeGet ---
!
! !INTERFACE:

   subroutine MAM_AerosolModeGet(self, name,               &
                                       long_name,          &
                                       sigma,              &
                                       size_default,       &
                                       size_min,           &
                                       size_max,           &
                                       rh_deliquescence,   &
                                       rh_crystallization, &
                                       f_conv_scavenging,  &
                                       f_wet,              &
                                       n_species,          &
                                       species,            &
                                       rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       character(len=*), optional, intent(inout)   :: name
       character(len=*), optional, intent(inout)   :: long_name

       real,             optional, intent(inout)   :: sigma

       real,             optional, intent(inout)   :: size_default
       real,             optional, intent(inout)   :: size_min
       real,             optional, intent(inout)   :: size_max

       real,             optional, intent(inout)   :: rh_deliquescence
       real,             optional, intent(inout)   :: rh_crystallization

       real,             optional, intent(inout)   :: f_conv_scavenging
       real,             optional, intent(inout)   :: f_wet

       integer,          optional, intent(inout)   :: n_species
       type(MAM_AerosolSpecies), optional, pointer :: species(:)


! !INPUT PARAMETERS:   
       type(MAM_AerosolMode),      intent(in)      :: self

! !OUTPUT PARAMETERS:
       integer,          optional, intent(out)     :: rc

! !DESCRIPTION: Sets aerosol species properties
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('MAM_AerosolModeGet')

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR
        
       if (present(name))&
           name = trim(self%name)

       if (present(long_name)) &
           long_name = trim(self%long_name)

       if (present(sigma)) &
           sigma = self%sigma

       if (present(size_default)) &
           size_default = self%size
   
       if (present(size_min)) &
           size_min = self%size_range%low

       if (present(size_max)) &
           size_max = self%size_range%up

       if (present(rh_deliquescence)) &
           rh_deliquescence = self%rh_deliquescence

       if (present(rh_crystallization)) &
           rh_crystallization = self%rh_crystallization

       if (present(f_conv_scavenging)) &
           f_conv_scavenging = self%f_conv_scav

       if (present(f_wet)) &
           f_wet = self%f_wet

       if (present(n_species)) &
           n_species = self%n_species

       if (present(species)) &
           species => self%species

       if (present(rc)) rc = MAM_SUCCESS
       
   end subroutine MAM_AerosolModeGet


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_AerosolModePrint ---
!
! !INTERFACE:

   subroutine MAM_AerosolModePrint(self, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:   
       type(MAM_AerosolMode), intent(in) :: self

! !OUTPUT PARAMETERS:
       integer, optional, intent(out)    :: rc

! !DESCRIPTION: Sets aerosol species properties
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                   __Iam__('MAM_AerosolModeGet')

       character(len=MAM_MAXSTR)        :: name
       character(len=MAM_MAXSTR)        :: long_name

       real                             :: sigma

       real                             :: size_default
       real                             :: size_min
       real                             :: size_max

       real                             :: rh_deliquescence
       real                             :: rh_crystallization 

       integer                          :: n, n_species


       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       call MAM_AerosolModeGet(self, name               = name,               &
                                     long_name          = long_name,          &
                                     sigma              = sigma,              &
                                     size_default       = size_default,       &
                                     size_min           = size_min,           & 
                                     size_max           = size_max,           &
                                     rh_deliquescence   = rh_deliquescence,   &
                                     rh_crystallization = rh_crystallization, &
                                     n_species          = n_species, __RC__)


       write (*, *) 'Aerosol Mode: '
       write (*, *) '    name           = ', trim(name)
       write (*, *) '    long_name      = ', trim(long_name) 
       write (*, *) '    sigma          = ', sigma
       write (*, *) '    default size   = ', size_default
       write (*, *) '    min/max size   = ', size_min, size_max
       write (*, *) '    RH deliquesc.  = ', rh_deliquescence
       write (*, *) '    RH crystalliz. = ', rh_crystallization
       write (*, *) '    species        = ', n_species
       write (*, *)  
       
       do n = 1, n_species
           call MAM_AerosolSpeciesPrint(self%species(n), __RC__)
       end do

       write (*, *)

       if (present(rc)) rc = MAM_SUCCESS
       
   end subroutine MAM_AerosolModePrint


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_SchemeValidate --- 
!
! !INTERFACE:

   subroutine MAM_SchemeValidate(self, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_Scheme), intent(inout) :: self        ! MAM scheme/configuration

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:
       integer, optional, intent(out)  :: rc          ! return code 

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  03Dec2013  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                     __Iam__('MAM_SchemeValidate')
  

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       select case(self%id)

           case (MAM3_SCHEME)
               STATUS = MAM_NOT_IMPLEMENTED_ERROR
               __raise__(MAM_NOT_IMPLEMENTED_ERROR, 'MAM3 scheme is not implemented yet.')

           case (MAM7_SCHEME)
               STATUS = MAM_SUCCESS

           case default
               STATUS = MAM_UNKNOWN_SCHEME_ERROR
                __raise__(MAM_UNKNOWN_SCHEME_ERROR, 'Unsupported MAM configuration.')

       end select
       
       if (present(rc)) rc = STATUS
       
   end subroutine MAM_SchemeValidate


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_SchemeInit --- 
!
! !INTERFACE:

   subroutine MAM_SchemeInit(self, scheme_id, verbose, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_Scheme), intent(inout) :: self        ! MAM scheme/configuration
       logical, optional, intent(in)   :: verbose     ! verbosity flag

! !INPUT PARAMETERS:
       integer, intent(in)             :: scheme_id   ! MAM model scheme -- MAM3/MAM7

! !OUTPUT PARAMETERS:
       integer, optional, intent(out)  :: rc          ! return code

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------

                     __Iam__('MAM_SchemeInit')
   
       logical :: verbose_

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(verbose)) then
           verbose_ = verbose
       else
           verbose_ = .False.
       end if


       select case(scheme_id)
           case (MAM3_SCHEME)
               call MAM3_Init(self, verbose_, __RC__)

           case (MAM7_SCHEME)
               call MAM7_Init(self, verbose_, __RC__)

           case default
               __raise__ (MAM_UNKNOWN_SCHEME_ERROR, 'Unsupported MAM model.')
       end select
      
     if (verbose_ .and. MAPL_AM_I_ROOT()) then
         call MAM_SchemePrint(self, __RC__)
     end if

     if (present(rc)) rc = STATUS

   end subroutine MAM_SchemeInit


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM7_Init --- 
!
! !INTERFACE:

   subroutine MAM7_Init(self, verbose, rc)
! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_Scheme),  intent(inout) :: self        ! MAM scheme/configuration
       logical, optional, intent(in)    :: verbose     ! verbosity flag
       integer, optional, intent(inout) :: rc          ! return code

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Set up MAM7 machinery
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                 __Iam__('MAM7_Init')
   

       type(MAM_AerosolComponent), pointer :: aerosol_component
       type(MAM_AerosolSpecies),   pointer :: aerosol_species
       type(MAM_AerosolMode),      pointer :: mode

       type(MAM_Range)                     :: emission_range
        
       logical :: verbose_
       integer :: m, s, c


       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(verbose)) then
           verbose_ = verbose
       else
           verbose_ = .false.
       end if    

       ! set the model scheme ID and name     
       self%id        = MAM7_SCHEME
       self%name      = 'MAM7'
       self%long_name = 'GEOS5/MAM7 Aerosol Model'

       ! initialize the aerosol components
       call MAM7_AerosolComponentsInit(self, verbose_, __RC__)


       ! initialize the aerosol modes
       call MAM7_AerosolModesInit(self, verbose_, __RC__)
   

       ! populate AITKEN
       ! ---------------
       m = MAM_SchemeGetModeIndex(self, MAM7_AITKEN_MODE_NAME, __RC__)
       mode => self%mode(m)

       s = 0
       
       ! add SULFATE
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - sulfate'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SULFATE_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name     = MAM_SULFATE_CONSTITUENT_NAME, &
                                                   long_name = 'sulfate aerosol species',    &
                                                   component = aerosol_component, __RC__)
       
       ! add AMMONIUM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - ammonium'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_AMMONIUM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_AMMONIUM_CONSTITUENT_NAME, &
                                                   long_name = 'ammonium aerosol species',    &
                                                   component = aerosol_component, __RC__)
       
       ! add SOA
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - SOA'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SOA_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SOA_CONSTITUENT_NAME, &
                                                   long_name = 'SOA aerosol species',    &
                                                   component = aerosol_component, __RC__)
       
       ! add SEASALT
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - seasalt'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SEASALT_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_RangeSet(emission_range, lower_bound = MAM7_AIT_SS_D_CUTOFF(1), &
                                         upper_bound = MAM7_AIT_SS_D_CUTOFF(2), __RC__)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SEASALT_CONSTITUENT_NAME, &
                                                   long_name = 'seasalt aerosol species',    &
                                                   component = aerosol_component,            &
                                                   emission_size_range = emission_range, __RC__)
                                                 
       ! check if the right number of species were added
       if ( s /= MAM7_AITKEN_MODE_SPECIES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol species in Aitken mode.")
       end if 


       ! populate ACCUMULATION
       ! ---------------------
       m = MAM_SchemeGetModeIndex(self, MAM7_ACCUMULATION_MODE_NAME, __RC__)
       mode => self%mode(m)

       s = 0
       
       ! add SULFATE
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - sulfate'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SULFATE_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SULFATE_CONSTITUENT_NAME, &
                                                   long_name = 'sulfate aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! add AMMONIUM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - ammonium'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_AMMONIUM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_AMMONIUM_CONSTITUENT_NAME, &
                                                   long_name = 'ammonium aerosol species',    &
                                                   component = aerosol_component, __RC__)
       
       ! add SOA
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - SOA'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SOA_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SOA_CONSTITUENT_NAME, &
                                                   long_name = 'SOA aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! add POM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - POM'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_POM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_POM_CONSTITUENT_NAME, &
                                                   long_name = 'POM aerosol species',    &
                                                   component = aerosol_component, __RC__)                                          

       ! add BLACK CARBON
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - black carbon'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_BLACK_CARBON_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_BLACK_CARBON_CONSTITUENT_NAME, &
                                                   long_name = 'black carbon aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! add SEASALT
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - seasalt'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SEASALT_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_RangeSet(emission_range, lower_bound = MAM7_ACC_SS_D_CUTOFF(1), &
                                         upper_bound = MAM7_ACC_SS_D_CUTOFF(2), __RC__)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SEASALT_CONSTITUENT_NAME, &
                                                   long_name = 'seasalt aerosol species',     &
                                                   component = aerosol_component,            &
                                                   emission_size_range = emission_range, __RC__)
      
       ! check if the right number of species were added
       if ( s /= MAM7_ACCUMULATION_MODE_SPECIES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol species in Accumulation mode.")
       end if

      
       ! populate PRIMARY CARBON MODE
       ! ----------------------------
       m = MAM_SchemeGetModeIndex(self, MAM7_PRIMARY_CARBON_MODE_NAME, __RC__)
       mode => self%mode(m)

       s = 0
       
       ! add POM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - POM'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_POM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_POM_CONSTITUENT_NAME, &
                                                   long_name = 'POM aerosol species',    &
                                                   component = aerosol_component, __RC__)                                          

       ! add BLACK CARBON
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - black carbon'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_BLACK_CARBON_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_BLACK_CARBON_CONSTITUENT_NAME, &
                                                   long_name = 'black carbon aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! check if the right number of species were added
       if ( s /= MAM7_PRIMARY_CARBON_MODE_SPECIES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol species in Primary organic mode.")
       end if


       ! populate FINE SEASALT
       ! ---------------------
       m = MAM_SchemeGetModeIndex(self, MAM7_FINE_SEASALT_MODE_NAME, __RC__)
       mode => self%mode(m)

       s = 0
       
       ! add SULFATE
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - sulfate'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SULFATE_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SULFATE_CONSTITUENT_NAME, &
                                                   long_name = 'sulfate aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! add AMMONIUM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - ammonium'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_AMMONIUM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_AMMONIUM_CONSTITUENT_NAME, &
                                                   long_name = 'ammonium aerosol species',     &
                                                   component = aerosol_component, __RC__)
       
       ! add SEASALT
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - seasalt'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SEASALT_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_RangeSet(emission_range, lower_bound = MAM7_FSS_SS_D_CUTOFF(1), &
                                         upper_bound = MAM7_FSS_SS_D_CUTOFF(2), __RC__)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SEASALT_CONSTITUENT_NAME, &
                                                   long_name = 'seasalt aerosol species',    &
                                                   component = aerosol_component,            &
                                                   emission_size_range = emission_range, __RC__)
      
       ! check if the right number of species were added
       if ( s /= MAM7_FINE_SEASALT_MODE_SPECIES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol species in Fine seasalt mode.")
       end if
 

       ! populate FINE DUST
       ! ------------------
       m = MAM_SchemeGetModeIndex(self, MAM7_FINE_DUST_MODE_NAME, __RC__)
       mode => self%mode(m)

       s = 0
       
       ! add SULFATE
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - sulfate'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SULFATE_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SULFATE_CONSTITUENT_NAME, &
                                                   long_name = 'sulfate aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! add AMMONIUM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - ammonium'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_AMMONIUM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_AMMONIUM_CONSTITUENT_NAME, &
                                                   long_name = 'ammonium aerosol species',     &
                                                   component = aerosol_component, __RC__)
       
       ! add DUST
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - dust'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_DUST_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_RangeSet(emission_range, lower_bound = MAM7_FDU_DU_D_CUTOFF(1), &
                                         upper_bound = MAM7_FDU_DU_D_CUTOFF(2), __RC__)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_DUST_CONSTITUENT_NAME, &
                                                   long_name = 'dust aerosol species',    &
                                                   component = aerosol_component,            &
                                                   emission_size_range = emission_range, __RC__)
      
       ! check if the right number of species were added
       if ( s /= MAM7_FINE_DUST_MODE_SPECIES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol species in Fine dust mode.")
       end if


       ! populate COARSE SEASALT
       ! ---------------------
       m = MAM_SchemeGetModeIndex(self, MAM7_COARSE_SEASALT_MODE_NAME, __RC__)
       mode => self%mode(m)

       s = 0
       
       ! add SULFATE
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - sulfate'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SULFATE_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SULFATE_CONSTITUENT_NAME, &
                                                   long_name = 'sulfate aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! add AMMONIUM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - ammonium'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_AMMONIUM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_AMMONIUM_CONSTITUENT_NAME, &
                                                   long_name = 'ammonium aerosol species',     &
                                                   component = aerosol_component, __RC__)
       
       ! add SEASALT
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - seasalt'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SEASALT_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_RangeSet(emission_range, lower_bound = MAM7_CSS_SS_D_CUTOFF(1), &
                                         upper_bound = MAM7_CSS_SS_D_CUTOFF(2), __RC__)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SEASALT_CONSTITUENT_NAME, &
                                                   long_name = 'seasalt aerosol species',    &
                                                   component = aerosol_component,            &
                                                   emission_size_range = emission_range, __RC__)
      
       ! check if the right number of species were added
       if ( s /= MAM7_COARSE_SEASALT_MODE_SPECIES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol species in Coarse seasalt mode.")
       end if
 

       ! populate COARSE DUST
       ! ------------------
       m = MAM_SchemeGetModeIndex(self, MAM7_COARSE_DUST_MODE_NAME, __RC__)
       mode => self%mode(m)

       s = 0
       
       ! add SULFATE
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - sulfate'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_SULFATE_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_SULFATE_CONSTITUENT_NAME, &
                                                   long_name = 'sulfate aerosol species',    &
                                                   component = aerosol_component, __RC__)

       ! add AMMONIUM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - ammonium'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_AMMONIUM_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_AMMONIUM_CONSTITUENT_NAME, &
                                                   long_name = 'ammonium aerosol species',    &
                                                   component = aerosol_component, __RC__)
       
       ! add DUST
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding species - dust'
       s = s + 1
       c = MAM_SchemeGetAerosolComponentIndex(self, MAM_DUST_COMPONENT_NAME, __RC__)
       aerosol_component => self%aerosol_component(c)
       call MAM_RangeSet(emission_range, lower_bound = MAM7_CDU_DU_D_CUTOFF(1), &
                                         upper_bound = MAM7_CDU_DU_D_CUTOFF(2), __RC__)
       call MAM_AerosolSpeciesSet(mode%species(s), name      = MAM_DUST_CONSTITUENT_NAME, &
                                                   long_name = 'dust aerosol species',    &
                                                   component = aerosol_component,         &
                                                   emission_size_range = emission_range, __RC__)
      
       ! check if the right number of species were added
       if ( s /= MAM7_COARSE_DUST_MODE_SPECIES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol species in Coarse dust mode.")
       end if


       if (present(rc)) rc = STATUS

   end subroutine MAM7_Init


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM7_AerosolComponentsInit --- 
!
! !INTERFACE:

   subroutine MAM7_AerosolComponentsInit(self, verbose, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_Scheme),  intent(inout) :: self        ! MAM scheme/configuration
       logical, optional, intent(in)    :: verbose     ! verbosity flag
       integer, optional, intent(inout) :: rc          ! return code

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:
   

! !DESCRIPTION: Set up MAM7 machinery
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                 __Iam__('MAM7_AerosolComponentsInit')
   

       type(MAM_AerosolComponent), pointer :: aero_component

       logical :: verbose_
       integer :: n

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(verbose)) then
           verbose_ = verbose
       else
           verbose_ = .false.
       end if

       ! allocate memory for the aerosol components
       allocate(self%aerosol_component(MAM7_AEROSOL_COMPONENTS), __STAT__)

       
       ! initialize the counter
       n = 0

       ! SULFATE
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding component - sulfate'
       n = n + 1
       aero_component => self%aerosol_component(n)
       call MAM_AerosolComponentSet(aero_component,                                                &
                                    name           = trim(MAM_SULFATE_COMPONENT_NAME),             &
                                    long_name      = trim(MAM_SULFATE_COMPONENT_NAME)//' aerosol', &
                                    density        = MAM_SULFATE_COMPONENT_DENSITY,                &
                                    hygroscopicity = MAM_SULFATE_COMPONENT_HYGROSCOPICITY,         &
                                    solubility     = MAM_SULFATE_COMPONENT_SOLUBILITY, __RC__)

       ! AMMONIUM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding component - ammonium'
       n = n + 1
       aero_component => self%aerosol_component(n)
       call MAM_AerosolComponentSet(aero_component,                                                 &
                                    name           = trim(MAM_AMMONIUM_COMPONENT_NAME),             &
                                    long_name      = trim(MAM_AMMONIUM_COMPONENT_NAME)//' aerosol', &
                                    density        = MAM_AMMONIUM_COMPONENT_DENSITY,                &
                                    hygroscopicity = MAM_AMMONIUM_COMPONENT_HYGROSCOPICITY,         &
                                    solubility     = MAM_AMMONIUM_COMPONENT_SOLUBILITY, __RC__)
       
       ! BLACK CARBON
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding component - black carbon'
       n = n + 1
       aero_component => self%aerosol_component(n)
       call MAM_AerosolComponentSet(aero_component,                                                &
                                    name = MAM_BLACK_CARBON_COMPONENT_NAME,                        &
                                    long_name = trim(MAM_BLACK_CARBON_COMPONENT_NAME)//' aerosol', &
                                    density        = MAM_BLACK_CARBON_COMPONENT_DENSITY,           &
                                    hygroscopicity = MAM_BLACK_CARBON_COMPONENT_HYGROSCOPICITY,    &
                                    solubility     = MAM_BLACK_CARBON_COMPONENT_SOLUBILITY, __RC__)
       
       ! SOA
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding component - SOA'
       n = n + 1
       aero_component => self%aerosol_component(n)
       call MAM_AerosolComponentSet(aero_component,                                            &
                                    name           = MAM_SOA_COMPONENT_NAME,                   &
                                    long_name      = trim(MAM_SOA_COMPONENT_NAME)//' aerosol', &
                                    density        = MAM_SOA_COMPONENT_DENSITY,                &
                                    hygroscopicity = MAM_SOA_COMPONENT_HYGROSCOPICITY,         &
                                    solubility     = MAM_SOA_COMPONENT_SOLUBILITY, __RC__)
       
       ! POM
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding component - POM'
       n = n + 1
       aero_component => self%aerosol_component(n)
       call MAM_AerosolComponentSet(aero_component,                                            &
                                    name           = MAM_POM_COMPONENT_NAME,                   &
                                    long_name      = trim(MAM_POM_COMPONENT_NAME)//' aerosol', &
                                    density        = MAM_POM_COMPONENT_DENSITY,                &
                                    hygroscopicity = MAM_POM_COMPONENT_HYGROSCOPICITY,         &
                                    solubility     = MAM_POM_COMPONENT_SOLUBILITY, __RC__)
       
       ! DUST
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding component - dust'
       n = n + 1
       aero_component => self%aerosol_component(n)
       call MAM_AerosolComponentSet(aero_component,                                             &
                                    name           = MAM_DUST_COMPONENT_NAME,                   &
                                    long_name      = trim(MAM_DUST_COMPONENT_NAME)//' aerosol', &
                                    density        = MAM_DUST_COMPONENT_DENSITY,                &
                                    hygroscopicity = MAM_DUST_COMPONENT_HYGROSCOPICITY,         &
                                    solubility     = MAM_DUST_COMPONENT_SOLUBILITY, __RC__)
       
       ! SEASALT
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding component - seasalt'
       n = n + 1
       aero_component => self%aerosol_component(n)
       call MAM_AerosolComponentSet(aero_component,                                                &
                                    name           = MAM_SEASALT_COMPONENT_NAME,                   &
                                    long_name      = trim(MAM_SEASALT_COMPONENT_NAME)//' aerosol', &
                                    density        = MAM_SEASALT_COMPONENT_DENSITY,                &
                                    hygroscopicity = MAM_SEASALT_COMPONENT_HYGROSCOPICITY,         &
                                    solubility     = MAM_SEASALT_COMPONENT_SOLUBILITY, __RC__)

       ! check if the number of components is right
       self%n_aerosol_components = n
       if ( self%n_aerosol_components /= MAM7_AEROSOL_COMPONENTS ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol components.")
       end if

       if (present(rc)) then
           rc = MAM_SUCCESS
       end if

   end subroutine MAM7_AerosolComponentsInit


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM7_AerosolModesInit --- 
!
! !INTERFACE:

   subroutine MAM7_AerosolModesInit(self, verbose, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_Scheme),  intent(inout) :: self        ! MAM scheme/configuration
       logical, optional, intent(in)    :: verbose     ! verbosity flag
       integer, optional, intent(inout) :: rc          ! return code

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:


! !DESCRIPTION: Set up MAM7 machinery
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                 __Iam__('MAM7_AerosolModesInit')
   

       type(MAM_AerosolMode), pointer :: mode

       logical :: verbose_
       integer :: n

 
       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR
       
       if (present(verbose)) then
           verbose_ = verbose
       else
           verbose_ = .false.
       end if

       ! allocate memory for the aerosol modes
       allocate(self%mode(MAM7_MODES), __STAT__)


       n = 0

       ! AITKEN
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding mode - aitken'
       n = n + 1
       mode => self%mode(n)
       call MAM_AerosolModeSet(mode, name               = MAM7_AITKEN_MODE_NAME,               &
                                     long_name          = 'aitken mode',                       &
                                     sigma              = MAM7_AITKEN_MODE_SIGMA,              &
                                     size_default       = MAM7_AITKEN_MODE_SIZE,               &
                                     size_min           = MAM7_AITKEN_MODE_SIZE_MIN,           &
                                     size_max           = MAM7_AITKEN_MODE_SIZE_MAX,           &
                                     rh_deliquescence   = MAM7_AITKEN_MODE_RH_DELIQUESCENCE,   &
                                     rh_crystallization = MAM7_AITKEN_MODE_RH_CRYSTALLIZATION, &
                                     n_species          = MAM7_AITKEN_MODE_SPECIES, __RC__)

       ! ACCUMULATION
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding mode - accumulation'
       n = n + 1
       mode => self%mode(n)
       call MAM_AerosolModeSet(mode, name               = MAM7_ACCUMULATION_MODE_NAME,               &
                                     long_name          = 'accumulation mode',                       &
                                     sigma              = MAM7_ACCUMULATION_MODE_SIGMA,              &
                                     size_default       = MAM7_ACCUMULATION_MODE_SIZE,               &
                                     size_min           = MAM7_ACCUMULATION_MODE_SIZE_MIN,           &
                                     size_max           = MAM7_ACCUMULATION_MODE_SIZE_MAX,           &
                                     rh_deliquescence   = MAM7_ACCUMULATION_MODE_RH_DELIQUESCENCE,   &
                                     rh_crystallization = MAM7_ACCUMULATION_MODE_RH_CRYSTALLIZATION, &
                                     n_species          = MAM7_ACCUMULATION_MODE_SPECIES, __RC__)

       ! PRIMARY CARBON
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding mode - primary carbon'
       n = n + 1
       mode => self%mode(n)
       call MAM_AerosolModeSet(mode, name               = MAM7_PRIMARY_CARBON_MODE_NAME,               &
                                     long_name          = 'primary carbon mode',                       &
                                     sigma              = MAM7_PRIMARY_CARBON_MODE_SIGMA,              &
                                     size_default       = MAM7_PRIMARY_CARBON_MODE_SIZE,               &
                                     size_min           = MAM7_PRIMARY_CARBON_MODE_SIZE_MIN,           &
                                     size_max           = MAM7_PRIMARY_CARBON_MODE_SIZE_MAX,           &
                                     rh_deliquescence   = MAM7_PRIMARY_CARBON_MODE_RH_DELIQUESCENCE,   &
                                     rh_crystallization = MAM7_PRIMARY_CARBON_MODE_RH_CRYSTALLIZATION, &
                                     n_species          = MAM7_PRIMARY_CARBON_MODE_SPECIES, __RC__)

       ! FINE SEASALT
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding mode - fine seasalt'
       n = n + 1
       mode => self%mode(n)
       call MAM_AerosolModeSet(mode, name               = MAM7_FINE_SEASALT_MODE_NAME,               &
                                     long_name          = 'fine seasalt mode',                       &
                                     sigma              = MAM7_FINE_SEASALT_MODE_SIGMA,              &
                                     size_default       = MAM7_FINE_SEASALT_MODE_SIZE,               &
                                     size_min           = MAM7_FINE_SEASALT_MODE_SIZE_MIN,           &
                                     size_max           = MAM7_FINE_SEASALT_MODE_SIZE_MAX,           &
                                     rh_deliquescence   = MAM7_FINE_SEASALT_MODE_RH_DELIQUESCENCE,   &
                                     rh_crystallization = MAM7_FINE_SEASALT_MODE_RH_CRYSTALLIZATION, &
                                     n_species          = MAM7_FINE_SEASALT_MODE_SPECIES, __RC__)

       ! FINE DUST
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding mode - fine dust'
       n = n + 1
       mode => self%mode(n)
       call MAM_AerosolModeSet(mode, name               = MAM7_FINE_DUST_MODE_NAME,               &
                                     long_name          = 'fine dust mode',                       &
                                     sigma              = MAM7_FINE_DUST_MODE_SIGMA,              &
                                     size_default       = MAM7_FINE_DUST_MODE_SIZE,               &
                                     size_min           = MAM7_FINE_DUST_MODE_SIZE_MIN,           &
                                     size_max           = MAM7_FINE_DUST_MODE_SIZE_MAX,           &
                                     rh_deliquescence   = MAM7_FINE_DUST_MODE_RH_DELIQUESCENCE,   &
                                     rh_crystallization = MAM7_FINE_DUST_MODE_RH_CRYSTALLIZATION, &
                                     n_species          = MAM7_FINE_DUST_MODE_SPECIES, __RC__)

       ! COARSE SEASALT
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding mode - coarse seasalt'
       n = n + 1
       mode => self%mode(n)
       call MAM_AerosolModeSet(mode, name               = MAM7_COARSE_SEASALT_MODE_NAME,               &
                                     long_name          = 'coarse seasalt mode',                       &
                                     sigma              = MAM7_COARSE_SEASALT_MODE_SIGMA,              &
                                     size_default       = MAM7_COARSE_SEASALT_MODE_SIZE,               &
                                     size_min           = MAM7_COARSE_SEASALT_MODE_SIZE_MIN,           &
                                     size_max           = MAM7_COARSE_SEASALT_MODE_SIZE_MAX,           &
                                     rh_deliquescence   = MAM7_COARSE_SEASALT_MODE_RH_DELIQUESCENCE,   &
                                     rh_crystallization = MAM7_COARSE_SEASALT_MODE_RH_CRYSTALLIZATION, &
                                     n_species          = MAM7_COARSE_SEASALT_MODE_SPECIES, __RC__)

       ! COARSE DUST
       if (verbose_ .and. MAPL_AM_I_ROOT()) write (*,*) 'adding mode - coarse dust'
       n = n + 1
       mode => self%mode(n)
       call MAM_AerosolModeSet(mode, name               = MAM7_COARSE_DUST_MODE_NAME,               &
                                     long_name          = 'coarse dust mode',                       &
                                     sigma              = MAM7_COARSE_DUST_MODE_SIGMA,              &
                                     size_default       = MAM7_COARSE_DUST_MODE_SIZE,               &
                                     size_min           = MAM7_COARSE_DUST_MODE_SIZE_MIN,           &
                                     size_max           = MAM7_COARSE_DUST_MODE_SIZE_MAX,           &
                                     rh_deliquescence   = MAM7_COARSE_DUST_MODE_RH_DELIQUESCENCE,   &
                                     rh_crystallization = MAM7_COARSE_DUST_MODE_RH_CRYSTALLIZATION, &
                                     n_species          = MAM7_COARSE_DUST_MODE_SPECIES, __RC__)

       ! check if the number of modes is right
       self%n_modes = n
       if ( self%n_modes /= MAM7_MODES ) then
           __raise__ (MAM_GENERAL_ERROR, "Incorrect number of aerosol modes.")
       end if

       if (present(rc)) then
           rc = STATUS
       end if

   end subroutine MAM7_AerosolModesInit


  
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM3_Init --- 
!
! !INTERFACE:

   subroutine MAM3_Init(scheme, verbose, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       type(MAM_Scheme), intent(inout)  :: scheme      ! MAM scheme/configuration 

! !INPUT PARAMETERS:
        logical, optional, intent(in)   :: verbose     ! verbosity flag

! !OUTPUT PARAMETERS:
       integer, optional, intent(inout) :: rc          ! return code

! !DESCRIPTION: Set up MAM3 machinery
!
! !REVISION HISTORY:
!
!  18Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                     __Iam__('MAM3_Init')
  
       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       __raise__ (MAM_NOT_IMPLEMENTED_ERROR, 'MAM3 initialization is not implemented yet.')
  
   end subroutine MAM3_Init

   
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_SchemeGetAerosolComponentIndexFromName --- 
!
! !INTERFACE:

   function MAM_SchemeGetAerosolComponentIndexFromName(self, name, rc) result(ix)

! !USES:

       implicit NONE

       integer :: ix

! !INPUT/OUTPUT PARAMETERS:
       integer, optional, intent(inout) :: rc          ! return code

! !INPUT PARAMETERS:
       type(MAM_Scheme), intent(in)     :: self        ! MAM scheme/configuration
       character(len=*), intent(in)     :: name        ! name of the component

! !OUTPUT PARAMETERS:

! !DESCRIPTION: Returns the index of aerosol component
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
             __Iam__('MAM_SchemeGetAerosolComponentIndexFromName')
  
       integer :: i

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       ix  = 0

       do i = 1, self%n_aerosol_components
           if (self%aerosol_component(i)%name == trim(name)) then
               ix  = i
               exit
           end if
       end do

       if (ix == 0) then
           STATUS = MAM_UNKNOWN_AEROSOL_COMPONENT_ERROR
           VERIFY_(STATUS)
       else
           STATUS = MAM_SUCCESS
           VERIFY_(STATUS)
       end if

       if (present(rc)) then
           rc = STATUS
       end if
   
   end function MAM_SchemeGetAerosolComponentIndexFromName


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_SchemeGetmodeIndexFromName --- 
!
! !INTERFACE:

   function MAM_SchemeGetModeIndexFromName(self, name, rc) result(ix)

! !USES:

       implicit NONE

       integer :: ix

! !INPUT/OUTPUT PARAMETERS:
       integer, optional, intent(inout) :: rc          ! return code

! !INPUT PARAMETERS:
       type(MAM_Scheme), intent(in)     :: self        ! MAM scheme/configuration
       character(len=*), intent(in)     :: name        ! name of the component

! !OUTPUT PARAMETERS:

! !DESCRIPTION: Returns the index of aerosol mode
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                 __Iam__('MAM_SchemeGetModeIndexFromName')
  
       integer :: i

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       ix  = 0

       do i = 1, self%n_modes
           if (self%mode(i)%name == trim(name)) then
               ix = i
               exit
           end if
       end do

       if (ix == 0) then
           STATUS = MAM_GENERAL_ERROR
           VERIFY_(STATUS)
       else
           STATUS = MAM_SUCCESS
           VERIFY_(STATUS)
       end if

       if (present(rc)) then
           rc = STATUS
       end if
   
   end function MAM_SchemeGetModeIndexFromName


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_SchemePrint --- 
!
! !INTERFACE:

   subroutine MAM_SchemePrint(self, rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:
       type(MAM_Scheme), intent(in)   :: self      ! MAM scheme/configuration

! !OUTPUT PARAMETERS:
       integer, optional, intent(out) :: rc        ! return code

! !DESCRIPTION: Returns the index of aerosol mode
!
! !REVISION HISTORY:
!
!  21Nov2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                 __Iam__('MAM_SchemeGetModeIndexFromName')
  
       integer :: c
       integer :: m

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR 
       
       do c = 1, self%n_aerosol_components
           call MAM_AerosolComponentPrint(self%aerosol_component(c), __RC__)
       end do

       do m = 1, self%n_modes
           call MAM_AerosolModePrint(self%mode(m), __RC__)
       end do

       if (present(rc)) rc = STATUS

   end subroutine MAM_SchemePrint


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MAM_SchemeGet --- 
!
! !INTERFACE:

   subroutine MAM_SchemeGet(self, name,                 &
                                  long_name,            &
                                  n_aerosol_components, &
                                  aerosol_components,   &
                                  n_modes,              &
                                  modes,                &
                                  rc)

! !USES:

       implicit NONE

! !INPUT/OUTPUT PARAMETERS:
       character(len=*), optional, intent(inout) :: name
       character(len=*), optional, intent(inout) :: long_name

       integer,          optional, intent(inout) :: n_aerosol_components
       integer,          optional, intent(inout) :: n_modes

       type(MAM_AerosolComponent), optional, &
            pointer, dimension(:), intent(inout) :: aerosol_components
       type(MAM_AerosolMode),      optional, &
            pointer, dimension(:), intent(inout) :: modes

! !INPUT PARAMETERS:
       type(MAM_Scheme),           intent(in)    :: self    ! MAM scheme/configuration

! !OUTPUT PARAMETERS:
       integer,          optional, intent(out)   :: rc      ! return code 
    

! !DESCRIPTION: Quires a MAM_Scheme instance.
!
! !REVISION HISTORY:
!
!  5Dec2011  A. Darmenov  First crack.
!
!EOP
!-------------------------------------------------------------------------
                 __Iam__('MAM_SchemeGet')
  

       if (present(rc)) rc = MAM_NOT_IMPLEMENTED_ERROR

       if (present(name)) then
           name = trim(self%name)
       end if

       if (present(long_name)) then
           long_name = trim(self%long_name)
       end if

       if (present(n_aerosol_components)) then
           n_aerosol_components = self%n_aerosol_components
       end if

       if (present(n_modes)) then
           n_modes = self%n_modes
       end if

       if (present(aerosol_components)) then
           aerosol_components => self%aerosol_component
       end if

       if (present(modes)) then
           modes => self%mode
       end if

       if (present(rc)) rc = MAM_SUCCESS

   end subroutine MAM_SchemeGet


   end module MAM_BaseMod


