#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!      NASA/GSFC, Atmospheric Chemistry and Dynamics Lab, Code 614       !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Runtime_RegistryMod --- Chemistry Registry Class
!
! !INTERFACE:
!

   module  Runtime_RegistryMod

! !USES:

   USE ESMF
   USE MAPL

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  Runtime_Registry  ! A listing of species, units and long names
                           
!
! !PUBLIIC MEMBER FUNCTIONS:
!
   PUBLIC  Runtime_RegistryCreate   ! Constructor from RC file
   PUBLIC  Runtime_RegistryDestroy  ! Destructor
   PUBLIC  Runtime_RegistryPrint    ! Prints a summary of the registry

!
! !DESCRIPTION:
!
!  This module implements a simple registry.
!
!
! !REVISION HISTORY:
!
!  2022.08.26  Manyin  First crack
!
!EOP
!-------------------------------------------------------------------------

  integer, parameter :: REGISTER_NAME_LENGTH     = 32
  integer, parameter :: REGISTER_UNITS_LENGTH    = 32
  integer, parameter :: REGISTER_LONGNAME_LENGTH = 64

! integer, parameter :: RC_DATA_LINE    = 0
! integer, parameter :: RC_END_OF_TABLE = 1
! integer, parameter :: RC_END_OF_FILE  = 2

! Registry
! --------
  type Runtime_Registry

     integer :: nq    ! Total number of tracers 

     character(len=REGISTER_NAME_LENGTH),     pointer :: vname(:)   ! (nq), variable short name
     character(len=REGISTER_UNITS_LENGTH),    pointer :: vunits(:)  ! (nq), variable units
     character(len=REGISTER_LONGNAME_LENGTH), pointer :: vtitle(:)  ! (nq), variable long  name

  end type Runtime_Registry

CONTAINS

!------------------------------------------------------------------------
!     NASA/GSFC, Atmospheric Chemistry and Dynamics Lab, Code 614       !
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Runtime_RegistryCreate --- Construct Chemistry Registry
!
! !INTERFACE:
!

  Function Runtime_RegistryCreate ( rcfile, table_name, rc )

  implicit none
  type(Runtime_Registry) Runtime_RegistryCreate 

! !USES:

! !INPUT PARAMETERS:

   character(len=*), intent(in) :: rcfile      ! Resource file name
   character(len=*), intent(in) :: table_name  ! includes '::'

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc            ! Error return code:
                                          !  0 - all is well

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
!  2022.08.26  Manyin  First crack
!
!EOP
!-------------------------------------------------------------------------

  ! Possible return codes from get_line:
  integer, parameter :: RC_DATA_LINE    = 1
  integer, parameter :: RC_END_OF_TABLE = 2
  integer, parameter :: RC_END_OF_FILE  = 3

  integer, parameter :: TOKEN_LENGTH  = 256

    __Iam__('Runtime_RegistryCreate')   ! NOTE: this macro declares STATUS
                                        ! ALSO: Never set Iam = TRIM(Iam) // suffix
                                        !       because Iam is a SAVED varaible

   character(len=*), parameter ::  myname = 'Runtime_RegistryCreate'

   character(len=TOKEN_LENGTH) ::  str_arr(3)

   type(ESMF_Config)  ::  cf

   type(Runtime_Registry) :: this
   integer :: nq, item_count, retcode
   integer :: i

   rc = 0
   str_arr(:) = ""
                
   cf = ESMF_ConfigCreate(__RC__)

   call ESMF_ConfigLoadFile(cf, rcfile, rc=rc)
   _ASSERT(rc==0, TRIM(Iam)//': Cannot load RC file '//TRIM(rcfile))

   call ESMF_ConfigGetAttribute(cf, nq, label='tracer_count:', rc=rc)
   _ASSERT(rc==0, TRIM(Iam)//': Cannot find tracer_count in file '//TRIM(rcfile))

!  Allocate memory in registry
!  ---------------------------
   this%nq = nq
   allocate ( this%vname(nq), this%vunits(nq), this%vtitle(nq), __STAT__ )

   call ESMF_ConfigFindLabel(cf, table_name, rc=rc)
   _ASSERT(rc==0, TRIM(Iam)//': Cannot find '//TRIM(table_name)//' in file '//TRIM(rcfile))

   do i=1,nq
      call get_line ( cf, 3, str_arr, item_count, retcode )
      select case( retcode )
        case( RC_END_OF_FILE  )
          _ASSERT(.FALSE., TRIM(Iam)//': early EOF in file '//TRIM(rcfile))
        case( RC_END_OF_TABLE )
          _ASSERT(.FALSE., TRIM(Iam)//': table too short '//TRIM(table_name)//' in file '//TRIM(rcfile))
        case( RC_DATA_LINE    )
          _ASSERT(item_count==3, TRIM(Iam)//': fewer than 3 entries in '//TRIM(table_name)//' in file '//TRIM(rcfile))
          this%vname(i)  = str_arr(1)
          this%vunits(i) = str_arr(2)
          this%vtitle(i) = str_arr(3)
      end select
   end do

   call ESMF_ConfigDestroy(cf, __RC__)

!  All done
!  --------
   Runtime_RegistryCreate = this
   
   return 

!                 -----------------
!                 Internal Routines
!                 -----------------

   CONTAINS

!     -------------------
!     GET_LINE
!     Advance one line and then try to read <expected_entries> items
!     Retcode will be set to one of these values:
!       RC_END_OF_FILE    - cannot advance a line
!       RC_END_OF_TABLE   - first item is '::'
!       RC_DATA_LINE      - at least one entry has been put into str_arr
!     Note that ESMF automatically skips over blank lines and comment lines
!     -------------------
      subroutine get_line ( cf, expected_entries, str_arr, item_count, retcode )
!     -------------------
      type(ESMF_Config),           intent(inout)  :: cf
      integer,                     intent(in)     :: expected_entries  ! read this many items
      character(len=TOKEN_LENGTH), intent(out)    :: str_arr(*)   ! space for one or more items
      integer,                     intent(out)    :: item_count   ! how many were successfully read
      integer,                     intent(out)    :: retcode      ! see possible values above

      integer :: i

      call ESMF_ConfigNextLine(cf, rc=rc)
      if ( rc/=0 ) then
        retcode = RC_END_OF_FILE
        return
      end if

      ! Because ESMF skips over blank lines, rc should always be 0:
      call ESMF_ConfigGetAttribute(cf, str_arr(1), rc=rc)
      if ( rc/=0 ) then
        retcode = RC_END_OF_FILE
        return
      end if
      if ( INDEX(str_arr(1), '::' ) == 1 ) then
        retcode = RC_END_OF_TABLE
        return
      end if

      retcode = RC_DATA_LINE
      item_count = 1
      do i=2,expected_entries
        call ESMF_ConfigGetAttribute(cf, str_arr(i), rc=rc)
        if (rc==0) item_count = item_count + 1
      end do

      return
      
      end subroutine get_line

 end Function Runtime_RegistryCreate

!------------------------------------------------------------------------
!     NASA/GSFC, Atmospheric Chemistry and Dynamics Lab, Code 614       !
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Runtime_RegistryDestroy --- Destruct Chemisty Registry
!
! !INTERFACE:
!
  subroutine Runtime_RegistryDestroy ( this, rc )

! !USES:

  implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(Runtime_Registry), intent(inout) :: this

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - 

! !DESCRIPTION: Destructor for Registry object.
!
! !REVISION HISTORY:
!
!  2022.08.26  Manyin  First crack
!
!EOP
!-------------------------------------------------------------------------

   __Iam__('Runtime_RegistryDestroy')   ! NOTE: this macro declares STATUS

   rc = 0
   this%nq = -1 
   deallocate ( this%vname, this%vunits, this%vtitle, __STAT__ )

end subroutine Runtime_RegistryDestroy 


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Runtime_RegistryPrint --- Print summary of Chemistry Registry
!
! !INTERFACE:
!
   SUBROUTINE Runtime_RegistryPrint( reg, mod_name )

! !USES:

! !INPUT PARAMETERS:

   IMPLICIT none
   TYPE(Runtime_Registry) :: reg
   character(len=*)       :: mod_name  ! name of the calling module

! !OUTPUT PARAMETERS:


! !DESCRIPTION:
!
!   Prints summary of Chemistry Registry
!
! !REVISION HISTORY:
!
!  2022.08.26  Manyin  First crack
!
!EOP
!-------------------------------------------------------------------------

   PRINT *
   PRINT *,'****'
   PRINT *,'****       Summary of the '//mod_name//' Registry'
   PRINT *,'****            from Runtime_RegistryPrint'
   PRINT *,'****'
   WRITE(*,FMT="(' ','       Number of species: ',I3)") reg%nq

   CALL reg_prt_( mod_name, reg%nq )

   PRINT *

  101 FORMAT(/,'       Number of species: ',I3)

   RETURN
   
   CONTAINS
   
      SUBROUTINE reg_prt_ ( compName, n )
      
      IMPLICIT none
      CHARACTER(LEN=*), INTENT(IN) :: compName
      INTEGER, INTENT(IN) :: n
      INTEGER :: i
      CHARACTER(LEN=7) :: string
      
      string = 'species'
      
      WRITE(*,101) TRIM(compName),n,string
      DO i = 1, n
       WRITE(*,201) i,TRIM(reg%vname(i)),TRIM(reg%vunits(i)),TRIM(reg%vtitle(i))
      END DO

  101 FORMAT(/,' Component ',A,' has ',I3,' ',A7,/, &
	     ' No ',2X,'  Name  ',2X,'   Units  ',2X,'Description',/, &
	     ' ---',2X,'--------',2X,'----------',2X,'-----------')
  201 FORMAT(' ',I3,2X,A8,2X,A10,2X,A)
          
  END SUBROUTINE reg_prt_
  
END SUBROUTINE Runtime_RegistryPrint

 end module Runtime_RegistryMod

