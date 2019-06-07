module GmiCheckNamelistFile_mod

      implicit none

      private
      public  :: CheckNamelistOptionRange
      public  :: CheckNamelistRead
      public  :: PositionNamelistPointer

      contains

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CheckNamelistOptionRange
!
! DESCRIPTION
!   This routine checks to see if an option is within the range given.
!   If not, it prints an error message and exits.
!
! ARGUMENTS
!   option_char_name : option character name
!   option_to_check  : option to check range of
!   loval            : low  value of range
!   hival            : high value of range
!
!-----------------------------------------------------------------------------

      subroutine CheckNamelistOptionRange  &
     &  (option_char_name, option_to_check, loval, hival)

      use GmiPrintError_mod, only : GmiPrintError

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      character (len=*) :: option_char_name
      integer :: option_to_check
      integer :: loval
      integer :: hival

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: clen

!     ----------------
!     Begin Execution.
!     ----------------

      if ((option_to_check < loval) .or. (option_to_check > hival)) then
         clen = Scan (option_char_name, '_opt', back = .true.)
         err_msg =  &
     &            'Bad value in Check_Opt_Range for option:  ' //  &
     &            option_char_name(1:clen) // ' =>'
         call GmiPrintError (err_msg, .true., 1, option_to_check, 0, 0, 0.0d0, 0.0d0)
      end if

      return

      end subroutine CheckNamelistOptionRange

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CheckNamelistRead
!
! DESCRIPTION
!   This routine handles an error that occurs while reading a section of the
!   namelist file.
!
! ARGUMENTS
!   namelist_section : character name of namelist section
!   ios              : I/O status flag
!
!-----------------------------------------------------------------------------

      subroutine CheckNamelistRead (namelist_section, ios)

      use GmiPrintError_mod, only : GmiPrintError

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      character (len=*) :: namelist_section
      integer :: ios

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

!     ----------------
!     Begin Execution.
!     ----------------

      err_msg = 'Namelist error in section:  ' // namelist_section

      call GmiPrintError (err_msg, .true., 1, ios, 0, 0, 0.0d0, 0.0d0)

      return

      end subroutine CheckNamelistRead

!-----------------------------------------------------------------------------
!
! ROUTINE
!   PositionNamelistPointer
!
! DESCRIPTION
!   This routine positions the namelist file pointer at a particular
!   section of the namelist.
!
! ARGUMENTS
!   nllun           : namelist logical unit number
!   nlsect_name_len : length of nlsect_name
!   nlsect_name     : name of the namelist section to set the namelist file
!                     pointer to
!
!-----------------------------------------------------------------------------

      subroutine PositionNamelistPointer  &
     &  (nllun, nlsect_name_len, nlsect_name)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nllun
      integer :: nlsect_name_len
      character (len=*) :: nlsect_name

!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=80) :: aline

!     ----------------
!     Begin execution.
!     ----------------

      Rewind (nllun)

      nl_loop1: do
         Read (nllun, *) aline
         aline = Adjustl (aline)

         if (aline(1:nlsect_name_len) == nlsect_name) then
            Backspace (nllun)
            exit nl_loop1
         end if

      end do nl_loop1

      return

      end subroutine PositionNamelistPointer
     
end module GmiCheckNamelistFile_mod
