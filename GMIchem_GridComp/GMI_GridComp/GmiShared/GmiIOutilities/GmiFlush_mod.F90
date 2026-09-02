
      module GmiFlush_mod

      implicit none

      private
      public  :: GmiFlush

!-----------------------------------------------------------------------------

      CONTAINS

!-----------------------------------------------------------------------------
!
! ROUTINE
!   GmiFlush
!
! DESCRIPTION
!   This routine is a port-a-potty; it flushes anywhere.
!
! ARGUMENTS
!    unit_no : logical unit number to flush
!
!-----------------------------------------------------------------------------

      subroutine GmiFlush (unit_no)

      implicit none

#     include "gem_sys_options.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: unit_no

!     ----------------
!     Begin execution.
!     ----------------

#if (ARCH_OPTION == ARCH_IBM_SP)

      call Flush_ (unit_no)

#else

      call Flush (unit_no)

#endif

      return

      end subroutine GmiFlush

      end module GmiFlush_mod
