!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiStringManipulation_mod
!
! !INTERFACE:
!
      module GmiStringManipulation_mod
!
      implicit none
!
#     include "GmiParameters.h"
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: stringUpperCase, stringLowerCase
      public  :: upperCase, lowerCase
      public  :: constructListNames
!
! !DESCRIPTION:
!  Routines/functions to manipulate strings.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  16April2007 - Initial code.
!
!EOP
!-------------------------------------------------------------------------

   CONTAINS

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: constructListNames
!
! !INTERFACE:

     subroutine constructListNames(speciesNames, names)
!
     implicit none
!
! OUTPUT PARAMETERS:
     character (len=*), intent(  out) :: speciesNames(*)
!
! !INPUT/OUTPUT PARAMETERS:
     character (len=*), intent(inout) :: names
!
! !DESCRIPTION:
!  This routine takes the long string "names" (containing a list
!  of species, separated by commas) and construct a new list
!  where each species name is a string.
!
! !LOCAL VARIABLES:
     integer             :: loc, namesLen, i
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

     i = 1
     loc = index(names, ',')
     if (loc == 0) then
        speciesNames(i) = names
        speciesNames(i) = adjustl(speciesNames(i))
        speciesNames(i) = trim   (speciesNames(i))
     else
        do while(loc > 0)
           speciesNames(i) = names(1:loc-1)
           speciesNames(i) = adjustl(speciesNames(i))
           speciesNames(i) = trim   (speciesNames(i))
           namesLen = len(names)
           names = names(loc+1:namesLen)
           i = i + 1
           loc = index(names, ',')
        end do
        speciesNames(i) = names
        speciesNames(i) = adjustl(speciesNames(i))
        speciesNames(i) = trim   (speciesNames(i))
     end if

     return

     end subroutine constructListNames

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: upperCase
! 
! !INTERFACE:
!
      subroutine upperCase(instr, outstr)
!
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in)  :: instr
!
! !OUTPUT PARAMETERS:
    character(len=*), intent(out) :: outstr
!
! !DESCRIPTION:
!  This routine turns a string to upper case.
!  instr and outstr can be the same variable.
!
! !LOCAL VARIABLES:
    integer :: i, j, ascii, ll, la, ua
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
    la = iachar('a')
    ua = iachar('A')

    ll = len_trim(outstr)
    do i = 1, min(len_trim(instr),ll)
       ascii = iachar( instr(i:i) )
       if (ascii >= la .and. ascii < la+26) then ! in [a,z]
          outstr(i:i) = achar( ascii - la + ua )
       else
          outstr(i:i) = instr(i:i)
       endif
    enddo
    do j = i, ll ! pad with blanks
       outstr(j:j) = ' '
    enddo
    return
  end subroutine upperCase

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lowerCase
!
! !INTERFACE:

  subroutine lowerCase(instr, outstr)
!
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in)  :: instr
!
! !INPUT PARAMETERS:
    character(len=*), intent(out) :: outstr
!
! !DESCRIPTION:
!  Turns a string to lower case.
!  instr and outstr can be the same variable.
!
! !LOCAL VARIABLES:
    integer :: i, j, ascii, ll, la, ua
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!      
!EOP   
!-------------------------------------------------------------------------
!BOC

    la = iachar('a')
    ua = iachar('A')

    ll = len_trim(outstr)
    do i = 1, min(len_trim(instr),ll)
       ascii = iachar( instr(i:i) )
       if (ascii >= ua .and. ascii < ua+26) then ! in [A,Z]
          outstr(i:i) = achar( ascii - ua + la )
       else
          outstr(i:i) = instr(i:i)
       endif
    enddo
    do j = i, ll ! pad with blanks
       outstr(j:j) = ' '
    enddo
    return
  end subroutine lowerCase

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: stringUpperCase
!
! !INTERFACE:

  function stringUpperCase(instr) result(outstr)
!
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in)  :: instr
!
! !RETURN VALUE:
    character(len=MAX_LENGTH_LABELS)    :: outstr
!
! !DESCRIPTION:
!  Turns a character string to upper case.
!
! !LOCAL VARIABLES:
    integer :: i, ascii, la, ua
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!      
!EOP   
!-------------------------------------------------------------------------
!BOC

    la = iachar('a')
    ua = iachar('A')

    outstr = instr
    do i = 1, min(len_trim(instr),len_trim(outstr))
       ascii = iachar( instr(i:i) )
       if (ascii >= la .and. ascii < la+26) then ! in [a,z]
          outstr(i:i) = achar( ascii - la + ua )
       endif
    enddo

    return
  end function stringUpperCase

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: stringLowerCase
!
! !INTERFACE:

  function stringLowerCase(instr) result(outstr)
!
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in)  :: instr
!
! !RETURN VALUE:
    character(len=MAX_LENGTH_LABELS)    :: outstr
!
! !DESCRIPTION:
!  Turns a character string to lower case.
!
! !LOCAL VARIABLES:
    integer :: i, ascii, la, ua
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!      
!EOP   
!-------------------------------------------------------------------------
!BOC
    la = iachar('a')
    ua = iachar('A')

    outstr = instr
    do i = 1, min(len_trim(instr),len_trim(outstr))
       ascii = iachar( instr(i:i) )
       if (ascii >= ua .and. ascii < ua+26) then ! in [A,Z]
          outstr(i:i) = achar( ascii - ua + la )
       endif
    enddo

    return
  end function stringLowerCase

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: string\_i
!
! !INTERFACE:

  function string_i(arg, format) result(str)
!
      implicit none
!
! !INPUT PARAMETERS:
    integer :: arg
    character(len=*),   optional :: format
!
! !RETURN VALUE:
    character(len=48)           :: str
!
! !DESCRIPTION:
!  Writes out an integer variable at the prescribed format.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!      
!EOP   
!-------------------------------------------------------------------------
!BOC
    if (present(format)) then
       write(str,format) arg
    else
       write(str,*) arg
    endif
    return
  end function string_i

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: string\_s
!
! !INTERFACE:

  function string_s(arg, format) result(str)
!
      implicit none
!
! !INPUT PARAMETERS:
    real*4 :: arg
    character(len=*),   optional :: format
!
! !RETURN VALUE:
    character(len=48)           :: str
!
! !DESCRIPTION:
!  Writes out a single precision real variable
!  at the prescribed format.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!      
!EOP   
!-------------------------------------------------------------------------
!BOC
    if (present(format)) then
       write(str,format) arg
    else
       write(str,*) arg
    endif
    return
  end function string_s

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: string\_d
!
! !INTERFACE:

  function string_d(arg, format) result(str)
!
      implicit none
!
! !INPUT PARAMETERS:
    real*8          , intent(in) :: arg
    character(len=*),   optional :: format
!
! !RETURN VALUE:
    character(len=48)           :: str
!
! !DESCRIPTION:
!  Writes out a double precision real variable
!  at the prescribed format.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!      
!EOP   
!-------------------------------------------------------------------------
!BOC

    if (present(format)) then
       write(str,format) arg
    else
       write(str,*) arg
    endif
    return
  end function string_d

!EOC
!-------------------------------------------------------------------------

      end module GmiStringManipulation_mod
