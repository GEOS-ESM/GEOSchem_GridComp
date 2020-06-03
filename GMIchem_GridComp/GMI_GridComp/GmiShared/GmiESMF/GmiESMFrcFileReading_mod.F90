#include "MAPL_ErrLog.h"
!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiESMFrcFileReading_mod
!
! !INTERFACE:
!
      module GmiESMFrcFileReading_mod
!
! !USES:
      use ESMF
      use MAPL
!
      implicit none
!
      private
!
! !PUBLIC MEMBER FUNCTIONS:
      public  :: rcEsmfReadTable
!     public  :: rcEsmfReadLogical
      public  :: reconstructPhrase

      interface rcEsmfReadTable
        module procedure rcEsmfReadTable2ListInt
        module procedure rcEsmfReadTable2ListReal
        module procedure rcEsmfReadTable2String
        module procedure rcEsmfReadTable2ListWords
      end interface

#     include "GmiParameters.h"
!
! !DESCRIPTION:
! Basic routines for reading resource file using ESMF calls.
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcEsmfReadTable2ListInt
!
! !INTERFACE:
!
      subroutine rcEsmfReadTable2ListInt(config, value, label, rc)
!
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in) :: label
!
! !OUTPUT PARAMETERS:
    integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
    integer          , intent(inOut) :: value(1:)
    type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a table of integers from a resource file.
!
! !LOCAL VARIABLES:
      logical :: endTable
      integer :: STATUS, counter
      integer :: temp
      character(len=ESMF_MAXSTR), parameter :: IAm = "rcEsmfReadTable2ListInt"
      logical :: isPresent
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_ConfigFindLabel(config, label=label, isPresent=isPresent, rc=STATUS )
      VERIFY_(STATUS)

      counter = 0

      if (isPresent) then
         call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
         VERIFY_(STATUS)

         do while (.not. endTable)
            counter = counter + 1

            call ESMF_ConfigGetAttribute(config, temp, rc=STATUS )
            VERIFY_(STATUS)

            value(counter) = temp

            call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
            VERIFY_(STATUS)
         end do
      end if

      if (present(rc)) rc = STATUS

      return

      end subroutine rcEsmfReadTable2ListInt
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcEsmfReadTable2ListReal
!
! !INTERFACE:
!
      subroutine rcEsmfReadTable2ListReal(config, value, label, rc)
!
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in) :: label
!
! !OUTPUT PARAMETERS:
    integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
    real*8, intent(inOut) :: value(1:)
    type(ESMF_Config), intent(inOut) :: config
!
! !DESCRIPTION:
! Reads in a table of real numbers from a resource file.
!
! !LOCAL VARIABLES:
      logical :: endTable
      integer :: STATUS, counter
      real(ESMF_KIND_R8) :: temp
      character(len=ESMF_MAXSTR), parameter :: IAm = "rcEsmfReadTable2ListReal"
      logical :: isPresent
!EOP
!------------------------------------------------------------------------------
!BOC
      call ESMF_ConfigFindLabel(config, label=label, isPresent=isPresent, rc=STATUS )
      VERIFY_(STATUS)

      counter = 0

      if (isPresent) then
         call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
         VERIFY_(STATUS)

         do while (.not. endTable)
            counter = counter + 1

            call ESMF_ConfigGetAttribute(config, temp, rc=STATUS )
            VERIFY_(STATUS)

            value(counter) = temp

            call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
            VERIFY_(STATUS)
         end do
      end if

      if (present(rc)) rc = STATUS

      return

      end subroutine rcEsmfReadTable2ListReal
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcEsmfReadTable2ListWords
!
! !INTERFACE:
!
      subroutine rcEsmfReadTable2ListWords(config, value, label, rc)
! 
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in) :: label
! 
! !OUTPUT PARAMETERS:
    integer, optional, intent(out) :: rc 
! 
! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(inOut) :: value(1:)
    type(ESMF_Config), intent(inOut) :: config
!     
! !DESCRIPTION:
! Read in a table of words from a resource file and construct a list of word.
!
! !LOCAL VARIABLES:
      integer :: STATUS, counter
      logical :: endTable
      character(len=1) :: cValue
      character(len=ESMF_MAXSTR) :: temp
      character(len=ESMF_MAXSTR), parameter :: IAm = "rcEsmfReadTable2ListWords"
      logical :: isPresent
!EOP
!------------------------------------------------------------------------------
!BOC     
      call ESMF_ConfigFindLabel(config, label=label, isPresent=isPresent, rc=STATUS )
      VERIFY_(STATUS)

      counter = 0

      if (isPresent) then
         call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
         VERIFY_(STATUS)

         do while (.not. endTable)
            counter = counter + 1

            call ESMF_ConfigGetAttribute(config, temp, rc=STATUS )
            VERIFY_(STATUS)

            call reconstructPhrase(temp) 

            value(counter) = temp

            call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
            VERIFY_(STATUS)
         end do
      end if

      if (present(rc)) rc = STATUS
            
      return

      end subroutine rcEsmfReadTable2ListWords
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcEsmfReadTable2String
!
! !INTERFACE:
!
      subroutine rcEsmfReadTable2String(config, value, label, rc)
! 
      implicit none
!
! !INPUT PARAMETERS:
    character(len=*), intent(in) :: label
! 
! !OUTPUT PARAMETERS:
    character(len=MAX_STRING_LENGTH), intent(out) :: value
    integer, optional, intent(out) :: rc 
! 
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Config), intent(inOut) :: config
!     
! !DESCRIPTION:
! Reads in a table of words from a resource file.
!
! !LOCAL VARIABLES:
      integer :: STATUS
      logical :: firstIter, endTable
      character(len=1) :: cValue
      character(len=ESMF_MAXSTR) :: tempWord
      character(len=ESMF_MAXSTR), parameter :: IAm = "rcEsmfReadTable2String"
      logical :: isPresent
!EOP
!------------------------------------------------------------------------------
!BOC     
      firstIter = .true.

      call ESMF_ConfigFindLabel(config, label=label, isPresent=isPresent, rc=STATUS )
      VERIFY_(STATUS)
            
      value = ''

      if (isPresent) then
         call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
         VERIFY_(STATUS)
            
         do while (.not. endTable)
            call ESMF_ConfigGetAttribute(config, tempWord, rc=STATUS )
            VERIFY_(STATUS)

            call reconstructPhrase(tempWord)
         
            if (firstIter) then
               firstIter = .false.
               value = trim(tempWord)
            else
               value = trim(value)//', '//trim(tempWord)
            end if

            call ESMF_ConfigNextLine  (config, tableEnd=endTable, rc=STATUS )
            VERIFY_(STATUS)
         end do
      end if

      if (present(rc)) rc = STATUS

      return

      end subroutine rcEsmfReadTable2String
!EOC
!------------------------------------------------------------------------------
!!BOP
!
!! !IROUTINE: rcEsmfReadLogical     INSTEAD USE: ESMF_ConfigGetAttribute
!!
!! !INTERFACE:
!!
!      subroutine rcEsmfReadLogical(config, value, label, default, rc)
!!
!      implicit none
!!
!! !INPUT PARAMETERS:
!    character(len=*) , intent(in) :: label
!    logical, optional, intent(in) :: default
!!
!! !OUTPUT PARAMETERS:
!    logical, intent(out) :: value
!    integer, optional, intent(out) :: rc
!!
!! !INPUT/OUTPUT PARAMETERS:
!    type(ESMF_Config), intent(inOut) :: config
!!
!! !DESCRIPTION:
!! Reads in a logical variable from a resource file.
!! Note that ESMF does not have a routine to read logical variables.
!! We assume that the variable in the resource file is character and we
!! do the conversion after the reading.
!!
!! !LOCAL VARIABLES:
!      integer :: STATUS
!      character(len=1) :: cDefault, cValue
!      character(len=ESMF_MAXSTR), parameter :: IAm = "rcEsmfReadLogical"
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!      cDefault = 'F'
!      if (present(default)) then
!         if (default) cDefault = 'T'
!      end if
!
!      call ESMF_ConfigGetChar(config, cValue, label=label, default=cDefault, rc=STATUS )
!      VERIFY_(STATUS)
!
!      if (present(rc)) rc = STATUS
!
!      if ((cValue == 'T') .or. (cValue == 't')) then
!         value = .true.
!      else
!         value = .false.
!      end if
!      
!      return
!      
!      end subroutine rcEsmfReadLogical
!!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: reconstructPhrase
!
! !INTERFACE:

     subroutine reconstructPhrase(phrase)
!
     implicit none
!
! !INPUT/OUTPUT PARAMETERS:
     character (len=*), intent(inOut) :: phrase
!
! !DESCRIPTION:
!  Takes a long string and reconstruct the phrase/name it represents.
!
! !LOCAL VARIABLES:
     integer                     :: loc, namesLen
     logical                     :: first
     character (len=ESMF_MAXSTR) :: names
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
     first = .true.
     names = phrase

     loc = index(names, '*')
     if (loc == 0) then
        phrase = names
        phrase = adjustl(phrase)
        phrase = trim   (phrase)
     else
        do while(loc > 0)
           if (first) then
              first  = .false.
              phrase = names(1:loc-1)
           else
              phrase = trim(phrase)//' '//names(1:loc-1)
           end if

           phrase = adjustl(phrase)
           phrase = trim   (phrase)
           namesLen = len(names)
           names = names(loc+2:namesLen)
           loc = index(names, '*')
        end do
        phrase = trim(phrase)//' '//names
        phrase = adjustl(phrase)
        phrase = trim   (phrase)
     end if

     return

     end subroutine reconstructPhrase

!EOC
!-------------------------------------------------------------------------
      end module GmiESMFrcFileReading_mod
