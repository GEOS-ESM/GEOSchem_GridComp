
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_sad_constants.h
!
! DESCRIPTION
!   This include file contains the constants that map a particular type of
!   aerosol SAD (i.e., surface area density) to a numeric value for array
!   indexing.
!
!=============================================================================


      integer, parameter ::  &
     &  ILBSSAD  = 1,    & ! index for liquid binary sulfate       SADs
     &  ISTSSAD  = 2,    & ! index for supercooled ternary sulfate SADs
     &  INATSAD  = 3,    & ! index for "NAT" SADs
     &  IICESAD  = 4,    & ! index for ice   SADs
     &  ISOOTSAD = 5   ! index for soot  SADs

