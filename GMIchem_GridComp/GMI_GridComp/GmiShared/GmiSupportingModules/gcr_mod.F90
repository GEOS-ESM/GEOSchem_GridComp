      module gcr_mod
!------------------------------------------------------------------
!  Arrays and dimensions needed for calculation of the NOx emission
!   from Galactic Cosmic Rays from the sun. Depends on number of
!   sunspots and a linear fit to the emission.
!------------------------------------------------------------------
!  NOTE: This code depends on GMI routines that call netCDF;
!        The necessary routines are included in this module,
!        but it would be better if they lived in their own module.
!        See GMI source code, under Shared/NcUtils_Single/
!------------------------------------------------------------------

!
! !USES:
   USE ESMF
   USE MAPL

!  For the netCDF routines:
   use GmiPrintError_mod, only : GmiPrintError

!  For the netCDF routines:
#define MAX_LENGTH_ERROR_MSG 128


!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  ::   READ_GCR_FILE
  public  ::   USE_GCR_DATA
  public  ::   SET_GCR_EMISS
  public  ::   GET_GCR_EMISS
  public  ::   INIT_GCR_DIAG
  public  ::   Finalize_GCR

!... variables for calculation of the Galactic Cosmic Ray production
!...  of NOx (read in)
!...  ( Items we read directly from the file: )
  INTEGER,              SAVE :: GCR_LAT_COUNT, GCR_LEV_COUNT, GCR_TIME_COUNT ! Dimensions
  INTEGER, ALLOCATABLE, SAVE :: GCR_YEARS (:)
  REAL*8,  ALLOCATABLE, SAVE :: GCR_LATS (:)
  REAL*8,  ALLOCATABLE, SAVE :: GCR_SUNSPOT_COEFFS (:)
  REAL*8,  ALLOCATABLE, SAVE :: GCR_AINTCP_COEFFS (:,:)
  REAL*8,  ALLOCATABLE, SAVE :: GCR_SLOPE_COEFFS (:,:)
  REAL*8,  ALLOCATABLE, SAVE :: GCR_LOGPLEVS (:)

  REAL*8,  ALLOCATABLE, SAVE :: GCR_EMISS (:,:,:)  ! For GMI to use in diagnostics

!-------------------------------------------------------------------------

  CONTAINS

!-------------------------------------------------------------------------
  subroutine READ_GCR_FILE ( GCR_FILENAME , rc )
    IMPLICIT NONE

#     include "netcdf.inc"

    CHARACTER (LEN=*), INTENT(IN) :: GCR_FILENAME

! !DESCRIPTION:
! Reads a 2D parameters of Galactic Cosmic Ray NOx emissions.
!   Data from Charley Jackman via the Goddard CTM.
!

! !LOCAL VARIABLES:
    integer :: FILE_ID
    integer :: count1d (1)
    integer :: start1d(1)
    integer :: count2d (2)
    integer :: start2d(2)
    integer :: ierr, varid

! !OUTPUT PARAMETERS:
    integer, intent(out) :: rc

    rc = 0 

!   call Ncop_Rd (FILE_ID, GCR_FILENAME)
    ierr = Nf_Open (GCR_FILENAME, NF_NOWRITE, FILE_ID)
    if (ierr /= NF_NOERR) then
      print*,'In Ncop_Rd, cannot open:  ' // Trim (GCR_FILENAME)
      rc = 99
      return 
    end if

!... set dimensions of netCDF file for GCR parameterization
    call Ncget_Dimlen (FILE_ID, 'latitude_dim', GCR_LAT_COUNT  )
    call Ncget_Dimlen (FILE_ID, 'pressure_dim', GCR_LEV_COUNT  )
    call Ncget_Dimlen (FILE_ID,     'time_dim', GCR_TIME_COUNT )

!... read in data
    start1d(:) = (/ 1 /)
    count1d(1) = GCR_TIME_COUNT

!... get years of input data
    allocate( GCR_YEARS(GCR_TIME_COUNT) )
    call Ncrd_1d_Int (GCR_YEARS, FILE_ID, 'time_dim', start1d, count1d)

!... get sunspot coeff of input data
    allocate( GCR_SUNSPOT_COEFFS(GCR_TIME_COUNT) )
    call Ncrd_1d (GCR_SUNSPOT_COEFFS, FILE_ID, 'GCR_Sunspot_coeff', start1d, count1d)

!... get pressures of input data (hPa)
    allocate( GCR_LOGPLEVS(GCR_LEV_COUNT ) )
    count1d(1) = GCR_LEV_COUNT
    call Ncrd_1d (GCR_LOGPLEVS, FILE_ID, 'pressure_dim', start1d, count1d)
!... take log for linear-log(p) interpolations
    GCR_LOGPLEVS = log(GCR_LOGPLEVS)

!... get latitudes of input data
    allocate( GCR_LATS(GCR_LAT_COUNT ) )
    count1d(1) = GCR_LAT_COUNT
    call Ncrd_1d (GCR_LATS, FILE_ID, 'latitude_dim', start1d, count1d)

!... read in the GCR parameters
    start2d(:) = (/ 1, 1 /)
    count2d(:) = (/ GCR_LAT_COUNT , GCR_LEV_COUNT  /)

    allocate( GCR_SLOPE_COEFFS(GCR_LAT_COUNT ,GCR_LEV_COUNT ) )
    call Ncrd_2d (GCR_SLOPE_COEFFS, FILE_ID, 'GCR_Slope_coeff', start2d, count2d)

    allocate( GCR_AINTCP_COEFFS(GCR_LAT_COUNT ,GCR_LEV_COUNT ) )
    call Ncrd_2d (GCR_AINTCP_COEFFS, FILE_ID, 'GCR_Aintcp_coeff', start2d, count2d)

!   call Nccl (FILE_ID)
    if (ierr .EQ. NF_NOERR)  ierr = Nf_Close (FILE_ID)

    if (ierr /= NF_NOERR) then
      print*,'Trouble with netCDF file ' // Trim (GCR_FILENAME)
      rc = 99 
      return 
    end if

    RETURN
  end subroutine READ_GCR_FILE
  
!-------------------------------------------------------------------------
  subroutine USE_GCR_DATA ( year, pressure, latdeg, &
                            sunspot, slope, aintcp, &
                            i1, i2, ju1, j2, k1, k2 )

!...  This routine was not up to the task;
!     use GmiInterpolation_mod     , only : Interp_Bilinear
!...  I wrote a new version which allows A to be outside
!...  the range of vector V.

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: year
      REAL*8,  INTENT(IN)  :: pressure(i1:i2, ju1:j2, k1:k2) ! hPa
      REAL*8 , INTENT(IN)  :: latdeg (ju1:j2)
      REAL*8,  INTENT(OUT) :: sunspot
      REAL*8,  INTENT(OUT) ::    slope(i1:i2, ju1:j2, k1:k2)
      REAL*8,  INTENT(OUT) ::   aintcp(i1:i2, ju1:j2, k1:k2)
      INTEGER, INTENT(IN)  :: i1, i2, ju1, j2, k1, k2


! !DESCRIPTION:
!   We have the values already read in from the file.
!   The caller provides PRESSURE array, and the year
!   The caller also gives us arrays to fill:
!     slope
!     aintcp
!

! !LOCAL VARIABLES:
  integer :: i, j, k
  integer :: idxyr
  real*8  :: logp

!... calculate index for date we want
    idxyr = year-GCR_YEARS(1)+1
!... make sure index is in range
    do while (idxyr < 1)
      idxyr = idxyr+11
    enddo
    do while (idxyr > GCR_TIME_COUNT)
      idxyr = idxyr-11
    enddo

    sunspot = GCR_SUNSPOT_COEFFS(idxyr)

!... interpolate input data to current model grid
  do k = k1, k2
    do j = ju1, j2
      do i = i1, i2

        logp = LOG( PRESSURE(i,j,k) ) ! hPa

!... interpolate gcr_slope to current grid
        call bi_linear_interp(latdeg(j), logp, slope(i,j,k),  &
                 GCR_LATS, GCR_LOGPLEVS, GCR_LAT_COUNT , GCR_LEV_COUNT , GCR_SLOPE_COEFFS)

!... interpolate gcr_aintcp to current grid
        call bi_linear_interp(latdeg(j), logp, aintcp(i,j,k),  &
                 GCR_LATS, GCR_LOGPLEVS, GCR_LAT_COUNT , GCR_LEV_COUNT , GCR_AINTCP_COEFFS)
      enddo
    enddo
  enddo

    RETURN
  end subroutine USE_GCR_DATA
!-------------------------------------------------------------------------
  subroutine   INIT_GCR_DIAG ( i1, i2, ju1, j2, k1, k2 )
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: i1, i2, ju1, j2, k1, k2

    ALLOCATE( GCR_EMISS(i1:i2,ju1:j2,k1:k2) )
    GCR_EMISS = 0.0

    RETURN
  end subroutine   INIT_GCR_DIAG
!-------------------------------------------------------------------------
  subroutine   SET_GCR_EMISS ( incoming, &
                               i1, i2, ju1, j2, k1, k2 )
      IMPLICIT NONE

      REAL*8,  INTENT(IN)  ::  incoming(i1:i2, ju1:j2, k1:k2)  ! should be molec/cm3/s
      INTEGER, INTENT(IN)  :: i1, i2, ju1, j2, k1, k2

    GCR_EMISS = incoming
    RETURN
  end subroutine   SET_GCR_EMISS
!-------------------------------------------------------------------------
  subroutine GET_GCR_EMISS ( outgoing, &
                             i1, i2, ju1, j2, k1, k2 )
      IMPLICIT NONE

      REAL*8,  INTENT(OUT) :: outgoing(i1:i2, ju1:j2, k1:k2)
      INTEGER, INTENT(IN)  :: i1, i2, ju1, j2, k1, k2

    outgoing = GCR_EMISS
    RETURN
  end subroutine GET_GCR_EMISS
  
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d
!
! !INTERFACE:
!
      subroutine Ncrd_1d (varrd_1d, ncid, varname, start1d, count1d)
!
      implicit none
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start1d   : vector specifying the index in varrd_1d where
!!               the first of the data values will be read
!!    count1d    : varrd_1d dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start1d(1)
      integer          , intent(in)   :: count1d (1)
!
! !OUTPUT PARAMETERS:
!!    varrd_1d : array to fill
      real*8           , intent(out)  :: varrd_1d(count1d(1))
!
! !DESCRIPTION:
!  Reads in a 1D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4  :: varrd_1d_tmp(count1d(1))

!
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, varname, varid)
!
      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d #1:  ' // Trim (varname) // &
                   ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if
!
      ierr =  Nf_Get_Vara_Real   (ncid, varid, start1d, count1d, varrd_1d_tmp)
!
      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if
!
      varrd_1d(:) = varrd_1d_tmp(:)
!
      return
!
      end subroutine Ncrd_1d

!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncrd_1d_Int
!
! !INTERFACE:
!
      subroutine Ncrd_1d_Int (varrd_1di, ncid, varname, start1d, count1d)
!
      implicit none
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start1d   : vector specifying the index in varrd_1di where
!!               the first of the data values will be read
!!    count1d    : varrd_1di dimension
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start1d(1)
      integer          , intent(in)   :: count1d (1)
!
! !OUTPUT PARAMETERS:
!!    varrd_1di : intger array to fill
      integer          , intent(out)  :: varrd_1di(count1d(1))
!
! !DESCRIPTION:
!  Reads in a 1D NetCDF integer array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
!
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, varname, varid)
!
      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d_Int #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if
!
!
      ierr = Nf_Get_Vara_Int (ncid, varid, start1d, count1d, varrd_1di)
!
      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_1d_Int #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if
!
      return
!
      end subroutine Ncrd_1d_Int
!EOC
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine Ncrd_2d (varrd_2d, ncid, varname, start2d, count2d)
!
      implicit none
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!    ncid     : NetCDF file id to read array input data from
!!    varname  : NetCDF variable name for array
!!    start2d   : vector specifying the index in varrd_2d where
!!               the first of the data values will be read
!!    count2d    : varrd_2d dimensions
      integer          , intent(in)   :: ncid
      character (len=*), intent(in)   :: varname
      integer          , intent(in)   :: start2d(2)
      integer          , intent(in)   :: count2d (2)
!
! !OUTPUT PARAMETERS:
!!    varrd_2d : array to fill
      real*8           , intent(out)  :: varrd_2d(count2d(1), count2d(2))
!
! !DESCRIPTION:
!  Reads in a 2D NetCDF real array and does some error checking.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: ierr
      integer             :: varid
      real*4              :: varrd_2d_tmp(count2d(1), count2d(2))
!
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOC
      ierr = Nf_Inq_Varid (ncid, varname, varid)
!
      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d #1:  ' // Trim (varname) // &
                  ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if
!
      ierr = Nf_Get_Vara_Real   (ncid, varid, start2d, count2d, varrd_2d_tmp)
!
      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncrd_2d #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, varid, 0, 0.0d0, 0.0d0)
      end if
!
      varrd_2d(:,:) = varrd_2d_tmp(:,:)
!
      return
!
      end subroutine Ncrd_2d
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncget_Dimlen
!
! !INTERFACE:
!
      subroutine Ncget_Dimlen (ncid, dim_name, dim_len)
!
! !USES:
!
      use GmiPrintError_mod, only : GmiPrintError
!
      implicit none
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!!  dim_name : netCDF dimension name
!!  ncid     : netCDF file id
      character (len=*), intent(in) :: dim_name
      integer,           intent(in) :: ncid
!
! !OUTPUT PARAMETERS:
!!  dim_len: NetCDF dimension length
      integer,           intent(out)   :: dim_len
!
! !DESCRIPTION:
!  Return the length of a given NetCDF dimension.
!
! !LOCAL VARIABLES:
      character (len=MAX_LENGTH_ERROR_MSG) :: err_msg
      integer             :: dimid
      integer             :: ierr
!
! !AUTHOR: 
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      ierr = Nf_Inq_Dimid  (ncid, dim_name, dimid)

      if (ierr /= NF_NOERR) then 
        err_msg = 'In Ncget_Dimlen #1:  ' // Trim (dim_name) // &
                   ', ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
      end if

      ierr = Nf_Inq_Dimlen (ncid, dimid, dim_len)

      if (ierr /= NF_NOERR) then
        err_msg = 'In Ncget_Dimlen #2:  ' // Nf_Strerror (ierr)
        call GmiPrintError (err_msg, .true., 2, ncid, dimid, 0, 0.0d0, 0.0d0)
      end if

      return
      end subroutine Ncget_Dimlen
!EOC
!-------------------------------------------------------------------------

! A local version of Interp_Bilinear :

      subroutine bi_linear_interp  &
     &  (a1, a2, vinterp, v1, v2, d1, d2, data)

      implicit none

      integer, INTENT(in)  :: d1, d2  ! dimensions

      real*8,  INTENT(in)  :: a1, a2
      real*8,  INTENT(out) :: vinterp
      real*8,  INTENT(in)  :: v1(d1), v2(d2)
      real*8,  INTENT(in)  :: data(d1, d2)


      integer :: i1a, i1b, i2a, i2b
      integer :: index_of_min, index_of_max

      real*8  :: t, u

!     ----------------------------------------------------
!     Find the indices in the table where v1>a1 and v2>a2.
!     ----------------------------------------------------

!  First vector v and value a

      if (v1(1) < v1(d1)) then
        index_of_min=1
        index_of_max=d1
!  Compute indices that apply in non-extreme cases:
        i1b = Minval (Minloc (v1, mask=(v1 > a1)))   ! v1(i1b) >  a1
        i1a = i1b - 1                                ! v1(i1a) <= a1
      else
        index_of_min=d1
        index_of_max=1
!  Compute indices that apply in non-extreme cases:
        i1b = Minval (Minloc (v1, mask=(v1 > a1)))   ! v1(i1b) >  a1
        i1a = i1b + 1                                ! v1(i1a) <= a1
      end if 

      if ( a1 .LE. v1(index_of_min) ) then
        t=1.0
        i1a = index_of_min   ! used in the term that gets cancelled out
        i1b = index_of_min
      else if ( a1 .GE. v1(index_of_max) ) then
        t=1.0
        i1a = index_of_max   ! used in the term that gets cancelled out
        i1b = index_of_max
      else
        t = (a1 - v1(i1a)) / ( v1(i1b) - v1(i1a) )
      end if


!  Second vector v and value a

      if (v2(1) < v2(d2)) then
        index_of_min=1
        index_of_max=d2
!  Compute indices that apply in non-extreme cases:
        i2b = Minval (Minloc (v2, mask=(v2 > a2)))   ! v2(i2b) >  a2
        i2a = i2b - 1                                ! v2(i2a) <= a2
      else
        index_of_min=d2
        index_of_max=1
!  Compute indices that apply in non-extreme cases:
        i2b = Minval (Minloc (v2, mask=(v2 > a2)))   ! v2(i2b) >  a2
        i2a = i2b + 1                                ! v2(i2a) <= a2
      end if 

      if ( a2 .LE. v2(index_of_min) ) then
        t=1.0
        i2a = index_of_min   ! used in the term that gets cancelled out
        i2b = index_of_min
      else if ( a2 .GE. v2(index_of_max) ) then
        t=1.0
        i2a = index_of_max   ! used in the term that gets cancelled out
        i2b = index_of_max
      else
        u = (a2 - v2(i2a)) / ( v2(i2b) - v2(i2a) )
      end if

!     ----------------------------------------------------------------
!     Interpolate the data to find the value vinterp at point (a1,a2).
!     ----------------------------------------------------------------

      vinterp = (1.-t) * (1.-u) * data(i1a,i2a) +  &
     &              t  * (1.-u) * data(i1b,i2a) +  &
     &              t  *     u  * data(i1b,i2b) +  &
     &          (1.-t) *     u  * data(i1a,i2b)

      return

      end subroutine bi_linear_interp
!-------------------------------------------------------------------------
  subroutine   Finalize_GCR ()
      IMPLICIT NONE

    DEALLOCATE( GCR_EMISS )
    DEALLOCATE( GCR_YEARS )
    DEALLOCATE( GCR_LATS  )
    DEALLOCATE( GCR_SUNSPOT_COEFFS )
    DEALLOCATE( GCR_AINTCP_COEFFS  )
    DEALLOCATE( GCR_SLOPE_COEFFS   )
    DEALLOCATE( GCR_LOGPLEVS       )

    RETURN
  end subroutine   Finalize_GCR
!-------------------------------------------------------------------------
      end module gcr_mod      
