
      ! Precision for real variables that passed to ESMF states

      integer,  parameter :: r8     =  8
      integer,  parameter :: r4     =  4
      integer,  parameter :: rPrec  = r8

!     -----------------
!     Dummy I/O values.
!     -----------------

      character (len=3), parameter :: CHAR_DUM_VALUE = 'NIL'

      integer          , parameter :: INT_DUM_VALUE = -999


!     -----------------------------
!     NetCDF header field constant.
!     -----------------------------

      integer, parameter :: NETCDF_HDF = 3

!     ----------
!     Constants.
!     ----------

      integer, parameter :: GMI_MAXSTR              = 128

      integer, parameter :: MAX_LENGTH_LABELS       = 80
      integer, parameter :: MAX_NUMBER_SPECIES      = 185
      integer, parameter :: MAX_LENGTH_SITE_NAME    = 16
      integer, parameter :: MAX_LENGTH_MET_NAME     = GMI_MAXSTR !16
      integer, parameter :: MAX_LENGTH_VAR_NAME     = GMI_MAXSTR !32
      integer, parameter :: MAX_LENGTH_FILE_NAME    = 256 !128
      integer, parameter :: MAX_LENGTH_ERROR_MSG    = GMI_MAXSTR !128
      integer, parameter :: MAX_LENGTH_SPECIES_NAME = GMI_MAXSTR !16

      integer, parameter :: SIZCNSML =  255 ! small  common block size
      integer, parameter :: SIZCNMED = 1023 ! medium common block size
      integer, parameter :: SIZCNBIG = 4095 ! big    common block size


!     --------------------------------------------------------
!     The following MIN/MAX parameters need to match the range
!     used for ACTM messages in msg_numbers.h!
!     --------------------------------------------------------

      integer, parameter :: MIN_MSG_NUM = 10000 ! minimum (i.e., smallest) message number
      integer, parameter :: MAX_MSG_NUM = 14999 ! maximum (i.e., biggest)  message number

!     -------------------------------------------------------------
!     IF ANY OF THE FOLLOWING PARAMETERS IS INCREASED, BE SURE THAT
!     THE COMMON BLOCK SIZE CONSTANTS ABOVE ARE LARGE ENOUGH!
!
!     IF MAX_NUM_CONST IS INCREASED, BE SURE THAT THERE ARE ENOUGH
!     MESSAGE NUMBERS IN msg_numbers.h!
!     -------------------------------------------------------------

      integer, parameter :: MAX_NUM_CGRP    =    2 ! maximum number of chemical groups
      integer, parameter :: MAX_NUM_CONST   =  185 ! maximum number of const species
      integer, parameter :: MAX_NUM_RXNR    = 1000 ! maximum number of reaction rates
      integer, parameter :: MAX_NUM_SMARRAY =   10 ! maximum number of small array elements

      integer, parameter :: SIZGIOSML       =  255 ! small  common block size
      integer, parameter :: SIZGIOMED       = 1023 ! medium common block size
      integer, parameter :: SIZGIOBIG       = 4095 ! big    common block size
!.sds
!... needed to make this bigger when expanded number of MAX_INFILE_NAMES
!...  should equal scalar vbles+arrays+extra 128 byte character vrbls
!...  or 40 + MAX_INFILE_NAMES + MAX_COL_DIAG_SITES/5 + MAX_NUM_CONST_GIO/8
!...       + MAX_NUM_QJ + MAX_NUM_QK + (extra for expansion)
!... right now that is 40+1826+40+24+200+500+100 = 2730*128 characters
!.sds
      integer, parameter :: SIZGIOBG2       = 4095*127 ! bigger character common block size

!     ----------------------------------------------------------------
!     IF ANY OF THE FOLLOWING PARAMETERS IS INCREASED, BE SURE THAT
!     THE COMMON BLOCK SIZE CONSTANTS ABOVE ARE LARGE ENOUGH!
!
!     IF MAX_NUM_CONST_GIO IS INCREASED, BE SURE THAT THERE ARE ENOUGH
!     MESSAGE NUMBERS IN msg_numbers.h!
!     ----------------------------------------------------------------

      integer, parameter :: MAX_COL_DIAG_PRES  =   50 ! maximum number of pressures for column
                                                      ! diagnostic output
      integer, parameter :: MAX_COL_DIAG_SITES = 300 ! maximum number of column diagnostic sites
      integer, parameter :: MAX_INFILE_NAMES   = 1826 ! maximum number of input file names
      integer, parameter :: MAX_NUM_CONST_GIO  =  185 ! maximum number of const gen. I/O species
      integer, parameter :: MAX_NUM_QJ         =  200 ! maximum number of qj (photolysis) rxns
      integer, parameter :: MAX_NUM_QK         =  500 ! maximum number of qk (thermal)    rxns

      integer, parameter :: num_AerDust        = 20   ! Number of entries for aerosol/dust 
                                                      ! diagnostics
      integer, parameter :: num_CM_AerDust     =  5   ! Number of column mass entries for 
                                                      ! aerosol/dust diagnostics

      integer, parameter :: MAX_STRING_LENGTH = MAX_LENGTH_SPECIES_NAME*MAX_COL_DIAG_SITES

