
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Mike Wehner, et.al., LLNL)
!   jrt@llnl.gov
!
! FILE
!   gem_msg_numbers.h
!
! DESCRIPTION
!   This file contains the message numbers used to identify specific
!   messages in the message passing calls.
!
!   The convention for numbering messages in the the LLNL ESM model is:
!       1 -   99 ==> general ESM messages (including intermodule messages)
!     100 - 1999 ==> ACTM messages
!
! HISTORY
!  - May 10, 2006 - Bigyani Das
!    Added message numbers for relative humidity, gridbox height, and mass
!    to be output in frequency output files
!  - March 23, 2006 - Bigyani Das
!    Added message numbers for psf and kel to be output in noon file:
!    ACTM_SG_KEL_NOON_NC, ACTM_SG_PSF_NOON_NC
!  - December 8, 2005 - Bigyani Das
!    Added message numbers for 4 output frequencies such as 
!    ACTM_SG_FREQ11_NC, ACTM_SG_FREQ12_NC, ACTM_SG_FREQ1C1_NC, 
!    ACTM_SG_FREQ1C2_NC,  ACTM_SG_KEL_FREQ1_NC, ACTM_SG_PSF_FREQ1_NC 
!    * Removed message numbers related to daily and hourly outputs
!=============================================================================


!     ====================================
!     Start of ESM message numbers (1-99).
!     ====================================

!     --------------------------------------------------------------------
!     Message numbers used for subdomain to subdomain communication of the
!     border data.
!
!     WARNING: Because the following border communication parameters are
!              added to other base values, the values of these parameters
!              must remain at 1,2,3,4 resp.  DO NOT CHANGE THEM.
!     --------------------------------------------------------------------

      integer, parameter :: WE_ESM = 1
      integer, parameter :: EW_ESM = 2
      integer, parameter :: SN_ESM = 3
      integer, parameter :: NS_ESM = 4


!     ---------------------
!     ESM control messages.
!     ---------------------

      integer, parameter :: ESM_CNTL_C  = 10
      integer, parameter :: ESM_CNTL_I  = 12
      integer, parameter :: ESM_CNTL_R  = 14


!     --------------------
!     ESM timing messages.
!     --------------------

      integer, parameter :: ESM_TIMER_I = 20
      integer, parameter :: ESM_TIMER_R = 22


!     ----------------------
!     For ESM processor map.
!     ----------------------

      integer, parameter :: PROC_ESM1   = 30
      integer, parameter :: PROC_ESM2   = 32


!     =========================================
!     Start of ACTM message numbers (100-1999).
!     =========================================

!     ------------------------
!     Control message numbers.
!     ------------------------

      integer, parameter :: ACTM_DM_L1   = 100
      integer, parameter :: ACTM_DM_I1   = 102
      integer, parameter :: ACTM_DM_I2   = 104
      integer, parameter :: ACTM_DM_I3   = 106
      integer, parameter :: ACTM_DM_I4   = 108
      integer, parameter :: ACTM_SD_I1   = 110

      integer, parameter :: ACTM_H2OCOND = 120
      integer, parameter :: ACTM_MCOR    = 122
      integer, parameter :: ACTM_RAREA   = 124

      integer, parameter :: ACTM_KEL     = 130
      integer, parameter :: ACTM_DKEL    = 132
      integer, parameter :: ACTM_PCTM1   = 134
      integer, parameter :: ACTM_PCTM2   = 136
      integer, parameter :: ACTM_PSX     = 138
      integer, parameter :: ACTM_DPSX    = 140
      integer, parameter :: ACTM_XMASS   = 142
      integer, parameter :: ACTM_YMASS   = 144

      integer, parameter :: ACTM_EMISS_AERO   = 150
      integer, parameter :: ACTM_EMISS_DUST   = 152

      integer, parameter :: ACTM_UVALBDATA    = 160

      integer, parameter :: ACTM_SFALBDATA(4) =  &
     &  (/ 170, 172, 174, 176 /)


!     ------------------------------------
!     Subdomain to Global message numbers.
!     ------------------------------------

      integer, parameter :: ACTM_SG_CONST_ASC    = 200

      integer, parameter :: ACTM_SG_MAX_ASC(3)   =  &
     &  (/ 210, 212, 214 /)

      integer, parameter :: ACTM_SG_MIN_ASC(3)   =  &
     &  (/ 220, 222, 224 /)


      integer, parameter :: ACTM_SG_AFLUX1_NC    = 230
      integer, parameter :: ACTM_SG_FLUX1_NC     = 232
      integer, parameter :: ACTM_SG_FLUX2_NC     = 234

      integer, parameter :: ACTM_SG_CONST1_NC    = 240
      integer, parameter :: ACTM_SG_CONST2_NC    = 242
      integer, parameter :: ACTM_SG_NOON1_NC     = 244
      integer, parameter :: ACTM_SG_NOON2_NC     = 246
      integer, parameter :: ACTM_SG_LOCAL1_NC     = 2444
      integer, parameter :: ACTM_SG_LOCAL2_NC     = 2446
      integer, parameter :: ACTM_SG_LOCALC1_NC     = 2448
      integer, parameter :: ACTM_SG_LOCALC2_NC     = 2450

      integer, parameter :: ACTM_SG_QJ1_NC       = 248
      integer, parameter :: ACTM_SG_QJ2_NC       = 250
      integer, parameter :: ACTM_SG_QK1_NC       = 252
      integer, parameter :: ACTM_SG_QK2_NC       = 254
      integer, parameter :: ACTM_SG_QQJ1_NC      = 256
      integer, parameter :: ACTM_SG_QQJ2_NC      = 258
      integer, parameter :: ACTM_SG_QQK1_NC      = 260
      integer, parameter :: ACTM_SG_QQK2_NC      = 262
      integer, parameter :: ACTM_SG_SAD1_NC      = 264
      integer, parameter :: ACTM_SG_SAD2_NC      = 266
      integer, parameter :: ACTM_SG_YDA1_NC      = 268
      integer, parameter :: ACTM_SG_YDA2_NC      = 270

      integer, parameter :: ACTM_SG_OptDepth_NC  = 272

      integer, parameter :: ACTM_SG_CM_AerDust   = 273
      integer, parameter :: ACTM_SG_GridHeight   = 274

      integer, parameter :: ACTM_SG_DMS_OH       = 290
      integer, parameter :: ACTM_SG_DMS_NO3      = 292
      integer, parameter :: ACTM_SG_SO2_OH       = 294
      integer, parameter :: ACTM_SG_SO2_H2O2     = 296
      integer, parameter :: ACTM_SG_SO2_O3       = 298

      integer, parameter :: ACTM_SG_DRY          = 300
      integer, parameter :: ACTM_SG_WET          = 302
      integer, parameter :: ACTM_SG_EMISS        = 304
      integer, parameter :: ACTM_SG_EMISS2       = 3042
      integer, parameter :: ACTM_SG_EMISS3D_NC   = 3044
      integer, parameter :: ACTM_SG_EMISS3D2_NC  = 3046

      integer, parameter :: ACTM_SG_H2OBACK_NC   = 310
      integer, parameter :: ACTM_SG_H2OCOND_NC   = 312
      integer, parameter :: ACTM_SG_HNO3COND_NC  = 314
      integer, parameter :: ACTM_SG_HNO3GAS_NC   = 316
      integer, parameter :: ACTM_SG_KEL_NC       = 318
      integer, parameter :: ACTM_SG_MASS_NC      = 320
      integer, parameter :: ACTM_SG_METWATER_NC  = 322
      integer, parameter :: ACTM_SG_PSF_NC       = 324
      integer, parameter :: ACTM_SG_TROPP_NC     = 325
      integer, parameter :: ACTM_SG_REFFICE_NC   = 326
      integer, parameter :: ACTM_SG_PV_NC        = 327
      integer, parameter :: ACTM_SG_REFFSTS_NC   = 328
      integer, parameter :: ACTM_SG_VFALL_NC     = 330

      integer, parameter :: ACTM_SG_H2OCOND_RST  = 340
      integer, parameter :: ACTM_SG_CONST1_RST   = 342
      integer, parameter :: ACTM_SG_CONST2_RST   = 344

      integer, parameter :: ACTM_SG_NCUMT1_NC    = 370
      integer, parameter :: ACTM_SG_NCUMT2_NC    = 372

      integer, parameter :: ACTM_SG_lightn_NC    = 374
      integer, parameter :: ACTM_SG_dtrn_NC      = 376
      integer, parameter :: ACTM_SG_cmf_NC       = 378
      integer, parameter :: ACTM_SG_flashr_NC    = 380
      integer, parameter :: ACTM_SG_cmiflg_NC    = 382

      integer, parameter :: ACTM_SG_KEL_NOON_NC       = 7318
      integer, parameter :: ACTM_SG_PSF_NOON_NC       = 7324
!     -----------------------------------------------------------
!     Boundary Condition message numbers.
!
!     BE SURE THAT ALL "ACTM_BC_" MESSAGE NUMBERS ARE AT LEAST 10
!     VALUES APART.
!     -----------------------------------------------------------

      integer, parameter :: ACTM_BC_CRX    = 400
      integer, parameter :: ACTM_BC_CRY    = 410
      integer, parameter :: ACTM_BC_DCX    = 420
      integer, parameter :: ACTM_BC_DCY    = 430
      integer, parameter :: ACTM_BC_FX     = 440
      integer, parameter :: ACTM_BC_PU     = 450
      integer, parameter :: ACTM_BC_PVADV2 = 460
      integer, parameter :: ACTM_BC_QQ1    = 470
      integer, parameter :: ACTM_BC_QQU    = 480
      integer, parameter :: ACTM_BC_QQV    = 490
      integer, parameter :: ACTM_BC_UA     = 500
      integer, parameter :: ACTM_BC_VA     = 510


!     --------------------------------
!     Large arrays of message numbers.
!     --------------------------------

      integer :: actm_val


      integer, parameter :: ACTM_GS_4D1(200) =  &
     &  (/ (actm_val, actm_val = 600, 799) /)

!
!     For aerosol/dust computations
!

      integer :: actm_val3, actm_val4, actm_val5

      integer, parameter :: ACTM_DAERSL(2) =  &
     &  (/ (actm_val3, actm_val3 = 810, 811) /)
      integer, parameter :: ACTM_WAERSL(5) =  &
     &  (/ (actm_val4, actm_val4 = 820, 824) /)
      integer, parameter :: ACTM_DUST(7)   =  &
     &  (/ (actm_val5, actm_val5 = 830, 836) /)
!
!     GOCART-GMI
!
      integer, parameter :: ACTM_erod(1)  = 840
      integer, parameter :: ACTM_erod_mod = 846
      integer, parameter :: ACTM_water    = 847

!    ==============================
!    For output frequency
!    ======================
      integer, parameter :: ACTM_SG_FREQ11_NC     = 1452
      integer, parameter :: ACTM_SG_FREQ12_NC     = 1454
      integer, parameter :: ACTM_SG_FREQ1C1_NC     = 1460
      integer, parameter :: ACTM_SG_FREQ1C1S1_NC     = 1462
      integer, parameter :: ACTM_SG_KEL_FREQ1_NC   = 1318
      integer, parameter :: ACTM_SG_PSF_FREQ1_NC    =1324
      integer, parameter :: ACTM_SG_MASS_FREQ1_NC    =1326
      integer, parameter :: ACTM_SG_GRID_HEIGHT_FREQ1_NC    =1328
      integer, parameter :: ACTM_SG_REL_HUM_FREQ1_NC    =1330
      integer, parameter :: ACTM_SG_METWATER_FREQ1_NC    =1464

!    ======================
      integer, parameter :: ACTM_SG_FREQ21_NC     = 3452
      integer, parameter :: ACTM_SG_FREQ22_NC     = 3454
      integer, parameter :: ACTM_SG_FREQ2C1_NC     = 3460
      integer, parameter :: ACTM_SG_FREQ2C1S1_NC     = 3462
      integer, parameter :: ACTM_SG_KEL_FREQ2_NC   = 3318
      integer, parameter :: ACTM_SG_PSF_FREQ2_NC    =3324
      integer, parameter :: ACTM_SG_MASS_FREQ2_NC    =3326
      integer, parameter :: ACTM_SG_GRID_HEIGHT_FREQ2_NC    =3328
      integer, parameter :: ACTM_SG_REL_HUM_FREQ2_NC    =3330
      integer, parameter :: ACTM_SG_METWATER_FREQ2_NC    =3464

!    ======================
      integer, parameter :: ACTM_SG_FREQ31_NC     = 4452
      integer, parameter :: ACTM_SG_FREQ32_NC     = 4454
      integer, parameter :: ACTM_SG_FREQ3C1_NC     = 4460
      integer, parameter :: ACTM_SG_FREQ3C1S1_NC     = 4462
      integer, parameter :: ACTM_SG_KEL_FREQ3_NC   = 4318
      integer, parameter :: ACTM_SG_PSF_FREQ3_NC    =4324
      integer, parameter :: ACTM_SG_MASS_FREQ3_NC    =4326
      integer, parameter :: ACTM_SG_GRID_HEIGHT_FREQ3_NC    =4328
      integer, parameter :: ACTM_SG_REL_HUM_FREQ3_NC    =4330
      integer, parameter :: ACTM_SG_METWATER_FREQ3_NC    =4464

!    ======================
      integer, parameter :: ACTM_SG_FREQ41_NC     = 5452
      integer, parameter :: ACTM_SG_FREQ42_NC     = 5454
      integer, parameter :: ACTM_SG_FREQ4C1_NC     = 5460
      integer, parameter :: ACTM_SG_FREQ4C1S1_NC     = 5462
      integer, parameter :: ACTM_SG_KEL_FREQ4_NC   = 5318
      integer, parameter :: ACTM_SG_PSF_FREQ4_NC    =5324
      integer, parameter :: ACTM_SG_MASS_FREQ4_NC    =5326
      integer, parameter :: ACTM_SG_GRID_HEIGHT_FREQ4_NC    =5328
      integer, parameter :: ACTM_SG_REL_HUM_FREQ4_NC    =5328
      integer, parameter :: ACTM_SG_METWATER_FREQ4_NC    =5464

!.sds for flux output    ======================
      integer, parameter :: ACTM_SG_FLUXm_NC     = 6450
      integer, parameter :: ACTM_SG_FLUXx_NC     = 6452
      integer, parameter :: ACTM_SG_FLUXy_NC     = 6454
      integer, parameter :: ACTM_SG_FLUXz_NC     = 6456
      integer, parameter :: ACTM_SG_FLUXx1_NC    = 6462
      integer, parameter :: ACTM_SG_FLUXy1_NC    = 6464
      integer, parameter :: ACTM_SG_FLUXz1_NC    = 6466
      integer, parameter :: ACTM_SG_PSF_FLUX_NC  = 6468
