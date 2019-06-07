      module module_data_mosaic_cloud

      use module_data_mosaic_kind, only:  r8

      implicit none

      integer, parameter :: nrxn_cld = 2

!------------------------------------------------------------------------
      integer, save ::   &
       iso4_c,      ino3_c,      icl_c,       inh4_c,      ioc_c,   &
       imsa_c,      ico2_c,      ina_c,       ica_c,       ibc_c,   &
       ioin_c,      iso2_c,      ihono_c,     ih2o2_c,     ich3ooh_c,   &
       ihcooh_c,    ircooh_c,    ihcho_c,     io3_c,       iho2_c,   &
       ino2_c,      ino3r_c,     in2o5_c

      integer, save ::   &
       jh_c,        jnh4_c,      jna_c,       jhso4_c,     jso4_c,   &
       jno3_c,      jcl_c,       jno2_c,      jhso3_c,     jso3_c,   &
       jhco3_c,     jco3_c,      jho2_c,      jhcoo_c,     jrcoo_c,   &
       jmsa_c,      joh_c,       jch2oh2_c

!------------------------------------------------------------------------

      end module module_data_mosaic_cloud
