
      module GmiFastJX_includeMod

      private
      PUBLIC  :: t_fastJXbundle

# include "parm_MIE_fastJX65.h"

TYPE t_fastJXbundle
!... solar cycle capability

       REAL*8 :: fjx_solar_cycle_param(W_)

!----------------------------------------------------------------------------

END TYPE t_fastJXbundle
      end module GmiFastJX_includeMod
