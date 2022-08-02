MODULE kpp_achem_gas_Model

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Completely defines the model kpp_achem_gas
!    by using all the associated modules
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE kpp_achem_gas_Precision
  USE kpp_achem_gas_Parameters
  USE kpp_achem_gas_Global
  USE kpp_achem_gas_Function
  USE kpp_achem_gas_Integrator
  USE kpp_achem_gas_Rates
  USE kpp_achem_gas_Jacobian
  USE kpp_achem_gas_Hessian
  USE kpp_achem_gas_Stoichiom
  USE kpp_achem_gas_LinearAlgebra
  USE kpp_achem_gas_Monitor
  USE kpp_achem_gas_Util

END MODULE kpp_achem_gas_Model

