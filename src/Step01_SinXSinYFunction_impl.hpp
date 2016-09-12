#ifndef __Step01_SinXSinYFunction_impl_hpp__
#define __Step01_SinXSinYFunction_impl_hpp__

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace user_app {

//**********************************************************************
template <typename EvalT,typename Traits>
SinXSinYFunction<EvalT,Traits>::SinXSinYFunction(const std::string & name,
                                             double xperiod,double yperiod,
                                             const panzer::IntegrationRule & ir)
  : xperiod_(xperiod) 
  , yperiod_(yperiod) 
  , ir_degree_(ir.cubature_degree)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;

  result = PHX::MDField<ScalarT,panzer::Cell,panzer::Point>(name, data_layout);
  this->addEvaluatedField(result);

  this->setName("SinXSinY Function("+name+")");
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SinXSinYFunction<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
                                                          PHX::FieldManager<Traits>& fm)
{

  this->utils.setFieldData(result,fm);

  ir_index_ = panzer::getIntegrationRuleIndex(ir_degree_,(*sd.worksets_)[0]);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SinXSinYFunction<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  for (panzer::index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int point = 0; point < result.extent_int(1); ++point) {

      const double& x = workset.int_rules[ir_index_]->ip_coordinates(cell,point,0);
      const double& y = workset.int_rules[ir_index_]->ip_coordinates(cell,point,1);

      result(cell,point) = sin(2*M_PI*xperiod_*x)*sin(2*M_PI*yperiod_*y);
    }
  }
}

//**********************************************************************
}

#endif
