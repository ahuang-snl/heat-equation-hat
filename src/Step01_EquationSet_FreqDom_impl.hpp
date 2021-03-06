// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Step01_EquationSet_FreqDom_impl_hpp__
#define __Step01_EquationSet_FreqDom_impl_hpp__

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesScalar.hpp"

// begin modification
#include "Panzer_Integrator_GradBasisDotVector.hpp"
// end modification

// begin HB mod
// equation set factory for time domain equation set
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"
// end HB mod

// implementation outline
// INPUT DATA: the time domain equation set
//             fundamental frequencies
//             truncation scheme: box, diamond, alpha (for the \ell^\alpha ball)
//             options: transient assisted,
// 0) this equation set should include headers for all FEM time domain equation sets
// 1) call the appropriate time domain equation set, invoking it to create DOFs and fields
//    this may essentially involve running one steady state calculation, saving the fields
//    (note that this solution can be used for the linearized problem,
//    since the steady state is the mode 0 response)
// 2) create the appropriate frequency domain fields from this information,
//    including the HB residual, HB DOF's, and HB Jacobian
// 3) solve the resulting non-linear HB system of equations

// what I need to learn how to do in order to pull this off:
// 1) call an equation set from here, to set up all the relevant fields
// 2) learn how to work with the fields created from that set-up run
// 3) define the HB residual from the residuals defined above
// 3) define the HB Jacobian
// 4) set up the non-linear system and solve it
// 5) revert to the time domain?


// ***********************************************************************

// begin HB mod
// instantiate the appropriate equation set object
PANZER_DECLARE_EQSET_TEMPLATE_BUILDER("FreqDom", user_app::EquationSet_Helmholtz,
				      EquationSet_Helmholtz)
// we can do this completely within the EquationSetFactory, instead
// end HB mod


template <typename EvalT>
user_app::EquationSet_FreqDom<EvalT>::
EquationSet_FreqDom(const Teuchos::RCP<Teuchos::ParameterList>& params,
		   const int& default_integration_order,
		   const panzer::CellData& cell_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(params,default_integration_order,cell_data,global_data,build_transient_support )
{
  // ********************
  // Validate the parameters for the chosen time domain equation set
  // ********************

  {    
    Teuchos::ParameterList valid_parameters;
    this->setDefaultValidParameters(valid_parameters);

    valid_parameters.set("Model ID","","Closure model id associated with this equaiton set");
    valid_parameters.set("Prefix","","Prefix for using multiple instatiations of thei equation set");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",-1,"Order of the integration rule");

  // begin HB mod
  // validate the FreqDom Options  
  // QUESTION: figure out why this is printed 3 times. 
    Teuchos::ParameterList& freqdom_opt = valid_parameters.sublist("FreqDom Options");  
    Teuchos::setStringToIntegralParameter<int>(
      "Time domain equation set",
      "Projection", // default gives an invalid type
      "Choose the time domain equation set to model in the frequency domain",
      Teuchos::tuple<std::string>("Helmholtz", "Projection"),
      &freqdom_opt
      );
    freqdom_opt.set("Truncation order",3,"Truncation order of the harmonic balance method.");
    // end HB mod

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  std::string model_id   = params->get<std::string>("Model ID");
  std::string prefix     = params->get<std::string>("Prefix");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order        = params->get<int>("Basis Order");
  int integration_order  = params->get<int>("Integration Order");

  // begin HB mod

  // grab the time domain equation set name
  std::string& time_domain_eqnset = params->sublist("FreqDom Options").get<std::string>("Time domain equation set");
  std::cout << "The time domain equation set we are setting up is: " << time_domain_eqnset + "." << std::endl;
  std::cout << "Do something to create multiple " + time_domain_eqnset + " equation set fields here." << std::endl;

  Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set=
    Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);
  bool found = false;

  // grab the frequency domain analysis parameters
  int truncation_order  = params->sublist("FreqDom Options").get<int>("Truncation order");

  // for now, we asume the time domain eqn set is Helmholtz
  PANZER_BUILD_EQSET_OBJECTS("FreqDom", user_app::EquationSet_Helmholtz, EquationSet_Helmholtz)
  std::cout << "Called PANZER_BUILD_EQSET_OBJECTS(\"FreqDom\", user_app::EquationSet_Helmholtz, EquationSet_Helmholtz)" << std::endl;
  std::cout << "We will print some info about what's created, here." << std::endl;
  // end HB mod

  // ********************
  // Setup DOFs and closure models
  // ********************

  {
    dof_name_ = prefix+"U";

    this->addDOF(dof_name_,basis_type,basis_order,integration_order);
    // begin modification
    this->addDOFGrad(dof_name_);
    // end modication
    // question: what other kinds of fields can be automatically created?

    // begin HB mod
    // we treat the default DOF's (named as per the time domain equation set) as the 0th mode
    // Need to add dof's for the higher order frequencies
    std::string harmonic;

    // for now, the total number of harmonics M is simply equal to the truncation order
    int M = truncation_order;

    for(int freq = 0 ; freq < M; freq++){    
        harmonic = dof_name_ + "_freq" + std::to_string(freq);
        std::cout << "Adding the " + std::to_string(freq) << "st/rd/th harmonic of the DOF." << std::endl;
	this->addDOF(harmonic, basis_type, basis_order, integration_order);
        this->addDOFGrad(harmonic);
    }
    // end HB mod

  }

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void user_app::EquationSet_FreqDom<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& fl,
				      const Teuchos::ParameterList& user_data) const
{
  // build the time domain equation set objects
  // grab the required inputs
  Teuchos::RCP<Teuchos::ParameterList> input_params = Teuchos::rcp(new Teuchos::ParameterList("User_App Parameters"));
  Teuchos::updateParametersFromXmlFile("input.xml", input_params.ptr());
  Teuchos::RCP<Teuchos::ParameterList> physics_blocks_pl   = Teuchos::rcp(new Teuchos::ParameterList(input_params->sublist("Physics Blocks")));
  Teuchos::RCP<Teuchos::ParameterList> domain_pl   = Teuchos::rcp(new Teuchos::ParameterList(physics_blocks_pl->sublist("domain")));
  Teuchos::RCP<Teuchos::ParameterList> params  = Teuchos::rcp(new Teuchos::ParameterList(domain_pl->sublist("child0")));
  Teuchos::RCP<Teuchos::ParameterList> freqdom_pl   = Teuchos::rcp(new Teuchos::ParameterList(params->sublist("FreqDom Options")));

  std::cout << "The EquationSet_FreqDom::buildAndRegisterEquationSetEvaluators() function was called!\n" 
            << "The time domain equation specified is: " << freqdom_pl->get<std::string>("Time domain equation set")
            << ". We will attempt to build its fields now." << std::endl;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // grabbing the harmonic balance options
  int truncation_order  = params->sublist("FreqDom Options").get<int>("Truncation order");
  // TODO: build and register the evaluators from the time domain equation set here
  // for now, assuming the Helmholtz equation set
  user_app::EquationSet_FreqDom<EvalT>::buildAndRegisterEquationSetEvaluators_Helmholtz(fm, fl, user_data);


  // define some special strings to use
  const std::string residual_timesfive_term     = "RESIDUAL_"+dof_name_+"_TIMESFIVE";
  // this must be satisfied by the closure model

  // ********************
  // Times five operator
  // ********************

  RCP<panzer::IntegrationRule> ir  = this->getIntRuleForDOF(dof_name_); 
  RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_name_); 

  {
    ParameterList p;
    p.set("Residual Name", residual_timesfive_term);
    p.set("Value Name",    dof_name_);
    p.set("Basis",         basis);
    p.set("IR",            ir);
    p.set("Multiplier",    0.0);

    RCP<PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));

    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian

  {
    std::vector<std::string> residual_operator_names;

    // add a bogus evaluated field (which is actually equal to 0)
    // this is the only field evaluated in this equation set
    residual_operator_names.push_back(residual_timesfive_term);

    // begin HB mod
    // trying to add a field which has an evaluator in a different equation set
    // some strings defined from the Helmholtz equation set (residual names)
    const std::string residual_projection_term     = "RESIDUAL_"+dof_name_+"_PROJECTION";
    const std::string residual_projection_src_term = "RESIDUAL_"+dof_name_+"_PROJECTION_SOURCE";
    const std::string projection_src_name = dof_name_+"_SOURCE";
    const std::string residual_laplacian_term      = "RESIDUAL_"+dof_name_+"_LAPLACIAN";
    // trying to add those evaluated fields into this residual
    residual_operator_names.push_back(residual_projection_term);
    residual_operator_names.push_back(residual_projection_src_term);
    residual_operator_names.push_back(residual_laplacian_term);
    // end HB mod

    // build a sum evaluator
    this->buildAndRegisterResidualSummationEvalautor(fm,dof_name_,residual_operator_names);

    // register residual evaluators for each harmonic
    // for now, as above, the total number of harmonics M is simply equal to the truncation order
    int M = truncation_order;
    for(int freq = 0 ; freq < M; freq++){
      this->buildAndRegisterResidualSummationEvalautor(fm,dof_name_+ "_freq" + std::to_string(freq),residual_operator_names);
      std::cout << "Adding the " + std::to_string(freq) << "st/nd/rd/th residual corresponding to harmonic DOF." << std::endl;
    }
  }

  fm.writeGraphvizFile<panzer::Traits::Residual>("graph_residual.dot");
  fm.writeGraphvizFile<panzer::Traits::Jacobian>("graph_jacobian.dot");

}

// ***********************************************************************




// *****************************************************************************
// begin HB mod
// for now, we simply copy and paste the evaluators from Step01_EquationSet_Helmholtz_impl.hpp
// in the future, we should call the buildAndRegisterEquationSetEvaluators function
// from the time domain equation set, instead of hacking it into EquationSet_FreqDom

// the goal is to minimally modify the residual evaluator from the time domain equation set
// note that the residual pushed back by the time domain equation set is NOT the complete
// 0th order harmonic residual

template <typename EvalT>
void user_app::EquationSet_FreqDom<EvalT>::
buildAndRegisterEquationSetEvaluators_Helmholtz(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& fl,
				      const Teuchos::ParameterList& user_data) const
{

  std::cout << "The EquationSet_FreqDom_impl function buildAndRegisterEquationSetEvaluators_Helmholtz was called." << std::endl;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // define some special strings to use
  const std::string residual_projection_term     = "RESIDUAL_"+dof_name_+"_PROJECTION";
  const std::string residual_projection_src_term = "RESIDUAL_"+dof_name_+"_PROJECTION_SOURCE";
  
  const std::string projection_src_name = dof_name_+"_SOURCE";
    // this must be satisfied by the closure model

  // begin modification
  const std::string residual_laplacian_term      = "RESIDUAL_"+dof_name_+"_LAPLACIAN";
  // end modification


  // ********************
  // Helmholtz Equation
  // ********************

  RCP<panzer::IntegrationRule> ir  = this->getIntRuleForDOF(dof_name_); 
  RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_name_); 

  // Projection operator (U,phi)
  {
    ParameterList p;
    p.set("Residual Name", residual_projection_term);
    p.set("Value Name",    dof_name_);
    p.set("Basis",         basis);
    p.set("IR",            ir);
    p.set("Multiplier",    1.0);

    RCP<PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Source operator -(u_source,phi)
  {
    ParameterList p;
    p.set("Residual Name", residual_projection_src_term);
    p.set("Value Name",    projection_src_name);
    p.set("Basis",         basis);
    p.set("IR",            ir);
    p.set("Multiplier",    -1.0);

    RCP<PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(p));
    
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // begin modification
  // Laplacian operator (grad u , grad basis)
  {
    ParameterList p;
    p.set("Residual Name", residual_laplacian_term);
    p.set("Flux Name", "GRAD_"+dof_name_);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }
  // end modification
  // note that we do not have to explicitly evaluate a "GRAD_"+dof_name_ field

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian

  // don't push back the residual here; instead, push it back in the FreqDom's buildAndRegisterEvalautors()
  /*
  {
    std::vector<std::string> residual_operator_names;

    residual_operator_names.push_back(residual_projection_term);
    residual_operator_names.push_back(residual_projection_src_term);

    // begin modification
    residual_operator_names.push_back(residual_laplacian_term);
    // end modification

    // build a sum evaluator
    this->buildAndRegisterResidualSummationEvalautor(fm,dof_name_,residual_operator_names);
  }
  */
}

// end HB mod
// *********************************************************************************

#endif
