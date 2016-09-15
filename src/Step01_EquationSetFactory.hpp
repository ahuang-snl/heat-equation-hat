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

#ifndef __Step01_EquationSetFactory_hpp__
#define __Step01_EquationSetFactory_hpp__

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

// Add my equation sets here
#include "Step01_EquationSet_Projection.hpp"

// begin modification
#include "Step01_EquationSet_Helmholtz.hpp"
// end modification

// begin HB mod
#include "Step01_EquationSet_FreqDom.hpp"
// end HB mod


namespace user_app {

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER("Projection", user_app::EquationSet_Projection,
					EquationSet_Projection)

  // begin modificatin
  // A macro that defines a class to make construction of the equation sets easier
  //   - The equation set is constructed over a list of automatic differention types
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER("Helmholtz", user_app::EquationSet_Helmholtz,
					EquationSet_Helmholtz)
  // end modification
  
  // begin HB mod
  // this provides a placeholder template, since the user_app::EquationSet_FreqDom
  // (defined in Step01_EquationSet_FreqDom_impl.hpp) can be used to grab the info
  // from the appropriate time domain equation set. The problem is that
  // we cannot grab the name of the time domain equation set at this point.
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER("FreqDom", user_app::EquationSet_FreqDom,
  					EquationSet_FreqDom)
  // end HB mod

  class EquationSetFactory : public panzer::EquationSetFactory {

  public:

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		     const int& default_integration_order,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     const bool build_transient_support) const
    {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set= 
	Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);
      
      bool found = false;

      PANZER_BUILD_EQSET_OBJECTS("Projection", user_app::EquationSet_Projection,
				 EquationSet_Projection)
      
      // begin modification
      // macro which builds the objects in the Helmholtz equation set
      // macro checks if(ies.name=="Helmholtz") then an EquationSet_Helmholtz object is constructed
      PANZER_BUILD_EQSET_OBJECTS("Helmholtz", user_app::EquationSet_Helmholtz,
				 EquationSet_Helmholtz)
      // end modification
      // question: what does this do? if it's required to be done when 
      //           the equations are being built anyway, why isn't this done by 
      //           the PANZER_DECLARE_EQSET_TEMPLATE_BUILDER?


     // begin HB mod
     // check that we are using the "FreqDom" equation set
     // check for presence of a "FreqDom Options" sublist (parallel to the "FreqDom" equation set specification
     // error if "FreqDom" sublist is missing
     // grab the time domain equation set name in time_domain_eqnset 
     if(params->get<std::string>("Type") == "FreqDom"){

        std::cout << "A frequency domain analysis is specified." << std::endl;
 
        if(params->isSublist("FreqDom Options")){

          // grab the time domain equation set here
	  std::cout << "Found a FreqDom Options sublist!" << std::endl;
	  std::string& time_domain_eqnset   = params->sublist("FreqDom Options").get<std::string>("Time domain equation set");
	  std::cout << "The time domain equation set you chose is: " + time_domain_eqnset << std::endl;

	PANZER_BUILD_EQSET_OBJECTS("FreqDom", user_app::EquationSet_FreqDom,EquationSet_FreqDom) 
        // alternatively, use the expanded macro, doesn't check key
	// EquationSet_FreqDom_TemplateBuilder builder(params, default_integration_order, cell_data, global_data, build_transient_support);
        // eq_set->buildObjects(builder);

        } 
        else if(!params->isSublist("FreqDom Options")) {
 
	   // error out if the "FreqDom Options" are missing 
      	  std::string msg = "Error - Equation set \"FreqDom\" chosen, but missing a \"FreqDom Options\" parameter sublist!.\n";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);            

        } 

	found = true;
     }
     // end HB mod


      if (!found) {
	std::string msg = "Error - the \"Equation Set\" with \"Type\" = \"" + params->get<std::string>("Type") +
	  "\" is not a valid equation set identifier. Please supply the correct factory.\n";
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
      
      return eq_set;
    }
    
  };

}

#endif
