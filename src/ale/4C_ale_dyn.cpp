// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_ale_dyn.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_ale_resulttest.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void dyn_ale_drt()
{
  Global::Problem* problem = Global::Problem::instance();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  std::shared_ptr<Core::FE::Discretization> actdis = problem->get_dis("ale");

  // -------------------------------------------------------------------
  // ask ALE::AleBaseAlgorithm for the ale time integrator
  // -------------------------------------------------------------------
  Adapter::AleBaseAlgorithm ale(*problem, problem->ale_dynamic_params(), actdis);
  std::shared_ptr<Adapter::Ale> aletimint = ale.ale_field();

  // -------------------------------------------------------------------
  // read the restart information, set vectors and variables if necessary
  // -------------------------------------------------------------------
  const int restart = problem->restart();
  if (restart) aletimint->read_restart(restart);

  // -------------------------------------------------------------------
  // call time loop
  // -------------------------------------------------------------------
  aletimint->create_system_matrix();
  aletimint->integrate();

  // -------------------------------------------------------------------
  // do the result test
  // -------------------------------------------------------------------
  // test results
  problem->add_field_test(aletimint->create_field_test());
  problem->test_all(actdis->get_comm());

  return;
}

FOUR_C_NAMESPACE_CLOSE
