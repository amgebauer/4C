// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_scatra_wrapper.hpp"

#include "4C_linalg_sparseoperator.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::AdapterScatraWrapper::AdapterScatraWrapper(Teuchos::RCP<ScatraInterface> scatra)
    : scatra_timint_(scatra)
{
  Teuchos::rcp_dynamic_cast<ScaTra::ScaTraTimIntImpl>(scatra_timint_, true)
      ->set_model_evaluatro_ptr(this);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::AdapterScatraWrapper::evaluate_additional_solution_depending_models(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,
    Teuchos::RCP<Core::LinAlg::Vector<double>> rhs)
{
  // do nothing so far
  return;
}

FOUR_C_NAMESPACE_CLOSE
