// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_serialdensesolver.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  void SerialDenseSolver::set_matrix(SerialDenseMatrix& matrix)
  {
    solver_.setMatrix(Teuchos::rcpFromRef(matrix.base()));
  }

  void SerialDenseSolver::set_vectors(SerialDenseVector& lhs, SerialDenseVector& rhs)
  {
    solver_.setVectors(Teuchos::rcpFromRef(lhs.base()), Teuchos::rcpFromRef(rhs.base()));
  }

  void SerialDenseSolver::set_vectors(SerialDenseMatrix& lhs, SerialDenseMatrix& rhs)
  {
    solver_.setVectors(Teuchos::rcpFromRef(lhs.base()), Teuchos::rcpFromRef(rhs.base()));
  }

  void SerialDenseSolver::factor_with_equilibration(bool enable)
  {
    solver_.factorWithEquilibration(enable);
  }

  void SerialDenseSolver::solve_to_refined_solution(bool enable)
  {
    solver_.solveToRefinedSolution(enable);
  }

  int SerialDenseSolver::factor() { return solver_.factor(); }

  int SerialDenseSolver::solve() { return solver_.solve(); }

  int SerialDenseSolver::invert() { return solver_.invert(); }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE
