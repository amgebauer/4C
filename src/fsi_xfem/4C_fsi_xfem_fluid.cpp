// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fsi_xfem_fluid.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_io_pstream.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidXFEMAlgorithm::FluidXFEMAlgorithm(MPI_Comm comm, Global::Problem& problem)
    : FluidMovingBoundaryBaseAlgorithm(problem, problem.fluid_dynamic_params(), "FSICoupling"),
      comm_(comm),
      problem_(problem)
{
  const Teuchos::ParameterList& fluiddyn = problem_.fluid_dynamic_params();

  step_ = 0;
  time_ = 0.;
  dt_ = fluiddyn.get<double>("TIMESTEP");
  nstep_ = fluiddyn.get<int>("NUMSTEP");
  maxtime_ = fluiddyn.get<double>("MAXTIME");
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::timeloop()
{
  if (problem_.get_problem_type() == Core::ProblemType::fluid_xfem)
  {
    if (Core::Communication::my_mpi_rank(get_comm()) == 0)
      std::cout << "Integrate routine for MOVING INTERFACES"
                << "\n"
                << std::endl;


    while (not_finished())
    {
      prepare_time_step();
      solve();
      update();
      output();
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::read_restart(int step)
{
  time_ = mb_fluid_field()->read_restart(step);
  step_ = step;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::prepare_time_step()
{
  step_ += 1;
  time_ += dt_;


  mb_fluid_field()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::solve() { mb_fluid_field()->nonlinear_solve(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::update() { mb_fluid_field()->update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::output() { mb_fluid_field()->output(); }

FOUR_C_NAMESPACE_CLOSE
