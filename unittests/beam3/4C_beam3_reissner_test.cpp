// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_beam3_reissner.hpp"

#include "4C_fem_general_element.hpp"

#include <array>

const double testTolerance = 1e-14;

namespace
{
  using namespace FourC;

  class Beam3r : public ::testing::Test
  {
   public:
    Beam3r()
    {
      testdis_ = std::make_shared<Core::FE::Discretization>("Beam3r", MPI_COMM_WORLD, 3);

      std::vector<std::vector<double>> xrefe{{-0.05, 0.05, 0.3}, {0.45, -0.05, 0.1}};
      std::vector<double> xrefe_full{-0.05, 0.05, 0.3, 0.45, -0.05, 0.1};

      for (int lid = 0; lid < 2; ++lid)
        testdis_->add_node(std::make_shared<Core::Nodes::Node>(lid, xrefe[lid], 0));

      testele_ = std::make_shared<Discret::Elements::Beam3r>(0, 0);
      std::array<int, 2> node_ids{0, 1};
      testele_->set_node_ids(2, node_ids.data());

      // create 1 element discretization
      testdis_->add_element(testele_);
      testdis_->fill_complete(false, false, false);

      // setup internal beam element parameters
      std::vector<double> rotrefe(9);
      rotrefe[0] = -2.135698785951414;
      rotrefe[1] = -1.1055190408131161;
      rotrefe[2] = -0.45792098016648797;
      rotrefe[3] = 0.09071600605476587;
      rotrefe[4] = -0.31314870676006484;
      rotrefe[5] = -0.5590172175309829;
      rotrefe[6] = -0.44757433200569813;
      rotrefe[7] = -0.14845112617443665;
      rotrefe[8] = -0.628849061811312;

      testele_->set_centerline_hermite(true);
      testele_->set_up_reference_geometry<3, 2, 2>(xrefe_full, rotrefe);
    }

   protected:
    //! dummy discretization for holding element and node pointers
    std::shared_ptr<Core::FE::Discretization> testdis_;
    //! the beam3r element to be tested
    std::shared_ptr<Discret::Elements::Beam3r> testele_;
  };

  /**
   * Test reference length calculation of Simo-Reissner beam
   */
  TEST_F(Beam3r, RefLength)
  {
    EXPECT_NEAR(testele_->ref_length(), 0.61920435714496047, testTolerance);
  }

}  // namespace
