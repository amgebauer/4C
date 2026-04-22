// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_PREDICT_PYTHON_WRAPPER_HPP
#define FOUR_C_STRUCTURE_NEW_PREDICT_PYTHON_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_structure_new_predict_generic.hpp"

#include <filesystem>

#ifdef FOUR_C_WITH_PYBIND11

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace Predict
  {
    class PythonWrapper : public Generic
    {
     public:
      //! constructor
      PythonWrapper();

      /*!
       * \brief Default destructor
       *
       * The destructor is defaulted to ensure proper cleanup of the python_implementation_ unique
       * pointer. Hence, it must be defined in the source file where the complete definition of the
       * Implementation class is available.
       */
      ~PythonWrapper() override;

      //! Performs class specific setup instructions once
      void setup() override;

      //! Performs the class specific predictor step
      void compute(::NOX::Abstract::Group& grp) override;

     private:
      //! File name of python file defining the surrogate predictor
      std::filesystem::path python_filename_;

      //! Forward declaration of the actual class handling the python implementation
      class Implementation;
      //! Pointer to the python implementation of the surrogate predictor
      std::unique_ptr<Implementation> python_implementation_;
    };  // class PythonWrapper
  }  // namespace Predict
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
#endif