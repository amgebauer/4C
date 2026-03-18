// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_config.hpp"

#include "4C_structure_new_predict_python_wrapper_utils.hpp"

#ifdef FOUR_C_WITH_PYBIND11

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
pybind11::module_ Solid::Predict::PythonWrapperUtils::ensure_python_wrapper_module()
{
  static pybind11::module_ module = []()
  {
    auto* def = new PyModuleDef{PyModuleDef_HEAD_INIT, "fourc_predict", nullptr, -1, nullptr,
        nullptr, nullptr, nullptr, nullptr};

    pybind11::module_ m = pybind11::module_::create_extension_module("fourc_predict", nullptr, def);

    pybind11::class_<Solid::Predict::PythonWrapperUtils::PythonMesh>(m, "GatheredMesh")
        .def(pybind11::init<>())
        .def_readwrite("node_ids", &Solid::Predict::PythonWrapperUtils::PythonMesh::node_ids)
        .def_readwrite("coordinates", &Solid::Predict::PythonWrapperUtils::PythonMesh::coordinates)
        .def_readwrite(
            "disp_dof_ids", &Solid::Predict::PythonWrapperUtils::PythonMesh::disp_dof_ids)
        .def_readwrite("element_ids", &Solid::Predict::PythonWrapperUtils::PythonMesh::element_ids)
        .def_readwrite(
            "connectivity", &Solid::Predict::PythonWrapperUtils::PythonMesh::connectivity);

    pybind11::class_<Solid::Predict::PythonWrapperUtils::PythonContext>(m, "Context")
        .def(pybind11::init<>())
        .def_readwrite("input_file", &Solid::Predict::PythonWrapperUtils::PythonContext::input_file)
        .def_readwrite(
            "problem_dim", &Solid::Predict::PythonWrapperUtils::PythonContext::problem_dim)
        .def_readwrite("step_max", &Solid::Predict::PythonWrapperUtils::PythonContext::step_max)
        .def_readwrite("time_max", &Solid::Predict::PythonWrapperUtils::PythonContext::time_max)
        .def_readwrite("mesh", &Solid::Predict::PythonWrapperUtils::PythonContext::mesh);

    pybind11::class_<Solid::Predict::PythonWrapperUtils::PythonState>(m, "State")
        .def(pybind11::init<>())
        .def_readwrite("step", &Solid::Predict::PythonWrapperUtils::PythonState::step)
        .def_readwrite("time", &Solid::Predict::PythonWrapperUtils::PythonState::time)
        .def_readwrite("dt", &Solid::Predict::PythonWrapperUtils::PythonState::dt)
        .def_readwrite(
            "num_global_dofs", &Solid::Predict::PythonWrapperUtils::PythonState::num_global_dofs)
        .def_readwrite("dis_n", &Solid::Predict::PythonWrapperUtils::PythonState::dis_n)
        .def_readwrite("dis_np", &Solid::Predict::PythonWrapperUtils::PythonState::dis_np)
        .def_readwrite("vel_n", &Solid::Predict::PythonWrapperUtils::PythonState::vel_n)
        .def_readwrite("vel_np", &Solid::Predict::PythonWrapperUtils::PythonState::vel_np)
        .def_readwrite("acc_n", &Solid::Predict::PythonWrapperUtils::PythonState::acc_n)
        .def_readwrite("acc_np", &Solid::Predict::PythonWrapperUtils::PythonState::acc_np);

    pybind11::module_::import("sys").attr("modules")["fourc_predict"] = m;
    return m;
  }();

  return module;
}

FOUR_C_NAMESPACE_CLOSE

#endif