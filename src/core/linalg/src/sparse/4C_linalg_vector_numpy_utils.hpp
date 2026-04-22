// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_VECTOR_NUMPY_UTILS_HPP
#define FOUR_C_LINALG_VECTOR_NUMPY_UTILS_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#ifdef FOUR_C_WITH_PYBIND11

#include <pybind11/numpy.h>

#include <span>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * \brief Function for creating a non-constant view of a contiguous block of memory represented by
   * a std::span<T> which will be wrapped into a python-compatible pybind11::array_t<double>.
   *
   * @param[in] source The contiguous block of memory whose data should be returned in a
   * python-compatible, mutable format.
   * \return The mutable python-compatible representation of the C++ data.
   */
  template <class T>
  pybind11::array_t<T> make_numpy_view(std::span<T> values)
  {
    auto* ptr = values.data();
    auto base = pybind11::capsule(static_cast<void*>(ptr), [](void*) {});
    return pybind11::array_t<T>({static_cast<pybind11::ssize_t>(values.size())},
        {static_cast<pybind11::ssize_t>(sizeof(T))}, ptr, base);
  }

  /*!
   * \brief Function for creating a constant view of a contiguous block of memory represented by a
   * std::span<T> which will be wrapped into a python-compatible pybind11::array_t<const double>.
   *
   * @param[in] source The contiguous block of memory whose data should be returned in a
   * python-compatible, read-only format.
   * \return The read-only python-compatible representation of the C++ data.
   */
  template <class T>
  pybind11::array_t<const T> make_numpy_view(std::span<const T> values)
  {
    auto* ptr = values.data();
    auto base = pybind11::capsule(const_cast<T*>(ptr), [](void*) {});
    pybind11::array_t<const T> array({static_cast<pybind11::ssize_t>(values.size())},
        {static_cast<pybind11::ssize_t>(sizeof(T))}, ptr, base);
    array.attr("flags").attr("writeable") = false;
    return array;
  }

  /*!
   * \brief Convenience function for creating a non-constant view of a Core::LinAlg::Vector<double>
   * which will be wrapped into a python-compatible pybind11::array_t<double>.
   *
   * @param[in] source The C++ vector instance whose data should be returned in a python-compatible,
   * mutable format.
   * \return The mutable python-compatible representation of the C++ data.
   */
  template <typename T>
  pybind11::array_t<T> make_numpy_view(Vector<T>& vector)
  {
    return make_numpy_view<T>(vector.local_values_as_span());
  }

  /*!
   * \brief Convenience function for creating a constant view of a Core::LinAlg::Vector<double>
   * which will be wrapped into a python-compatible pybind11::array_t<const double>.
   *
   * @param[in] source The C++ vector instance whose data should be returned in a python-compatible,
   * read-only format.
   * \return The read-only python-compatible representation of the C++ data.
   */
  template <typename T>
  pybind11::array_t<const T> make_const_numpy_view(const Vector<T>& vector)
  {
    return make_numpy_view<T>(vector.local_values_as_span());
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
#endif
