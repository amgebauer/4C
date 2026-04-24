// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_PORO_ELE_UTILS_HPP
#define FOUR_C_SOLID_PORO_ELE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_scatra_input.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{
  constexpr auto get_supported_impl_types()
  {
    return std::array{ScaTra::ImplType::impltype_advreac,
        ScaTra::ImplType::impltype_cardiac_monodomain, ScaTra::ImplType::impltype_chemo,
        ScaTra::ImplType::impltype_chemoreac, ScaTra::ImplType::impltype_loma,
        ScaTra::ImplType::impltype_poro, ScaTra::ImplType::impltype_pororeac,
        ScaTra::ImplType::impltype_pororeacECM, ScaTra::ImplType::impltype_multipororeac,
        ScaTra::ImplType::impltype_std, ScaTra::ImplType::impltype_undefined};
  }

  inline std::map<std::string, ScaTra::ImplType> get_impltype_inpar_map()
  {
    constexpr auto supported_impl_types = get_supported_impl_types();
    std::map<std::string, ScaTra::ImplType> impltype_map;

    for (const auto& impltype : supported_impl_types)
    {
      impltype_map[ScaTra::impltype_to_string(impltype)] = impltype;
    }

    return impltype_map;
  }

}  // namespace Discret::Elements

FOUR_C_NAMESPACE_CLOSE

#endif