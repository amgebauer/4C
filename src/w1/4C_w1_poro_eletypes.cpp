// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_w1_poro_eletypes.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_w1_poro.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                       |
 *----------------------------------------------------------------------*/

Discret::Elements::WallQuad4PoroType Discret::Elements::WallQuad4PoroType::instance_;

Discret::Elements::WallQuad4PoroType& Discret::Elements::WallQuad4PoroType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::WallQuad4PoroType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::Wall1Poro<Core::FE::CellType::quad4>(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallQuad4PoroType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ4PORO")
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::quad4>>(id, owner);
    return ele;
  }
  return nullptr;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallQuad4PoroType::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::quad4>>(id, owner);
  return ele;
}

void Discret::Elements::WallQuad4PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ4PORO"];

  defs["QUAD4"] = Input::LineDefinition::Builder(defs_wall["QUAD4"])
                      .add_optional_named_double_vector("POROANISODIR1", 2)
                      .add_optional_named_double_vector("POROANISODIR2", 2)
                      .add_optional_named_double_vector("POROANISONODALCOEFFS1", 4)
                      .add_optional_named_double_vector("POROANISONODALCOEFFS2", 4)
                      .build();
}

int Discret::Elements::WallQuad4PoroType::initialize(Core::FE::Discretization& dis)
{
  Discret::Elements::Wall1Type::initialize(dis);
  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;
    auto* actele = dynamic_cast<Discret::Elements::Wall1Poro<Core::FE::CellType::quad4>*>(
        dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->init_element();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                       |
 *----------------------------------------------------------------------*/
Discret::Elements::WallQuad9PoroType Discret::Elements::WallQuad9PoroType::instance_;

Discret::Elements::WallQuad9PoroType& Discret::Elements::WallQuad9PoroType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::WallQuad9PoroType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::Wall1Poro<Core::FE::CellType::quad9>(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallQuad9PoroType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLQ9PORO")
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::quad9>>(id, owner);
    return ele;
  }
  return nullptr;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallQuad9PoroType::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::quad9>>(id, owner);
  return ele;
}

void Discret::Elements::WallQuad9PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLQ9PORO"];

  defs["QUAD9"] = Input::LineDefinition::Builder(defs_wall["QUAD9"])
                      .add_optional_named_double_vector("POROANISODIR1", 2)
                      .add_optional_named_double_vector("POROANISODIR2", 2)
                      .build();
}

int Discret::Elements::WallQuad9PoroType::initialize(Core::FE::Discretization& dis)
{
  Discret::Elements::Wall1Type::initialize(dis);
  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;
    auto* actele = dynamic_cast<Discret::Elements::Wall1Poro<Core::FE::CellType::quad9>*>(
        dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->init_element();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  NURBS 4 Element                                       |
 *----------------------------------------------------------------------*/

Discret::Elements::WallNurbs4PoroType Discret::Elements::WallNurbs4PoroType::instance_;

Discret::Elements::WallNurbs4PoroType& Discret::Elements::WallNurbs4PoroType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::WallNurbs4PoroType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs4>(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallNurbs4PoroType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN4PORO")
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs4>>(id, owner);
    return ele;
  }
  return nullptr;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallNurbs4PoroType::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs4>>(id, owner);
  return ele;
}

void Discret::Elements::WallNurbs4PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLN4PORO"];

  defs["NURBS4"] = Input::LineDefinition::Builder(defs_wall["NURBS4"])
                       .add_optional_named_double_vector("POROANISODIR1", 2)
                       .add_optional_named_double_vector("POROANISODIR2", 2)
                       .build();
}

int Discret::Elements::WallNurbs4PoroType::initialize(Core::FE::Discretization& dis)
{
  Discret::Elements::Wall1Type::initialize(dis);
  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;
    auto* actele = dynamic_cast<Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs4>*>(
        dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->init_element();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  NURBS 9 Element                                       |
 *----------------------------------------------------------------------*/

Discret::Elements::WallNurbs9PoroType Discret::Elements::WallNurbs9PoroType::instance_;

Discret::Elements::WallNurbs9PoroType& Discret::Elements::WallNurbs9PoroType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::WallNurbs9PoroType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs9>(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallNurbs9PoroType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLN9PORO")
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs9>>(id, owner);
    return ele;
  }
  return nullptr;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallNurbs9PoroType::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs9>>(id, owner);
  return ele;
}

void Discret::Elements::WallNurbs9PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLN9PORO"];

  defs["NURBS9"] = Input::LineDefinition::Builder(defs_wall["NURBS9"])
                       .add_optional_named_double_vector("POROANISODIR1", 2)
                       .add_optional_named_double_vector("POROANISODIR2", 2)
                       .build();
}

int Discret::Elements::WallNurbs9PoroType::initialize(Core::FE::Discretization& dis)
{
  Discret::Elements::Wall1Type::initialize(dis);
  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;
    auto* actele = dynamic_cast<Discret::Elements::Wall1Poro<Core::FE::CellType::nurbs9>*>(
        dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->init_element();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                       |
 *----------------------------------------------------------------------*/

Discret::Elements::WallTri3PoroType Discret::Elements::WallTri3PoroType::instance_;

Discret::Elements::WallTri3PoroType& Discret::Elements::WallTri3PoroType::instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::Elements::WallTri3PoroType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::Wall1Poro<Core::FE::CellType::tri3>(-1, -1);
  object->unpack(buffer);
  return object;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallTri3PoroType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALLT3PORO")
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::tri3>>(id, owner);
    return ele;
  }
  return nullptr;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::WallTri3PoroType::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::Wall1Poro<Core::FE::CellType::tri3>>(id, owner);
  return ele;
}

void Discret::Elements::WallTri3PoroType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, std::map<std::string, Input::LineDefinition>> definitions_wall;
  Wall1Type::setup_element_definition(definitions_wall);

  std::map<std::string, Input::LineDefinition>& defs_wall = definitions_wall["WALL"];

  std::map<std::string, Input::LineDefinition>& defs = definitions["WALLT3PORO"];

  defs["TRI3"] = Input::LineDefinition::Builder(defs_wall["TRI3"])
                     .add_optional_named_double_vector("POROANISODIR1", 2)
                     .add_optional_named_double_vector("POROANISODIR2", 2)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS1", 3)
                     .add_optional_named_double_vector("POROANISONODALCOEFFS2", 3)
                     .build();
}

int Discret::Elements::WallTri3PoroType::initialize(Core::FE::Discretization& dis)
{
  Discret::Elements::Wall1Type::initialize(dis);
  for (int i = 0; i < dis.num_my_col_elements(); ++i)
  {
    if (dis.l_col_element(i)->element_type() != *this) continue;
    auto* actele =
        dynamic_cast<Discret::Elements::Wall1Poro<Core::FE::CellType::tri3>*>(dis.l_col_element(i));
    if (!actele) FOUR_C_THROW("cast to Wall1_Poro* failed");
    actele->init_element();
  }
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
