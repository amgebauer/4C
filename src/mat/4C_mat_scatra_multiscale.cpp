// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_scatra_multiscale.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::PAR::ScatraMultiScale::ScatraMultiScale(const Core::Mat::PAR::Parameter::Data& matdata)
    : ScatraMat(matdata),
      ScatraMicroMacroCoupling(matdata),
      porosity_(matdata.parameters.get<double>("POROSITY")),
      tortuosity_(matdata.parameters.get<double>("TORTUOSITY"))
{
  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMultiScale::create_material()
{
  return Teuchos::make_rcp<Mat::ScatraMultiScale>(this);
}


Mat::ScatraMultiScaleType Mat::ScatraMultiScaleType::instance_;


Core::Communication::ParObject* Mat::ScatraMultiScaleType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::ScatraMultiScale* ScatraMatMultiScale = new Mat::ScatraMultiScale();
  ScatraMatMultiScale->unpack(buffer);
  return ScatraMatMultiScale;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::ScatraMultiScale::ScatraMultiScale() : params_(nullptr) { return; }


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
Mat::ScatraMultiScale::ScatraMultiScale(Mat::PAR::ScatraMultiScale* params)
    : ScatraMat(params), params_(params)
{
  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::ScatraMultiScale::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // pack base class material
  ScatraMat::pack(data);

  return;
}


/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::ScatraMultiScale::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::ScatraMultiScale*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not match calling type %d!", mat->type(),
            material_type());
    }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  ScatraMat::unpack(basedata_buffer);

  // final safety check
  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");

  return;
}

FOUR_C_NAMESPACE_CLOSE
