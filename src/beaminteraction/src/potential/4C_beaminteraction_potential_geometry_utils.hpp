// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_POTENTIAL_GEOMETRY_UTILS_HPP
#define FOUR_C_BEAMINTERACTION_POTENTIAL_GEOMETRY_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN


namespace BeamInteraction::Potential
{
  namespace Geo
  {
    // point-to-curve projection: solve minimal distance problem
    // convergence criteria for local Newton's method
    const unsigned int POINT_TO_CURVE_PROJECTION_MAX_NUM_ITER = 50;
    const double POINT_TO_CURVE_PROJECTION_TOLERANCE_RESIDUUM = 1.0e-10;
    const double POINT_TO_CURVE_PROJECTION_TOLERANCE_INCREMENT = 1.0e-10;
    // threshold values for sanity checks
    const double POINT_TO_CURVE_PROJECTION_IDENTICAL_POINTS_TOLERANCE = 1.0e-12;
    const double POINT_TO_CURVE_PROJECTION_NONUNIQUE_MINIMAL_DISTANCE_TOLERANCE = 1.0e-12;

    /** \brief solves minimal distance problem to find the closest point on a 3D spatial curve
     *         (i.e. its curve parameter value) relative to a given point
     *         a.k.a 'unilateral' closest-point projection
     *
     */
    template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
    bool point_to_curve_projection(Core::LinAlg::Matrix<3, 1, T> const& r_source, T& xi_target,
        double const& xi_target_initial_guess,
        const Core::LinAlg::Matrix<3 * numnodes * numnodalvalues, 1, T>&
            target_centerline_dof_values,
        const Core::FE::CellType& target_distype, double target_ele_ref_length);

    /** \brief evaluates residual of orthogonality condition for so-called unilateral closest-point
     *         projection, i.e. a point-to-curve projection
     *
     */
    template <typename T>
    void evaluate_point_to_curve_orthogonality_condition(T& f,
        const Core::LinAlg::Matrix<3, 1, T>& delta_r, const double norm_delta_r,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target);

    /** \brief evaluates Jacobian of orthogonality condition for so-called unilateral closest-point
     *         projection, i.e. a point-to-curve projection
     *
     */
    template <typename T>
    bool evaluate_linearization_point_to_curve_orthogonality_condition(T& df,
        const Core::LinAlg::Matrix<3, 1, T>& delta_r, const double norm_delta_r,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& r_xixi_target);

    /** \brief compute linearization of parameter coordinate on target if determined by a
     *         point-to-curve projection
     *
     */
    template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
    void calc_linearization_point_to_curve_projection_parameter_coord_target(
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, T>& lin_xi_target_sourceDofs,
        Core::LinAlg::Matrix<1, 3 * numnodes * numnodalvalues, T>& lin_xi_target_targetDofs,
        const Core::LinAlg::Matrix<3, 1, T>& delta_r,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& r_xixi_target,
        const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, double>& N_source,
        const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, T>& N_target,
        const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, T>& N_xi_target);

    /** \brief point-to-curve projection:
     *         partial derivatives of the parameter coordinate on target xi_target with respect to
     *         centerline position of source point, target point and centerline tangent of target
     *
     */
    template <typename T>
    void calc_point_to_curve_projection_parameter_coord_target_partial_derivs(
        Core::LinAlg::Matrix<1, 3, T>& xi_target_partial_r_source,
        Core::LinAlg::Matrix<1, 3, T>& xi_target_partial_r_target,
        Core::LinAlg::Matrix<1, 3, T>& xi_target_partial_r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& delta_r,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& r_xixi_target);

    /** \brief point-to-curve projection:
     *         partial second derivatives of the parameter coordinate on target xi_target with
     *         respect to centerline position of source point, target point and centerline tangent
     * of target
     *
     */
    template <typename T>
    void calc_point_to_curve_projection_parameter_coord_target_partial2nd_derivs(
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_source_partial_r_source,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_source_partial_r_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_source_partial_r_xi_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_source_partial_r_xixi_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_target_partial_r_source,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_target_partial_r_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_target_partial_r_xi_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_target_partial_r_xixi_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_xi_target_partial_r_source,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_xi_target_partial_r_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_xi_target_partial_r_xi_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_xi_target_partial_r_xixi_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_xixi_target_partial_r_source,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_xixi_target_partial_r_target,
        Core::LinAlg::Matrix<3, 3, T>& xi_target_partial_r_xixi_target_partial_r_xi_target,
        const Core::LinAlg::Matrix<1, 3, T>& xi_target_partial_r_source,
        const Core::LinAlg::Matrix<1, 3, T>& xi_target_partial_r_target,
        const Core::LinAlg::Matrix<1, 3, T>& xi_target_partial_r_xi_target,
        const Core::LinAlg::Matrix<3, 3, T>& delta_r_deriv_r_source,
        const Core::LinAlg::Matrix<3, 3, T>& delta_r_deriv_r_target,
        const Core::LinAlg::Matrix<3, 3, T>& delta_r_deriv_r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& delta_r,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& r_xixi_target,
        const Core::LinAlg::Matrix<3, 1, T>& r_xixixi_target);

    /** \brief point-to-curve projection:
     *         partial derivative of the orthogonality condition with respect to parameter
     * coordinate on target xi_target
     *
     */
    template <typename T>
    void calc_ptc_projection_orthogonality_condition_partial_deriv_parameter_coord_target(
        T& orthogon_condition_partial_xi_target, const Core::LinAlg::Matrix<3, 1, T>& delta_r,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& r_xixi_target);

    /** \brief point-to-curve projection:
     *         partial derivative of the orthogonality condition with respect to centerline position
     *         on source
     *
     */
    template <typename T>
    void calc_ptc_projection_orthogonality_condition_partial_deriv_cl_pos_source(
        Core::LinAlg::Matrix<1, 3, T>& orthogon_condition_partial_r_source,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target);

    /** \brief point-to-curve projection:
     *         partial derivative of the orthogonality condition with respect to centerline position
     *         on target
     *
     */
    template <typename T>
    void calc_ptc_projection_orthogonality_condition_partial_deriv_cl_pos_target(
        Core::LinAlg::Matrix<1, 3, T>& orthogon_condition_partial_r_target,
        const Core::LinAlg::Matrix<3, 1, T>& r_xi_target);

    /** \brief point-to-curve projection:
     *         partial derivative of the orthogonality condition with respect to centerline tangent
     *         on target
     *
     */
    template <typename T>
    void calc_ptc_projection_orthogonality_condition_partial_deriv_cl_tangent_target(
        Core::LinAlg::Matrix<1, 3, T>& orthogon_condition_partial_r_xi_target,
        const Core::LinAlg::Matrix<3, 1, T>& delta_r);

    /** \brief calculate angle enclosed by two vectors a and b
     *
     */
    template <typename T>
    void calc_enclosed_angle(T& angle, T& cosine_angle, const Core::LinAlg::Matrix<3, 1, T>& a,
        const Core::LinAlg::Matrix<3, 1, T>& b);

  }  // namespace Geo
}  // namespace BeamInteraction::Potential

FOUR_C_NAMESPACE_CLOSE

#endif
