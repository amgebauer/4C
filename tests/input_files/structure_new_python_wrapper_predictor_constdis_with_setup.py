# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import numpy as np
from pathlib import Path

# global flag to check whether setup has been called
_setup_called = False

# global step counter to see how often compute() has been called
step_counter = 0


def setup(context):
    """
    Perform optional once-per-simulation initialization for the PythonWrapper predictor.

    Parameters
    ----------
    context : PythonContext
        Read-only object with static problem information assembled by 4C.
        The object currently contains:
        - input_file: path to the top-level 4C input file
        - problem_dim: spatial dimension of the problem
        - mesh: lightweight mesh representation gathered on rank 0
        - max_time: the final time in the simulation
        - max_step: the maximum number of time/load steps

        The mesh attribute is of type PythonMesh and contains:
        - node_ids: global node ids
        - coordinates: nodal coordinates
        - disp_dof_ids: displacement DOF ids associated with each node
        - element_ids: global element ids
        - connectivity: element-to-node connectivity

    Notes
    -----
    - This function is optional. If it is not provided, 4C will skip the setup step.
    - The function is called only once per simulation and only on MPI rank 0.
    - The provided data should be treated as read-only.
    - A return value is ignored.
    """
    # mark the setup function as called
    global _setup_called
    _setup_called = True
    # check the provided values
    assert (
        context.input_file
        == (
            Path(__file__).parent
            / "structure_new_python_wrapper_predictor_constdis_with_setup.4C.yaml"
        ).as_posix()
    )
    assert context.problem_dim == 3
    assert context.step_max == 2
    assert context.time_max == 1
    # check the mesh data
    assert context.mesh.node_ids.shape == (24,)
    assert context.mesh.coordinates.shape == (24, 3)
    assert context.mesh.disp_dof_ids.shape == (24, 3)
    assert context.mesh.element_ids.shape == (6,)
    assert len(context.mesh.connectivity) == 6
    for element in context.mesh.connectivity:
        assert len(element) == 8
    # check the node ids
    assert np.array_equal(context.mesh.node_ids, np.arange(0, 24))
    # check the coordinates
    expected_coordinates = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.0, 0.5, 1.0],
            [0.6666666666666666, 0.0, 1.0],
            [0.6666666666666667, 0.0, 0.5],
            [0.6666666666666667, 0.5, 0.5],
            [0.6666666666666667, 0.5, 1.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.6666666666666667, 0.0, 0.0],
            [0.6666666666666666, 0.5, 0.0],
            [1.3333333333333333, 0.0, 1.0],
            [1.3333333333333335, 0.0, 0.5],
            [1.3333333333333335, 0.5, 0.5],
            [1.3333333333333335, 0.5, 1.0],
            [1.3333333333333335, 0.0, 0.0],
            [1.3333333333333333, 0.5, 0.0],
            [2.0, 0.0, 1.0],
            [2.0, 0.0, 0.5],
            [2.0, 0.5, 0.5],
            [2.0, 0.5, 1.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.5, 0.0],
        ]
    )
    assert np.array_equal(context.mesh.coordinates, expected_coordinates)
    # check the displacement dof ids
    assert np.array_equal(context.mesh.disp_dof_ids, np.arange(0, 72).reshape(24, 3))
    # check the element ids
    assert np.array_equal(context.mesh.element_ids, np.arange(0, 6))
    # check the connectivity array
    expected_connectivity = np.array(
        [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 8, 9, 2, 5, 10, 11, 6],
            [4, 5, 6, 7, 12, 13, 14, 15],
            [5, 10, 11, 6, 13, 16, 17, 14],
            [12, 13, 14, 15, 18, 19, 20, 21],
            [13, 16, 17, 14, 19, 22, 23, 20],
        ]
    )
    assert np.array_equal(context.mesh.connectivity, expected_connectivity)


# simulate a ConstDis predictor
def compute(state):
    """
    Perform one predictor step for the PythonWrapper predictor.

    Parameters
    ----------
    state : PythonState
        Object with predictor data for the current step, assembled by 4C on rank 0.

        Scalar entries:
        - step: current load/time step index
        - time: current time value
        - dt: current time/load increment
        - num_global_dofs: total number of global DOFs

        Read-only array entries:
        - dis_n: displacement state at step n
        - vel_n: velocity state at step n
        - acc_n: acceleration state at step n

        Writable array entries:
        - dis_np: displacement state at step n+1
        - vel_np: velocity state at step n+1
        - acc_np: acceleration state at step n+1

    Contract
    --------
    - The arrays dis_np, vel_np, and acc_np must be updated in place.
    - Do not replace these dictionary entries with new arrays.
    - Supported patterns are for example:
          np.copyto(state.dis_np, ...)
          state.dis_np[:] = ...
          state.vel_np.fill(0.0)
    - Unsupported pattern:
          state.dis_np = new_array

    Notes
    -----
    - This function is called once per predictor step and only on MPI rank 0.
    - The arrays dis_n, vel_n, and acc_n must be treated as read-only.
    - 4C uses the modified contents of dis_np, vel_np, and acc_np after this function returns.
    - The function must not return a value; the predictor output is taken from the in-place modified arrays.
    """
    # check that the setup function has been called
    assert _setup_called == True
    # increment the global step counter
    global step_counter
    step_counter += 1

    # Make sure that the arrays have the correct properties
    assert not state.dis_n.flags.writeable
    assert not state.vel_n.flags.writeable
    assert not state.acc_n.flags.writeable
    assert state.dis_np.flags.writeable
    assert state.vel_np.flags.writeable
    assert state.acc_np.flags.writeable
    # check the provided values
    assert state.step == step_counter
    assert state.time == step_counter * 0.5
    assert state.dt == 0.5
    assert state.num_global_dofs == 72
    # check the correct shape of the provided vectors
    assert state.dis_n.shape == (72,)
    assert state.dis_np.shape == (72,)
    assert state.vel_n.shape == (72,)
    assert state.vel_np.shape == (72,)
    assert state.acc_n.shape == (72,)
    assert state.acc_np.shape == (72,)
    # make a crazy prediction that can be easily tested against in 4C
    # IMPORTANT:
    # Modify the arrays contained in the state dictionary in-place!
    state.dis_np[:] = np.array(state.dis_n, copy=True) + np.tile([1, -2, 3], 24)
    state.vel_np[:] = np.zeros_like(state.vel_np)
    state.acc_np[:] = np.zeros_like(state.acc_np)
