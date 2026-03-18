# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import numpy as np

# global step counter to see how often compute() has been called
step_counter = 0


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
    # try to modify the read-only arrays should raise an error
    try:
        state.dis_n[0] = 42
        raise AssertionError("dis_n should be read-only")
    except ValueError:
        # This should raise a ValueError because dis_n is read-only
        pass
    try:
        state.vel_n[0] = 42
        raise AssertionError("vel_n should be read-only")
    except ValueError:
        # This should raise a ValueError because vel_n is read-only
        pass
    try:
        state.acc_n[0] = 42
        raise AssertionError("acc_n should be read-only")
    except ValueError:
        # This should raise a ValueError because acc_n is read-only
        pass
    # make the same prediction as a ConstDis predictor would
    np.copyto(state.dis_np, np.array(state.dis_n, copy=True) + np.tile([1, -2, 3], 24))
    state.vel_np.fill(0)
    state.acc_np.fill(0)
