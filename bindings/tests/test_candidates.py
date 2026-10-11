import numpy as np
import pytest

import mxi_stills


def test_niggli_reduction_orders_the_axes():
    r_inv, n = mxi_stills.niggli_reduction(np.array([20.0, 10.0, 30.0, 90.0, 90.0, 90.0]))
    assert n >= 1
    # The new a is the old b.
    assert np.array_equal(r_inv[:, 0], [0, -1, 0])


def test_niggli_reduction_iteration_limit():
    with pytest.raises(mxi_stills.NiggliLimitExceeded):
        mxi_stills.niggli_reduction(np.array([20.0, 10.0, 30.0, 90.0, 90.0, 90.0]), iteration_limit=0)


def test_candidates_of_an_orthogonal_basis():
    bv = np.array([[10.0, 0, 0], [0, 20.0, 0], [0, 0, 30.0]])
    A, cell = mxi_stills.candidate_orientation_matrices(bv)
    assert A.shape == (1, 3, 3)
    assert np.allclose(cell[0], [10, 20, 30, 90, 90, 90])
    assert np.isclose(abs(np.linalg.det(A[0])), 1 / 6000)


def test_assign_indices_global():
    A = np.diag([0.1, 0.1, 0.1])[None]
    rlp = np.array([[0.1, 0.2, 0.0], [0.15, 0.0, 0.0], [0.5, 0.0, 0.0]])
    miller, crystal = mxi_stills.assign_indices_global(rlp, np.full(3, -1), A, d_min=3.0)
    assert miller.tolist() == [[1, 2, 0], [0, 0, 0], [0, 0, 0]]
    assert crystal.tolist() == [0, -1, -1]
