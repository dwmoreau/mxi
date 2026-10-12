import numpy as np
import pytest

import mxi_stills


def lattice_rlp():
    # A primitive orthorhombic lattice, 10 x 20 x 30 A: its reciprocal lattice
    # points with indices up to 4.
    a_star = np.diag([1 / 10.0, 1 / 20.0, 1 / 30.0])
    hkl = [(h, k, l) for h in range(-4, 5) for k in range(-4, 5) for l in range(-4, 5)
           if (h, k, l) != (0, 0, 0)]
    return np.array([a_star @ np.array(m, float) for m in hkl])


def test_fft1d_finds_a_basis_of_the_lattice():
    # The strongest vectors include three that span the 6000 A^3 cell (b may
    # come as a + b, since the grid need not hold b's direction itself).
    v = mxi_stills.fft1d_basis_vectors(lattice_rlp(), 60.0)
    assert v.shape[1] == 3 and 3 <= len(v) <= 30
    assert np.isclose(np.linalg.norm(v[0]), 10, rtol=0.02)
    volumes = [abs(np.linalg.det(v[[i, j, k]])) for i in range(6) for j in range(i + 1, 6)
               for k in range(j + 1, 6)]
    assert any(np.isclose(x, 6000, rtol=0.05) for x in volumes)


def test_fft1d_trace():
    v, trace = mxi_stills.fft1d_basis_vectors(lattice_rlp(), 60.0, trace=True)
    assert not trace["adaptive"]
    assert trace["n_directions"] == len(mxi_stills.hemisphere_grid(trace["incr"]))
    assert trace["refined"].shape == (len(trace["unrefined"]), 5)


def test_fft1d_refuses_seven_spots():
    with pytest.raises(mxi_stills.TooFewSpots):
        mxi_stills.fft1d_basis_vectors(lattice_rlp()[:7], 60.0)
