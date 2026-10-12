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
    assert all(trace[k] >= 0 for k in ("seconds_quick", "seconds_scan", "seconds_refine"))


def test_fft1d_grid_settings():
    # rstbx's grid is the default. The data would allow 0.07 rad here, so a
    # coarser max_grid searches fewer directions and still finds the 10 A axis.
    rlp = lattice_rlp()
    default = mxi_stills.fft1d_basis_vectors(rlp, 60.0)
    same = mxi_stills.fft1d_basis_vectors(rlp, 60.0, max_grid=0.029, quick_grid=0.016,
                                          refine_divisor=100.0)
    assert np.array_equal(default, same)
    v, trace = mxi_stills.fft1d_basis_vectors(rlp, 60.0, trace=True)
    assert trace["characteristic_grid"] == trace["incr"] == 0.029
    v, trace = mxi_stills.fft1d_basis_vectors(rlp, 60.0, trace=True, max_grid=0.05)
    assert trace["incr"] == 0.05
    assert trace["n_directions"] == len(mxi_stills.hemisphere_grid(0.05))
    assert np.isclose(np.linalg.norm(v[0]), 10, rtol=0.02)


def test_fft1d_refuses_seven_spots():
    with pytest.raises(mxi_stills.TooFewSpots):
        mxi_stills.fft1d_basis_vectors(lattice_rlp()[:7], 60.0)


def test_real_space_grid_search_finds_the_cell_axes():
    # The target cell's own lengths, searched against its own reciprocal
    # lattice: the strongest vectors include each axis at its length.
    s = mxi_stills.RealSpaceGridSearch((10.0, 20.0, 30.0, 90, 90, 90))
    assert len(s.search_vectors) == 3 * len(mxi_stills.hemisphere_grid(0.02))
    v = s.find_basis_vectors(lattice_rlp())
    assert 3 <= len(v) <= 30
    lengths = np.linalg.norm(v, axis=1)
    for axis, length in ((0, 10), (1, 20), (2, 30)):
        along = [x for x, n in zip(v, lengths)
                 if np.isclose(n, length) and abs(abs(x[axis]) - n) < 0.05 * n]
        assert along, (axis, length)
