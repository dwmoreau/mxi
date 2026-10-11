import numpy as np
import pytest

import mxi_stills

S0 = np.array([0.0, 0.0, -1.0])  # 1 A


def detector():
    return mxi_stills.Detector(
        origin=np.array([[-50.0, 50.0, -100.0]]),
        fast=np.array([[1.0, 0.0, 0.0]]),
        slow=np.array([[0.0, -1.0, 0.0]]),
        pixel_size=np.array([[0.1, 0.1]]),
        mu=np.zeros(1),
        t0=np.zeros(1),
    )


def rotation(axis, angle):
    u = np.asarray(axis, float) / np.linalg.norm(axis)
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def observed():
    """A cubic cell's reflections near the Ewald sphere, observed where the crystal
    predicts them."""
    U = rotation([1, 2, 3], 0.7)
    B = np.diag([1 / 60.0] * 3)
    h = np.array(
        [(a, b, c) for a in range(-12, 13) for b in range(-12, 13) for c in range(-12, 13)]
    )
    q = h @ (U @ B).T
    s0_len = 1.0
    excess = np.abs(np.linalg.norm(q + S0, axis=1) - s0_len)
    h = h[(excess < 2e-3) & (np.abs(h).sum(axis=1) > 0)]
    p = mxi_stills.predict_still(detector(), S0, U, B, h, np.zeros(len(h), np.int64))
    keep = (p["predicted"] == 1) & (np.abs(p["xyzcal.mm"][:, :2]).max(axis=1) < 50)
    return U, B, h[keep], p["xyzcal.mm"][keep]


def refine(U, B, h, xyz, **kw):
    n = len(h)
    return mxi_stills.refine_still(
        detector(),
        S0,
        U,
        B,
        " P 2 2 3",
        h,
        np.zeros(n, np.int64),
        xyz,
        np.tile([1e-4, 1e-4, 0.0], (n, 1)),
        np.zeros(n, np.uint64),
        **kw,
    )


def test_prediction_and_linearisation_shapes():
    U, B, h, xyz = observed()
    assert len(h) > 20
    lin = mxi_stills.linearise_still(
        detector(), S0, U, B, " P 2 2 3", h, np.zeros(len(h), np.int64), xyz,
        np.tile([1e-4, 1e-4, 0.0], (len(h), 1)), np.zeros(len(h), np.uint64),
    )
    assert lin["parameter_names"] == ["Crystal1Phi1", "Crystal1Phi2", "Crystal1Phi3",
                                      "Crystal1g_param_0"]
    assert lin["jacobian"].shape == (3 * len(h), 4)
    assert np.allclose(lin["residual"][:, :2], 0.0, atol=1e-12)


@pytest.mark.parametrize("mode", ["dials", "dials_limit", "analytic"])
def test_refinement_moves_towards_the_crystal_that_made_the_spots(mode):
    # The Delta-Psi restraint pulls the minimum slightly off that crystal, since its
    # spots are not exactly on the Ewald sphere.
    U, B, h, xyz = observed()
    start = rotation([0.2, 0.9, 0.1], np.radians(0.05)) @ U
    out = refine(start, B * 1.001, h, xyz, axis_derivative=mode)
    assert not out["failed"], out["message"]
    assert out["n_used"] == len(h)
    start_rmsd = out["journal"][0]["rmsd"]
    assert out["rmsd"][0] < 0.1 * start_rmsd[0] and out["rmsd"][1] < 0.1 * start_rmsd[1]
    assert out["journal"][-1]["objective"] < out["journal"][0]["objective"]


def test_an_unknown_axis_derivative_is_refused():
    U, B, h, xyz = observed()
    with pytest.raises(ValueError):
        refine(U, B, h, xyz, axis_derivative="finite")
