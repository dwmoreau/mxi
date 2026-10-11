import numpy as np

import mxi_stills


def panel(mu=0.0, t0=0.0):
    return mxi_stills.Detector(
        origin=np.array([[-50.0, 50.0, -100.0]]),
        fast=np.array([[1.0, 0.0, 0.0]]),
        slow=np.array([[0.0, -1.0, 0.0]]),
        pixel_size=np.array([[0.1, 0.1]]),
        mu=np.array([mu]),
        t0=np.array([t0]),
    )


def test_mm_to_px_without_parallax_is_the_pixel_size():
    xy = np.array([[12.3, 45.6, 0.0], [0.05, 99.0, 0.0]])
    px = mxi_stills.mm_to_px(panel(), np.zeros(2, np.int64), xy)
    assert np.array_equal(px, xy[:, :2] / 0.1)


def test_mm_to_px_with_parallax_inverts_px_to_mm():
    det = panel(mu=3.9, t0=0.32)
    px = np.array([[120.0, 340.0, 0.0], [900.0, 20.0, 0.0]])
    mm = mxi_stills.map_spots(
        det, np.array([0.0, 0.0, -1.0]), 1.0, np.zeros(2, np.int64), px, np.ones((2, 3))
    )["xyzobs.mm.value"]
    back = mxi_stills.mm_to_px(det, np.zeros(2, np.int64), mm)
    # dxtbx's inverse correction is an approximation, good to well under a pixel.
    assert np.abs(back - px[:, :2]).max() < 1e-2
    assert np.abs(back - px[:, :2]).max() > 0


def test_sauter_poon_rejects_a_far_spot():
    n = 60
    p = (np.arange(n) + 0.5) / n
    r = 0.05 * np.sqrt(-2 * np.log(1 - p))
    r[-1] = 1.0
    obs = np.column_stack([r, np.zeros(n)])
    out = mxi_stills.sauter_poon(obs, np.zeros((n, 2)))
    assert out["status"] == "ok"
    assert not out["kept"][-1] and out["kept"][: n // 2].all()
    assert abs(out["sigma"] - 0.05) < 0.005


def test_sauter_poon_needs_ten_spots():
    obs = np.column_stack([np.linspace(0.01, 0.09, 9), np.zeros(9)])
    assert mxi_stills.sauter_poon(obs, np.zeros((9, 2)))["status"] == "too_few"


def test_nave_on_a_band():
    rng = np.random.default_rng(0)
    d = rng.uniform(2.0, 10.0, 200)
    alpha, eta = 2 / 3000.0, 2e-3
    psi = rng.uniform(-1, 1, 200) * (d * alpha + eta) / 2
    two_theta = mxi_stills.two_theta_deg(d, 1.0)
    out = mxi_stills.nave_parameters(psi, d, two_theta, 1e6, 1.0)
    assert out["status"] in ("ok", "ewald_volume")
    assert out["accepted"].mean() > 0.9
    assert out["domain_size_ang"] > 0 and out["half_mosaicity_deg"] > 0
