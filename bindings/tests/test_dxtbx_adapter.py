"""The dxtbx adapter, on stand-in models: no dxtbx is needed to run these."""

import numpy as np
import pytest

import mxi_stills
from mxi_stills.dxtbx_adapter import Workspace, detector_arrays, reflection_arrays


class SimplePxMmStrategy:
    pass


class ParallaxCorrectedPxMmStrategy:
    def __init__(self, mu, t0):
        self._mu, self._t0 = mu, t0

    def mu(self):
        return self._mu

    def t0(self):
        return self._t0


class OffsetPxMmStrategy:
    pass


class Panel:
    def __init__(self, origin, strategy):
        self.origin, self.strategy = origin, strategy

    def get_origin(self):
        return self.origin

    def get_fast_axis(self):
        return (1.0, 0.0, 0.0)

    def get_slow_axis(self):
        return (0.0, -1.0, 0.0)

    def get_pixel_size(self):
        return (0.075, 0.075)

    def get_px_mm_strategy(self):
        return self.strategy


class Beam:
    def get_s0(self):
        return (0.0, 0.0, -1.0)

    def get_wavelength(self):
        return 1.0


class Experiment:
    def __init__(self, detector):
        self.detector, self.beam = detector, Beam()


def two_panels(z=-100.0):
    return [
        Panel((-40.0, 40.0, z), ParallaxCorrectedPxMmStrategy(3.9, 0.45)),
        Panel((0.0, 40.0, z), SimplePxMmStrategy()),
    ]


def spots():
    return {
        "panel": np.array([0, 1, 1]),
        "xyzobs.px.value": np.array(
            [[10.0, 20.0, 0.5], [300.0, 400.0, 0.5], [5.0, 5.0, 0.5]]
        ),
        "xyzobs.px.variance": np.full((3, 3), 0.25),
    }


def test_parallax_comes_from_the_strategy():
    a = detector_arrays(two_panels())
    assert a["mu"].tolist() == [3.9, 0.0]
    assert a["t0"].tolist() == [0.45, 0.0]


def test_an_unknown_strategy_is_refused():
    detector = two_panels()
    detector[1].strategy = OffsetPxMmStrategy()
    with pytest.raises(ValueError, match="OffsetPxMmStrategy"):
        detector_arrays(detector)


def test_missing_variance_is_one():
    table = spots()
    del table["xyzobs.px.variance"]
    _, _, var = reflection_arrays(table)
    assert (var == 1.0).all()


def test_the_workspace_rebuilds_only_on_change():
    ws = Workspace()
    first = ws.map_spots(Experiment(two_panels()), spots())
    ws.map_spots(Experiment(two_panels()), spots())
    assert ws.n_builds == 1
    moved = ws.map_spots(Experiment(two_panels(z=-110.0)), spots())
    assert ws.n_builds == 2
    assert not np.array_equal(first["rlp"], moved["rlp"])


def test_mapping_shapes_and_the_stub():
    ws = Workspace()
    out = ws.map_spots(Experiment(two_panels()), spots())
    for key in ("xyzobs.mm.value", "xyzobs.mm.variance", "s1", "rlp"):
        assert out[key].shape == (3, 3)
    assert (out["xyzobs.mm.value"][:, 2] == 0).all()
    # The simple panel is a plain multiplication by the pixel size.
    assert out["xyzobs.mm.value"][1, 0] == 300.0 * 0.075
    result = ws.index(Experiment(two_panels()), spots())
    assert not result.indexed and result.reason and result.n_spots == 3


def test_a_panel_past_the_end_is_refused():
    table = spots()
    table["panel"] = np.array([0, 1, 2])
    with pytest.raises(mxi_stills.GeometryError):
        Workspace().map_spots(Experiment(two_panels()), table)
