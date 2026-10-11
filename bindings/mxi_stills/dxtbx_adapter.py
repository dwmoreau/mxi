"""Plain arrays from live dxtbx models and dials reflection tables.

Nothing here imports dxtbx or dials: the functions call the models' own
methods, so the package stays importable without either.

Panel geometry is read from the live panel in the laboratory frame, with any
hierarchy composed (``get_origin``, ``get_fast_axis``, ``get_slow_axis``).
The parallax coefficients come from the panel's px/mm strategy, which is
where dxtbx keeps the values it uses; a panel's own ``get_mu`` can differ.
"""

from __future__ import annotations

import numpy as np

from ._mxi_stills import Detector, index_still, map_spots

SUPPORTED_STRATEGIES = ("SimplePxMmStrategy", "ParallaxCorrectedPxMmStrategy")


def detector_arrays(detector):
    """The arrays ``Detector`` is built from, one row per panel."""
    n = len(detector)
    out = {
        "origin": np.empty((n, 3)),
        "fast": np.empty((n, 3)),
        "slow": np.empty((n, 3)),
        "pixel_size": np.empty((n, 2)),
        "mu": np.zeros(n),
        "t0": np.zeros(n),
    }
    for i, panel in enumerate(detector):
        out["origin"][i] = panel.get_origin()
        out["fast"][i] = panel.get_fast_axis()
        out["slow"][i] = panel.get_slow_axis()
        out["pixel_size"][i] = panel.get_pixel_size()
        strategy = panel.get_px_mm_strategy()
        name = type(strategy).__name__
        if name not in SUPPORTED_STRATEGIES:
            raise ValueError(
                f"panel {i}: px/mm strategy {name} is not supported; "
                f"mxi_stills handles {', '.join(SUPPORTED_STRATEGIES)}"
            )
        if name == "ParallaxCorrectedPxMmStrategy":
            out["mu"][i] = strategy.mu()
            out["t0"][i] = strategy.t0()
    return out


def reflection_arrays(reflections):
    """panel, xyzobs.px.value and xyzobs.px.variance as numpy arrays.

    A table without a variance column gets (1, 1, 1), as dials assigns.
    """
    panel = np.asarray(reflections["panel"], dtype=np.int64)
    n = len(panel)
    xyz = np.asarray(reflections["xyzobs.px.value"], dtype=np.float64).reshape(n, 3)
    if "xyzobs.px.variance" in reflections:
        var = np.asarray(reflections["xyzobs.px.variance"], dtype=np.float64)
        var = var.reshape(n, 3)
    else:
        var = np.ones((n, 3))
    return panel, xyz, var


class Workspace:
    """The C++ detector, kept between stills.

    It is rebuilt only when the live detector's arrays change, which for a
    process using a reference geometry is never after the first still.
    """

    def __init__(self):
        self._arrays = None
        self.detector = None
        self.n_builds = 0

    def update(self, dxtbx_detector):
        arrays = detector_arrays(dxtbx_detector)
        if self._arrays is None or any(
            not np.array_equal(arrays[k], self._arrays[k]) for k in arrays
        ):
            self.detector = Detector(**arrays)
            self._arrays = arrays
            self.n_builds += 1
        return self.detector

    def map_spots(self, experiment, reflections):
        """dials' xyzobs.mm, s1 and rlp for one still, as numpy arrays."""
        detector = self.update(experiment.detector)
        beam = experiment.beam
        return map_spots(
            detector,
            np.asarray(beam.get_s0()),
            beam.get_wavelength(),
            *reflection_arrays(reflections),
        )

    def index(self, experiment, reflections):
        detector = self.update(experiment.detector)
        beam = experiment.beam
        return index_still(
            detector,
            np.asarray(beam.get_s0()),
            beam.get_wavelength(),
            *reflection_arrays(reflections),
        )
