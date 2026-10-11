# Stills indexing: where it departs from dials

`mxi_stills` reimplements, in new C++, the algorithms `dials.stills_process` uses to index a
still. Each entry here is a place where it does something different from dials, or from the
paper the algorithm comes from. An entry is written when the change is made, and gives its
reason: a measured cost, with the run that measured it, or a stated weakness of the model.

| # | date | what differs | from (dials code or paper) | reason | measured effect |
|---|---|---|---|---|---|

## What follows dials exactly

**Spot centroids to reciprocal space** (`src/stills/geometry.cc`). This is dials'
`centroid_px_to_mm_panel` and `map_centroids_to_reciprocal_space` for a still, with
dxtbx's `parallax_correction_inv` and `get_lab_coord`. The arithmetic is the same and
in the same order, and it is compiled without fused multiply-add, as dials is. On all
26,639 strong spots of the 460-frame corpus (detectors of 256, 32, 8 and 1 panels),
every millimetre position, variance, s1 and reciprocal lattice point is bit-identical to
dials. Evidence: `stills_process_profiling/claude/rewrite/runs/r5-plumbing/`.

Not yet departures, but named as candidates: replacing Nave's mosaicity with a stills
partiality model, and replacing Niggli reduction with Selling reduction if a profile shows that
reduction costs time.
