import numpy as np
import pytest

import mxi_stills

IDENTITY = ((12, 0, 0, 0, 12, 0, 0, 0, 12, 12, 0, 0, 0, 144),) * 2
C_A_B = ((0, 0, 12, 12, 0, 0, 0, 12, 0, 12, 0, 0, 0, 144),
         (0, 12, 0, 0, 0, 12, 12, 0, 0, 12, 0, 0, 0, 144))


def test_reference_setting_is_its_own_reference():
    number, cb_op = mxi_stills.sgtbx_type(" P 2y")
    assert number == 3
    assert cb_op == IDENTITY


def test_type_of_a_permuted_setting():
    # P 1 2 1 with its axes permuted: still number 3, and the change of basis
    # back to the reference setting is the inverse permutation.
    number, (c, c_inv) = mxi_stills.sgtbx_type(" P 2y", C_A_B)
    assert number == 3
    assert (c_inv, c) == C_A_B


def test_tidy_group_order_is_canonical():
    a = mxi_stills.sgtbx_group(" P 4 2", tidy=True)
    b = mxi_stills.sgtbx_group(" P 4 2", C_A_B)
    b_tidy = mxi_stills.sgtbx_group(" P 4 2", C_A_B, tidy=True)
    assert len(a[3]) == len(b[3]) == 8
    assert b_tidy != a  # a different setting, so different operations
    assert sorted(b[3]) == sorted(b_tidy[3])


def test_bravais_symbol():
    assert mxi_stills.bravais_symbol(" C 2y") == "mC"
    assert mxi_stills.bravais_symbol(" R 3 2\"") == "hR"


def test_known_symmetry_accepts_its_own_cell():
    k = mxi_stills.KnownSymmetry([78.0, 78.0, 37.0, 90, 90, 90], "-P 4 2")
    assert k.target_bravais == "tP"
    outcome, max_ang, best = k.filter(np.array([[78.0, 78.0, 37.0, 90, 90, 90],
                                                [10.0, 20.0, 30.0, 90, 90, 90]]))
    assert outcome.tolist() == [0, 1]
    assert max_ang[0] == 0.0 and np.isnan(max_ang[1])
    assert np.allclose(best[0], [78, 78, 37, 90, 90, 90])


def test_known_symmetry_errors():
    with pytest.raises(mxi_stills.SgtbxError):
        mxi_stills.sgtbx_group("Q 1")
