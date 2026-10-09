import mxi_stills


def test_build_info():
    assert mxi_stills.__version__
    assert mxi_stills.git_revision and mxi_stills.git_revision != "unknown"
    assert mxi_stills.compiler
