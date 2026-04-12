# type: ignore
import os

import pytest

pytest.importorskip("mdigest.core.correlation")

import biobb_dyna.dyna.dyncorr.dccm as dccm_mod
from biobb_dyna.dyna.dyncorr.dccm import Dccm, compute_dccm


class TestDccm:
    def test_module_exports(self):
        assert dccm_mod.Dccm is Dccm
        assert callable(compute_dccm)

    def test_dccm_integration(self, tmp_path):
        traj = os.environ.get("BIOB_DYNA_TEST_TRAJ")
        top = os.environ.get("BIOB_DYNA_TEST_TOP")
        if not traj or not top or not os.path.isfile(traj) or not os.path.isfile(top):
            pytest.skip(
                "Set BIOB_DYNA_TEST_TRAJ and BIOB_DYNA_TEST_TOP to run Dccm integration test."
            )
        out = tmp_path / "dccm.npz"
        rc = compute_dccm(
            traj,
            top,
            str(out),
            properties={
                "matrix_type": "correlation",
                "dynamics_type": "displacements",
                "stride_step": 50,
            },
        )
        assert rc == 0
        assert out.is_file()
