# type: ignore
import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("MDAnalysis")
pytest.importorskip("mdigest.core.auxiliary")

from biobb_dyna.dyna.network.create_graph import create_graph, CreateGraph


class TestCreateGraph:
    def test_create_graph_full_correlation(self, pdb_path, tmp_path):
        import MDAnalysis as mda

        u = mda.Universe(pdb_path)
        n = len(u.select_atoms("protein and name CA"))
        assert n > 1
        rng = np.random.default_rng(0)
        mat = rng.random((n, n))
        mat = (mat + mat.T) / 2.0
        np.fill_diagonal(mat, 1.0)
        npz = tmp_path / "m.npz"
        np.savez(npz, dccm=mat)
        out = tmp_path / "g.graphml"
        rc = create_graph(
            str(npz),
            pdb_path,
            str(out),
            properties={"prune_mode": "none", "corr_threshold": 0.0},
        )
        assert rc == 0
        assert out.is_file()
        text = out.read_text()
        assert "graphml" in text.lower()

    def test_create_graph_class_construct(self, pdb_path, tmp_path):
        npz = tmp_path / "m.npz"
        n = 4
        mat = np.eye(n)
        np.savez(npz, dccm=mat)
        out = tmp_path / "g.graphml"
        bb = CreateGraph(
            str(npz),
            pdb_path,
            str(out),
            properties={"prune_mode": "none"},
        )
        assert bb.corr_threshold == 0.0
