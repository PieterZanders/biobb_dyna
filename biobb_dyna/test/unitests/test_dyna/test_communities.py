# type: ignore
import json

import networkx as nx
import pytest

from biobb_dyna.dyna.analysis.communities import communities, Communities


def _write_simple_graphml(path):
    G = nx.karate_club_graph()
    for i, node in enumerate(G.nodes()):
        G.nodes[node]["resid"] = i + 1
    nx.write_graphml(G, path)


class TestCommunities:
    def test_communities_louvain(self, tmp_path):
        gml = tmp_path / "g.graphml"
        out = tmp_path / "c.json"
        _write_simple_graphml(gml)
        rc = communities(
            str(gml),
            str(out),
            properties={"algorithm": "louvain", "resolution": 1.0},
        )
        assert rc == 0
        data = json.loads(out.read_text())
        assert data["algorithm"] == "louvain"
        assert "communities" in data
        assert len(data["communities"]) >= 1

    def test_communities_label_propagation(self, tmp_path):
        gml = tmp_path / "g.graphml"
        out = tmp_path / "c.json"
        _write_simple_graphml(gml)
        rc = communities(
            str(gml),
            str(out),
            properties={"algorithm": "label_propagation"},
        )
        assert rc == 0

    def test_communities_unknown_algorithm(self, tmp_path):
        gml = tmp_path / "g.graphml"
        out = tmp_path / "c.json"
        _write_simple_graphml(gml)
        with pytest.raises(ValueError, match="Unknown algorithm"):
            Communities(
                str(gml),
                str(out),
                properties={"algorithm": "not_a_real_algo"},
            ).launch()
