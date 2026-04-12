# type: ignore
import json

import networkx as nx

from biobb_dyna.dyna.analysis.shortest_paths import shortest_paths


class TestShortestPaths:
    def test_shortest_paths_between_residues(self, tmp_path):
        G = nx.path_graph(6)
        for node in G.nodes():
            G.nodes[node]["resid"] = node + 1
        gml = tmp_path / "g.graphml"
        out = tmp_path / "sp.json"
        nx.write_graphml(G, gml)
        rc = shortest_paths(
            str(gml),
            str(out),
            properties={
                "source_residue": 1,
                "sink_residue": 6,
                "num_paths": 2,
            },
        )
        assert rc == 0
        data = json.loads(out.read_text())
        assert "0" in data
        assert "ResIDs" in data["0"]
