from biobb_common.generic.biobb_object import BiobbObject
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_common.tools.file_utils import launchlogger
import networkx as nx
import argparse
import json

class ShortestPaths(BiobbObject):
    """
    | biobb_analysis ShortestPaths
    | Calculates shortest paths in a network graph using NetworkX.
    | Calculates shortest paths in a network graph using NetworkX.

    Args:
        input_graph_path (str): Path to input graph file (GraphML). File type: input. Accepted formats: graphml.
        output_shortest_paths_path (str): Path to output shortest paths file (JSON). File type: output. Accepted formats: json.
        properties (dict - Python dictionary object containing the tool parameters, not input/output files):
            * **remove_tmp** (*bool*) - (True) Remove temporary files.
            * **restart** (*bool*) - (False) Do not execute if output files exist.

    Examples:
        from biobb_network.biobb_network.network_from_matrix import network_from_matrix
        prop = {'corr_threshold': 0.5, 'prune_mode': 'none'}
        network_from_matrix(input_matrix_path='/path/to/matrix.npz',
                            input_top_path='/path/to/top.pdb',
                            output_graph_path='/path/to/graph.graphml',
                            properties=prop)

    Info:
        * wrapped_software:
            * name: NetworkX
            * version: >=2.5
            * license: BSD 
            * url: https://networkx.org/
        * ontology:
            * name: EDAM
            * schema: http://edamontology.org/EDAM.owl
    """

    def __init__(self, input_graph_path, output_shortest_paths_path, properties=None, **kwargs) -> None:
        properties = properties or {}
        super().__init__(properties, **kwargs)
        self.locals_var_dict = locals().copy()
        self.io_dict = {
            "in": {"input_graph_path": input_graph_path},
            "out": {"output_shortest_paths_path": output_shortest_paths_path}
        }
        self.source_residue = properties.get('source_residue', None)
        self.sink_residue = properties.get('sink_residue', None)
        self.num_paths = properties.get('num_paths', 1)
        self.weight = properties.get('weight', 'weight')
        self.method = properties.get('method', 'dijkstra')
        self.properties = properties
        self.check_properties(properties)
        self.check_arguments()

    @launchlogger
    def launch(self) -> int:
        """Execute the shortest paths analysis using NetworkX."""
        fu.log(f'Executing ShortestPaths analysis', self.out_log, self.global_log)
        # Load graph
        G = nx.read_graphml(self.io_dict['in']['input_graph_path'])
        # Build mapping from node to residue number (resid)
        node_to_resid = {node: data.get('resid', node) for node, data in G.nodes(data=True)}
        resid_to_node = {data.get('resid', node): node for node, data in G.nodes(data=True)}
        # If source and sink are provided, find shortest paths between them
        if self.source_residue is not None and self.sink_residue is not None:
            source_node = resid_to_node.get(str(self.source_residue)) or resid_to_node.get(int(self.source_residue))
            sink_node = resid_to_node.get(str(self.sink_residue)) or resid_to_node.get(int(self.sink_residue))
            if source_node is None or sink_node is None:
                raise ValueError(f"Residue not found in graph nodes: source_residue={self.source_residue}, sink_residue={self.sink_residue}, method={self.method}")
            def calculate_path_weight(G, path):
                return sum(G[u][v].get(self.weight, 1.0) for u, v in zip(path[:-1], path[1:]))
            try:
                # Use shortest_simple_paths generator
                simple_paths_gen = nx.shortest_simple_paths(G, source=source_node, target=sink_node, weight=self.weight)
                all_simple_paths = []
                for i, path in enumerate(simple_paths_gen):
                    if i >= self.num_paths:
                        break
                    all_simple_paths.append(path)
            except nx.NetworkXNoPath:
                all_simple_paths = []
            # Calculate the weight of each path
            paths_with_weights = [([node_to_resid[node] for node in path], calculate_path_weight(G, path)) for path in all_simple_paths]
            # Sort paths by their total weight
            sorted_paths = sorted(paths_with_weights, key=lambda x: x[1])
            # Select the desired number of paths and store in the required format
            output_dict = {i: {"ResIDs": path, "Weight": weight} for i, (path, weight) in enumerate(sorted_paths[:self.num_paths])}
            with open(self.io_dict['out']['output_shortest_paths_path'], 'w') as f:
                json.dump(output_dict, f, indent=2)

        fu.log('ShortestPaths analysis complete.', self.out_log, self.global_log)
        return 0
    
def shortest_paths(input_graph_path: str, 
                   output_shortest_paths_path: str, 
                   properties: dict = None, 
                   **kwargs) -> int:
    """Create :class:`ShortestPaths <biobb_dyna.dyna.analysis.shortest_paths.ShortestPaths>` class and execute the :meth:`launch() <biobb_dyna.dyna.analysis.shortest_paths.ShortestPaths.launch>` method."""
    return ShortestPaths(input_graph_path=input_graph_path,
                        output_shortest_paths_path=output_shortest_paths_path,
                        properties=properties, **kwargs).launch()

def main():
    """Command line execution of this building block. Please check the command line documentation."""
    parser = argparse.ArgumentParser(
        description="Calculates shortest paths in a network graph.",
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, width=99999),
    )
    parser.add_argument(
        "-c",
        "--config",
        required=False,
        help="This file can be a YAML file, JSON file or JSON string",
    )
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-i",
        "--input_graph_path",
        required=True,
        help="Path to input graph file (GraphML). Accepted formats: graphml.",
    )
    required_args.add_argument(
        "-o",
        "--output_shortest_paths_path",
        required=True,
        help="Path to output shortest paths file (JSON). Accepted formats: json.",
    )
    args = parser.parse_args()
    config = args.config if args.config else None
    properties = settings.ConfReader(config=config).get_prop_dic()
    shortest_paths(input_graph_path=args.input_graph_path,
                  output_shortest_paths_path=args.output_shortest_paths_path,
                  properties=properties)

if __name__ == "__main__":
    main()
