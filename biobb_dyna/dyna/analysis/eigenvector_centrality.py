from biobb_common.generic.biobb_object import BiobbObject
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_common.tools.file_utils import launchlogger
import networkx as nx
import argparse
import json

class EigenvectorCentrality(BiobbObject):
    """
    | biobb_analysis EigenvectorCentrality
    | Calculates eigenvector centrality in a network graph using NetworkX.
    | Calculates eigenvector centrality in a network graph using NetworkX.

    Args:
        input_graph_path (str): Path to input graph file (GraphML). File type: input. Accepted formats: graphml.
        output_eigenvector_path (str): Path to output eigenvector centrality file (JSON). File type: output. Accepted formats: json.
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

    def __init__(self, input_graph_path, output_eigenvector_path, properties=None, **kwargs) -> None:
        properties = properties or {}
        super().__init__(properties, **kwargs)
        self.locals_var_dict = locals().copy()
        self.io_dict = {
            "in": {"input_graph_path": input_graph_path},
            "out": {"output_eigenvector_path": output_eigenvector_path}
        }
        self.properties = properties
        self.check_properties(properties)
        self.check_arguments()

    @launchlogger
    def launch(self) -> int:
        """Execute the eigenvector centrality analysis using NetworkX."""
        fu.log(f'Executing EigenvectorCentrality analysis', self.out_log, self.global_log)
        # Load graph
        G = nx.read_graphml(self.io_dict['in']['input_graph_path'])
        # Build mapping from node to residue number (resid)
        node_to_resid = {node: data.get('resid', node) for node, data in G.nodes(data=True)}
        # Calculate eigenvector centrality
        try:
            eigenvector = nx.eigenvector_centrality(G)
        except nx.NetworkXException as e:
            fu.log(f'Error in eigenvector centrality calculation: {e}', self.out_log, self.global_log)
            eigenvector = {}
        # Output as {resid: eigenvector, ...}
        output_dict = {str(node_to_resid[node]): value for node, value in eigenvector.items()}
        with open(self.io_dict['out']['output_eigenvector_path'], 'w') as f:
            json.dump(output_dict, f, indent=2)
        fu.log('EigenvectorCentrality analysis complete.', self.out_log, self.global_log)
        return 0

def eigenvector_centrality(input_graph_path: str, 
                           output_eigenvector_path: str, 
                           properties: dict = None, 
                           **kwargs) -> int:
    """Create :class:`EigenvectorCentrality <biobb_dyna.dyna.analysis.eigenvector_centrality.EigenvectorCentrality>` class and execute the :meth:`launch() <biobb_dyna.dyna.analysis.eigenvector_centrality.EigenvectorCentrality.launch>` method."""
    return EigenvectorCentrality(input_graph_path=input_graph_path,
                                output_eigenvector_path=output_eigenvector_path,
                                properties=properties, **kwargs).launch()

def main():
    """Command line execution of this building block. Please check the command line documentation."""
    parser = argparse.ArgumentParser(
        description="Calculates eigenvector centrality in a network graph.",
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
        "--output_eigenvector_path",
        required=True,
        help="Path to output eigenvector centrality file (JSON). Accepted formats: json.",
    )
    args = parser.parse_args()
    config = args.config if args.config else None
    properties = settings.ConfReader(config=config).get_prop_dic()
    eigenvector_centrality(input_graph_path=args.input_graph_path,
                          output_eigenvector_path=args.output_eigenvector_path,
                          properties=properties)

if __name__ == "__main__":
    main()
