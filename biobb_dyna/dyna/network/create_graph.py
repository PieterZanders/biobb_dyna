from biobb_common.generic.biobb_object import BiobbObject
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_common.tools.file_utils import launchlogger
from mdigest.core.auxiliary import prune_adjacency
import MDAnalysis as mda
import networkx as nx
import numpy as np
import argparse

class CreateGraph(BiobbObject):
    """
    | biobb_network CreateNetwork
    | Creates a network graph from a correlation matrix.
    | Creates a network graph from a correlation matrix using MDiGest utilities.

    Args:
        input_matrix_path (str): Path to input matrix file (e.g., NPZ from Dccm). File type: input. Accepted formats: npz (edam:format_1476).
        input_top_path (str): Path to topology file (e.g., PDB). File type: input. Accepted formats: pdb (edam:format_1476).
        input_contact_matrix_path (str) (optional): Path to contact matrix file (NPZ). File type: input. Accepted formats: npz (edam:format_1476).
        output_graph_path (str): Path to output graph file. File type: output. Accepted formats: graphml (edam:format_1476).
        properties (dict - Python dictionary object containing the tool parameters, not input/output files):
            * **selection** (*str*) - ("name CA") Atom selection (MDAnalysis syntax, e.g., "name CA" for alpha carbons).
            * **corr_threshold** (*float*) - (0.0) Minimum (absolute) correlation value to include an edge.
            * **use_abs** (*bool*) - (True) Use absolute values of the correlation matrix.
            * **loc_factor** (*float*) - (5.0) Distance threshold for pruning (in Angstroms).
            * **prune_mode** (*str*) - ("greater") Pruning mode: 'greater' (prune distances > loc_factor), 'lower' (prune < loc_factor), 'none' (no pruning).
            * **exclude_neighbors** (*int*) - (0) Exclude edges between nodes with index difference <= this value (e.g., for sequential residues).
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
            * name: MDiGest
            * version: >=0.1.0
            * license: Apache-2.0
        * ontology:
            * name: EDAM
            * schema: http://edamontology.org/EDAM.owl
    """

    def __init__(self, input_matrix_path, input_top_path, output_graph_path, input_contact_matrix_path=None, properties=None, **kwargs) -> None:
        properties = properties or {}
        super().__init__(properties, **kwargs)
        self.locals_var_dict = locals().copy()

        # Input/Output files
        self.io_dict = {
            "in": {"input_matrix_path": input_matrix_path, 
                   "input_top_path": input_top_path,
                   "input_contact_matrix_path": input_contact_matrix_path},
            "out": {"output_graph_path": output_graph_path}
        }

        # Properties
        self.selection = properties.get('selection', 'name CA')
        self.corr_threshold = properties.get('corr_threshold', 0.0)
        self.use_abs = properties.get('use_abs', True)
        self.loc_factor = properties.get('loc_factor', 5.0)
        self.prune_mode = properties.get('prune_mode', 'greater')
        self.exclude_neighbors = properties.get('exclude_neighbors', 0)
        self.properties = properties
        
        # Check the properties
        self.check_properties(properties)
        self.check_arguments()

    @launchlogger
    def launch(self) -> int:
        """Execute the network creation using MDiGest utilities."""
        fu.log(f'Executing NetworkFromMatrix', self.out_log, self.global_log)

        # Load matrix
        data = np.load(self.io_dict['in']['input_matrix_path'])
        matrix = data['dccm']

        # Optionally load contact matrix and multiply
        contact_matrix_path = self.io_dict['in'].get('input_contact_matrix_path')
        if contact_matrix_path:
            contact_data = np.load(contact_matrix_path)
            contact_matrix = contact_data['contact'] if 'contact' in contact_data else contact_data[list(contact_data.keys())[0]]
            matrix = matrix * contact_matrix

        # Load universe for distances and node info
        u = mda.Universe(self.io_dict['in']['input_top_path'])
        atoms = u.select_atoms(f"protein and {self.selection}")
        positions = atoms.positions
        resids = atoms.resids  # For node labels

        # Compute distance matrix
        from MDAnalysis.lib.distances import distance_array
        dist_mat = distance_array(positions, positions)

        # Process matrix
        if self.use_abs:
            matrix = np.abs(matrix)

        if self.corr_threshold > 0:
            matrix[matrix < self.corr_threshold] = 0

        if self.exclude_neighbors > 0:
            n = matrix.shape[0]
            i, j = np.indices((n, n))
            excl = np.abs(i - j) <= self.exclude_neighbors
            matrix[excl] = 0

        if self.prune_mode != 'none':
            greater = (self.prune_mode == 'greater')
            lower = (self.prune_mode == 'lower')
            matrix = prune_adjacency(matrix, dist_mat, loc_factor=self.loc_factor, greater=greater, lower=lower)

        # Create graph
        G = nx.Graph()
        n = matrix.shape[0]
        for i in range(n):
            G.add_node(i, resid=resids[i])
        for i in range(n):
            for j in range(i + 1, n):
                if matrix[i, j] != 0:
                    G.add_edge(i, j, weight=matrix[i, j])

        resid_to_node = {data['resid']: node for node, data in G.nodes(data=True)}

        # Save graph
        nx.write_graphml(G, self.io_dict['out']['output_graph_path'])

        fu.log('Network creation complete.', self.out_log, self.global_log)
        return 0
    
def create_graph(input_matrix_path: str, 
                 input_top_path: str, 
                 output_graph_path: str, 
                 input_contact_matrix_path: str = None,
                 properties: dict = None, **kwargs) -> int:
    """Create :class:`CreateGraph <biobb_dyna.dyna.network.create_graph.CreateGraph>` class and
    execute the :meth:`launch() <biobb_dyna.dyna.network.create_graph.CreateGraph.launch>` method."""
    return CreateGraph(input_matrix_path=input_matrix_path,
                       input_top_path=input_top_path,
                       output_graph_path=output_graph_path,
                       input_contact_matrix_path=input_contact_matrix_path,
                       properties=properties, **kwargs).launch()

def main():
    """Command line execution of this building block. Please check the command line documentation."""

    parser = argparse.ArgumentParser(
        description="",
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
        "--input_matrix_path",
        required=True,
        help="Path to input matrix file (e.g., NPZ from Dccm). Accepted formats: npz.",
    )
    required_args.add_argument(
        "-s",
        "--input_topology_path",
        required=True,
        help="Path to topology file (e.g., PDB). Accepted formats: pdb.",
    )
    required_args.add_argument(
        "-o",
        "--output_matrix_path",
        required=True,
        help="Path to output graph file. Accepted formats: graphml.",
    )
    parser.add_argument(
        "-ct",
        "--input_contact_matrix_path",
        required=False,
        help="Optional: Path to contact matrix file (NPZ). If given, used for pruning.",
    )

    args = parser.parse_args()
    config = args.config if args.config else None
    properties = settings.ConfReader(config=config).get_prop_dic()

    create_graph(input_matrix_path=args.input_matrix_path,
                  input_top_path=args.input_topology_path,
                  output_graph_path=args.output_matrix_path,
                  input_contact_matrix_path=args.input_contact_matrix_path,
                  properties=properties,
                  )

if __name__ == "__main__":
    main()
