from biobb_common.generic.biobb_object import BiobbObject
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_common.tools.file_utils import launchlogger
import networkx as nx
import json
import argparse

class Communities(BiobbObject):
    """
    | biobb_analysis Communities
    | Detects communities in a network graph using NetworkX algorithms.
    | Detects communities in a network graph using NetworkX algorithms.

    Args:
        input_graph_path (str): Path to input graph file (GraphML). File type: input. Accepted formats: graphml.
        output_communities_path (str): Path to output communities file (JSON). File type: output. Accepted formats: json.
        properties (dict - Python dictionary object containing the tool parameters, not input/output files):
            * **algorithm** (*str*) - ("louvain") Community detection algorithm. Options: "louvain", "girvan_newman", "greedy_modularity", "label_propagation".
            * **resolution** (*float*) - (1.0) Resolution parameter for modularity optimization in louvain and greedy_modularity (higher values lead to more communities).
            * **max_levels** (*int*) - (None) Maximum number of partition levels to compute for hierarchical algorithms (e.g., louvain, girvan_newman). If None, computes all.
            * **remove_tmp** (*bool*) - (True) Remove temporary files.
            * **restart** (*bool*) - (False) Do not execute if output files exist.

    Examples:
        from biobb_analysis.biobb_analysis.communities import communities
        prop = {'algorithm': 'louvain', 'resolution': 1.0, 'max_levels': 5}
        communities(input_graph_path='/path/to/graph.graphml',
                    output_communities_path='/path/to/communities.json',
                    properties=prop)

    Info:
        * wrapped_software:
            * name: NetworkX
            * version: >=3.0
            * license: BSD
        * ontology:
            * name: EDAM
            * schema: http://edamontology.org/EDAM.owl
    """

    def __init__(self, input_graph_path, output_communities_path, properties=None, **kwargs) -> None:
        properties = properties or {}
        super().__init__(properties, **kwargs)
        self.locals_var_dict = locals().copy()

        # Input/Output files
        self.io_dict = {
            "in": {"input_graph_path": input_graph_path},
            "out": {"output_communities_path": output_communities_path}
        }

        # Properties
        self.algorithm = properties.get('algorithm', 'louvain')
        self.resolution = properties.get('resolution', 1.0)
        self.max_levels = properties.get('max_levels', None)
        self.properties = properties

        # Check the properties
        self.check_properties(properties)
        self.check_arguments()

    @launchlogger
    def launch(self) -> int:
        """Execute the communities analysis using NetworkX."""
        fu.log(f'Executing Communities analysis with algorithm: {self.algorithm}', self.out_log, self.global_log)

        # Load graph
        G = nx.read_graphml(self.io_dict['in']['input_graph_path'])

        # Get sorted node list (assuming nodes are integers 0 to n-1)
        node_list = sorted(G.nodes())

        # Get resids if available
        resids = [G.nodes[node].get('resid', node) for node in node_list]

        # Compute partitions based on algorithm
        if self.algorithm == 'louvain':
            partitions_gen = nx.community.louvain_partitions(G, resolution=self.resolution)
            partitions = list(partitions_gen)
        elif self.algorithm == 'girvan_newman':
            partitions_gen = nx.community.girvan_newman(G)
            partitions = list(partitions_gen)
        elif self.algorithm == 'greedy_modularity':
            partitions = [nx.community.greedy_modularity_communities(G, resolution=self.resolution)]
        elif self.algorithm == 'label_propagation':
            partitions = [list(nx.community.label_propagation_communities(G))]
        else:
            raise ValueError(f"Unknown algorithm: {self.algorithm}")

        # Limit to max_levels if specified
        if self.max_levels is not None:
            partitions = partitions[:self.max_levels]

        # Process each partition level
        results = []
        for level, partition in enumerate(partitions):
            # Compute modularity
            mod_value = nx.community.modularity(G, partition)

            # Build membership list: community index for each node in node_list order
            membership = [-1] * len(node_list)
            for idx, comm in enumerate(partition):
                for node in comm:
                    membership[node_list.index(node)] = idx

            # Community sizes
            community_sizes = [len(comm) for comm in partition]

            # Community resids dict: {comm_index: [resids]}
            community_resids = {}
            for idx in range(len(partition)):
                community_resids[idx] = [resids[i] for i in range(len(membership)) if membership[i] == idx]

            # Append to results
            results.append({
                "level": level,
                "num_communities": len(partition),
                "modularity": mod_value,
                "community_sizes": community_sizes,
                "membership": membership,
                "communities": community_resids
            })

        # Output JSON with resids and communities
        output_data = {
            "algorithm": self.algorithm,
            "communities": results
        }
        with open(self.io_dict['out']['output_communities_path'], 'w') as f:
            json.dump(output_data, f, indent=2)

        fu.log('Communities analysis complete.', self.out_log, self.global_log)
        return 0

def communities(input_graph_path: str, 
                output_communities_path: str, 
                properties: dict = None, 
                **kwargs) -> int:
    """Create :class:`Communities <biobb_dyna.analysis.communities.Communities>` class and execute the :meth:`launch() <biobb_dyna.analysis.communities.Communities.launch>` method."""
    return Communities(input_graph_path=input_graph_path,
                      output_communities_path=output_communities_path,
                      properties=properties, **kwargs).launch()

def main():
    """Command line execution of this building block. Please check the command line documentation."""
    parser = argparse.ArgumentParser(
        description="Detects communities in a network graph.",
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
        "--output_communities_path",
        required=True,
        help="Path to output communities file (CSV/TXT/JSON). Accepted formats: csv, txt, json.",
    )
    parser.add_argument(
        "-c",
        "--config",
        required=False,
        help="This file can be a YAML file, JSON file or JSON string",
    )

    args = parser.parse_args()
    config = args.config if args.config else None
    properties = settings.ConfReader(config=config).get_prop_dic()

    communities(input_graph_path=args.input_graph_path,
                output_communities_path=args.output_communities_path,
                properties=properties)

if __name__ == "__main__":
    main()
