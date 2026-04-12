import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import BoundaryNorm
import MDAnalysis as mda
import json
from scipy.cluster.hierarchy import dendrogram, linkage

def plot_eigenvector_centrality(eigenvector_centrality_result_path: str):

    eigenvector_data = json.load(open(eigenvector_centrality_result_path))

    residues = np.array([int(res) for res in eigenvector_data.keys()])
    eigenvector_values = np.array(list(eigenvector_data.values()))

    # Sort by residue ID for a clearer plot
    sorted_indices = np.argsort(residues)
    sorted_residues = residues[sorted_indices]
    sorted_eigenvector = eigenvector_values[sorted_indices]

    # Create a bar plot
    plt.figure(figsize=(20, 6)) # Adjust figure size for better readability of many residues
    plt.bar(sorted_residues, sorted_eigenvector)
    plt.xlabel('Residue ID')
    plt.ylabel('Eigenvector Centrality')
    plt.title('Eigenvector Centrality per Residue')
    plt.xticks(sorted_residues[::10], rotation=90) # Show ticks for every 10th residue to avoid clutter
    plt.tight_layout() # Adjust layout to prevent labels overlapping
    plt.show()

def plot_betweenness(betweenness_result_path: str):

    betweenness_data = json.load(open(betweenness_result_path))

    residues = np.array([int(res) for res in betweenness_data.keys()])
    betweenness_values = np.array(list(betweenness_data.values()))

    # Sort by residue ID for a clearer plot
    sorted_indices = np.argsort(residues)
    sorted_residues = residues[sorted_indices]
    sorted_betweenness = betweenness_values[sorted_indices]

    # Create a bar plot
    plt.figure(figsize=(20, 6)) # Adjust figure size for better readability of many residues
    plt.bar(sorted_residues, sorted_betweenness)
    plt.xlabel('Residue ID')
    plt.ylabel('Betweenness Centrality')
    plt.title('Betweenness Centrality per Residue')
    plt.xticks(sorted_residues[::10], rotation=90) # Show ticks for every 10th residue to avoid clutter
    plt.tight_layout() # Adjust layout to prevent labels overlapping
    plt.show()


def plot_shortest_paths(graph_path, input_top_path, shortest_paths_path):

    shortest_paths = json.load(open(shortest_paths_path))

    # Load the graph from the graphml file
    graph = nx.read_graphml(graph_path)

    # Load the topology to get the 3D coordinates
    u = mda.Universe(input_top_path)

    # Get the atom selection used for the matrix/network
    selection = "name CA"
    atoms = u.select_atoms(f"protein and {selection}")

    # Create a dictionary mapping residue ID to index in the atoms selection
    resid_to_index = {atom.resid: i for i, atom in enumerate(atoms)}

    # Create a dictionary of node positions using the 3D coordinates
    pos = {}
    for node_id in graph.nodes():
        resid = int(graph.nodes[node_id]['resid'])
        atom_index_in_selection = resid_to_index[resid]
        pos[node_id] = atoms.positions[atom_index_in_selection, :2]  # Use x and y for 2D plot


    # Create a figure and axes
    plt.figure(figsize=(12, 10))
    ax = plt.gca()

    # Define a colormap for different paths
    cmap = plt.cm.get_cmap('viridis')
    path_indices = list(shortest_paths.keys())
    num_paths = len(path_indices)

    # Lists to store nodes and edges to draw
    all_path_nodes = set()
    all_path_edges = []
    path_colors = []
    source_nodes = set()
    sink_nodes = set()

    # Process each path
    for idx, path_index in enumerate(path_indices):
        path_resids = shortest_paths[path_index]['ResIDs']
        # Convert residue IDs to node IDs
        path_nodes = [node_id for node_id in graph.nodes() if int(graph.nodes[node_id]['resid']) in path_resids]
        # Ensure the nodes are in the order of the path
        ordered_nodes = []
        for resid in path_resids:
            for node_id in path_nodes:
                if int(graph.nodes[node_id]['resid']) == resid:
                    ordered_nodes.append(node_id)
                    break
        # Add source (first node) and sink (last node) to their respective sets
        if ordered_nodes:
            source_nodes.add(ordered_nodes[0])
            sink_nodes.add(ordered_nodes[-1])
        # Add nodes to the path nodes set
        all_path_nodes.update(ordered_nodes)
        # Create edges for the path
        path_edges = [(ordered_nodes[i], ordered_nodes[i+1]) for i in range(len(ordered_nodes)-1)]
        all_path_edges.extend(path_edges)
        # Assign a color for this path
        path_colors.extend([cmap(idx % 10)] * len(path_edges))  # Cycle through tab10 colors if >10 paths

    # Identify non-path nodes
    all_nodes = set(graph.nodes())
    non_path_nodes = all_nodes - all_path_nodes
    # Identify intermediate path nodes (excluding source and sink)
    intermediate_path_nodes = all_path_nodes - source_nodes - sink_nodes

    # Draw non-path nodes in grey
    nx.draw_networkx_nodes(graph, pos, nodelist=list(non_path_nodes), node_size=50, node_color='grey', alpha=0.6, ax=ax)

    # Draw intermediate path nodes in lightblue
    nx.draw_networkx_nodes(graph, pos, nodelist=list(intermediate_path_nodes), node_size=50, node_color='lightblue', alpha=0.6, ax=ax)

    # Draw source nodes in red
    nx.draw_networkx_nodes(graph, pos, nodelist=list(source_nodes), node_size=100, node_color='red', alpha=0.8, ax=ax)

    # Draw sink nodes in green
    nx.draw_networkx_nodes(graph, pos, nodelist=list(sink_nodes), node_size=100, node_color='green', alpha=0.8, ax=ax)

    # Draw edges for non path nodes
    nx.draw_networkx_edges(graph, pos, edgelist=list(graph.edges()), edge_color='gray', width=0.5, alpha=0.6, ax=ax)

    # Draw edges for the paths with different colors
    nx.draw_networkx_edges(graph, pos, edgelist=all_path_edges, edge_color=path_colors, width=2, alpha=0.8, ax=ax)

    # Add a colorbar to indicate path indices
    norm = BoundaryNorm([i - 0.5 for i in range(num_paths + 1)], cmap.N)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)

    # Set ticks with at most 10
    if num_paths <= 10:
        ticks = range(num_paths)
    else:
        ticks = np.linspace(0, num_paths-1, 10, dtype=int)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([i for i in ticks])
    cbar.set_label('Path Index')

    plt.title("Shortest Paths Representation")
    plt.show()

def plot_community_modularity(communities_results_path: str):

    communities_results = json.load(open(communities_results_path))

    modularities = {}
    for i, community in enumerate(communities_results['communities']):
        modularities[community['num_communities']] = community['modularity']
    max_modularity_id = max(modularities, key=modularities.get)

    plt.figure(figsize=(10, 6))
    plt.plot(modularities.keys(), modularities.values())
    plt.xlabel('Number of Communities')
    plt.ylabel('Modularity')
    plt.axvline(x=max_modularity_id, color='r', linestyle='--', label=f'Max Modularity: {max_modularity_id:.2f}')
    plt.legend()
    plt.show()

def plot_communities(
    graph_path: str,
    input_top_path: str,
    communities_results_path: str
):
    # Load the graph from the graphml file
    graph = nx.read_graphml(graph_path)
    communities_results = json.load(open(communities_results_path))

    # Load the topology to get the 3D coordinates
    u = mda.Universe(input_top_path)

    # communities
    max_modularity_id = max(
        communities_results['communities'],
        key=lambda x: x['modularity']
    )['num_communities'] - 1

    communities = communities_results['communities'][max_modularity_id]['communities']

    # Get the atom selection used for the matrix/network
    # This should match the 'selection' property used in the Dccm
    selection = "name CA"

    # Select the corresponding atoms in the universe
    atoms = u.select_atoms(f"protein and {selection}")

    # Create a dictionary mapping residue ID to index in the atoms selection
    resid_to_index = {atom.resid: i for i, atom in enumerate(atoms)}

    # Create a dictionary of node positions using the 3D coordinates
    pos = {}
    for node_id in graph.nodes():
        # Get the residue ID from the node's attributes
        resid = int(graph.nodes[node_id]['resid'])

        # Find the corresponding atom index in the selection using the residue ID
        atom_index_in_selection = resid_to_index[resid]
        pos[node_id] = atoms.positions[atom_index_in_selection, :2] # Use x and y for 2D plot
    
    # Create a list of community indices for coloring, in the same order as graph.nodes()
    community_colors = [int(i) for i in communities.keys() for _ in communities[i]]

    # Draw the graph using the spatial positions and scaled edge widths, colored by community
    plt.figure(figsize=(12, 10))
    ax = plt.gca() # Get the current axes
    unique_communities = sorted(list(set(community_colors)))

    # Create a BoundaryNorm
    boundaries = [x - 0.5 for x in unique_communities] + [unique_communities[-1] + 0.5]
    norm = BoundaryNorm(boundaries, plt.cm.viridis.N)

    # Apply the normalization to the community colors
    nx.draw(graph, pos, with_labels=False,
            node_size=50,
            node_color=community_colors,
            cmap=plt.cm.viridis,
            edge_color='gray',
            width=1,
            alpha=0.6,
            ax=ax)

    # Add a colorbar for community indices
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
    sm.set_array([]) # Important for the colorbar to work correctly
    cbar = plt.colorbar(sm, ax=ax, ticks=unique_communities) # Set ticks to unique community indices
    cbar.set_label('Community Index')

    plt.title("Network Community Representation")
    plt.show()

def plot_community_dendrogram(communities_results_path,
                              method='average'):

    # Load the JSON output
    data = json.load(open(communities_results_path))

    # Extract resids from a singleton level if available
    n = len(data['communities'][0]['membership'])  # Number of nodes
    resids = list(range(n))  # Default fallback to node indices
    singleton_level = None
    for level in data['communities']:
        if level['num_communities'] == n:
            singleton_level = level
            break
    if singleton_level:
        for node_idx in range(n):
            comm_idx = str(singleton_level['membership'][node_idx])
            if comm_idx in singleton_level['communities']:
                resid_list = singleton_level['communities'][comm_idx]
                if len(resid_list) == 1:
                    resids[node_idx] = resid_list[0]

    # create a similarity matrix (from average membership co-occurrence across levels)
    similarity = np.zeros((n, n))
    for level in data['communities']:
        membership = np.array(level['membership'])
        for i in range(n):
            for j in range(i + 1, n):
                if membership[i] == membership[j]:
                    similarity[i, j] += 1
                    similarity[j, i] += 1
    similarity /= len(data['communities']) # Normalize

    # Convert to condensed distance (1 - similarity; use 'average' method for clustering)
    dist = 1 - similarity
    condensed_dist = dist[np.triu_indices(n, k=1)]
    Z = linkage(condensed_dist, method=method) # Methods: 'single', 'complete', 'average', etc.; 'ward' for variance minimization

    # Plot dendrogram
    plt.figure(figsize=(80, 40))
    dendrogram(Z, labels=[str(r) for r in resids])
    plt.xlabel('Residues', size=80)
    plt.ylabel('Distance', size=80)
    plt.xticks(rotation=90, size=8)
    plt.title('Community Partitioning Dendrogram', size=100)
    plt.show()

def plot_shortest_path_length_distribution(shortest_paths_path: str):

    shortest_paths = json.load(open(shortest_paths_path))

    path_lengths = [len(values['ResIDs']) for k, values in shortest_paths.items()]

    plt.figure(figsize=(10, 6))
    plt.hist(path_lengths, bins=range(min(path_lengths), max(path_lengths) + 2), align='left', rwidth=0.3)
    plt.xticks(range(min(path_lengths), max(path_lengths) + 1))
    plt.xlabel('Path Length')
    plt.ylabel('Frequency')
    plt.title('Path Length Distribution')
    plt.show()

def plot_path_frequency_per_residue(shortest_paths_path: str):

    shortest_paths = json.load(open(shortest_paths_path))

    # Count frequency of each residue in all paths
    residue_frequency = {}
    for path in shortest_paths.values():
        for resid in path['ResIDs']:
            if resid not in residue_frequency:
                residue_frequency[resid] = 0
            residue_frequency[resid] += 1

    # Prepare data for plotting
    residues = sorted(residue_frequency.keys())
    frequencies = [residue_frequency[resid] for resid in residues]

    plt.figure(figsize=(20, 6))
    plt.bar(residues, frequencies)
    plt.xlabel('Residue ID')
    plt.ylabel('Number of paths')
    plt.title('Frequency of paths passing through each residue')
    plt.xticks(residues[::10], rotation=90)  # Show ticks for every 10th residue to avoid clutter
    plt.tight_layout()  # Adjust layout to prevent labels overlapping
    plt.show()

def plot_dccm(dccm_result_path: str):

    matrix = np.load(dccm_result_path)['dccm']

    plt.figure(figsize=(10, 10))
    plt.imshow(matrix)
    plt.colorbar()
    plt.title('Correlation Matrix')
    plt.show()

def plot_network(graph_path: str,
                 input_top_path: str):
    
    # Load the graph from the graphml file
    graph = nx.read_graphml(graph_path)

    # topology
    u = mda.Universe(input_top_path)

    # Get the atom selection used for the matrix/network
    selection = "name CA"
    atoms = u.select_atoms(f"protein and {selection}")
    # Create a dictionary mapping residue ID to index in the atoms selection
    resid_to_index = {atom.resid: i for i, atom in enumerate(atoms
    )}  
    # Create a dictionary of node positions using the 3D coordinates
    pos = {}
    for node_id in graph.nodes():
        resid = int(graph.nodes[node_id]['resid'])
        atom_index_in_selection = resid_to_index[resid]
        pos[node_id] = atoms.positions[atom_index_in_selection, :2]  # Use x and y for 2D plot

    # Get edge weights and scale them for visualization
    edge_weights = [abs(graph[u][v]['weight']) for u, v in graph.edges()]
    min_width = 0
    max_width = 1
    scaled_widths = [min_width + (weight * (max_width - min_width)) for weight in edge_weights]

    # Draw the graph using the spatial positions
    plt.figure(figsize=(12, 10))
    ax = plt.gca()
    nx.draw(graph, pos, with_labels=False,
            node_size=50,
            node_color='gray',
            edge_color='gray',
            width=scaled_widths,
            alpha=0.6,
            ax=ax)
    plt.title("Network Representation")
    plt.show()