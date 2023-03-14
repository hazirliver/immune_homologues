import argparse
import os
import pickle
import sys
from pathlib import Path
from typing import List, Tuple

from functools import reduce
import networkx as nx
from loguru import logger
from matplotlib import pyplot as plt
from networkx.classes.graph import Graph


def parse_args() -> Tuple[Path, Path]:
    """
    It parses the command line arguments and returns the path to the directory with stored graphs
    :return: The path to the directory where the graphs are stored.
    """
    parser = argparse.ArgumentParser(description='Read path with stored graphs')
    parser.add_argument('gpath', type=Path, help='Path to the strored graphs')
    parser.add_argument('ppath', type=Path, help='Path to the file containing the initial proteins of interest')
    args = parser.parse_args()

    if not sys.argv[1:]:
        raise ValueError("Please provide a path to the stored graphs.")
    if not args.gpath.exists():
        raise FileNotFoundError(f"Provided path '{args.gpath}' does not exist.")
    if not args.gpath.is_dir():
        raise NotADirectoryError(f"Provided path '{args.gpath}' is not a directory.")
    if not os.access(args.gpath, os.R_OK):
        raise PermissionError(f"You do not have read permissions for the file '{args.gpath}'.")

    graphs_path = args.gpath
    initial_proteins_path = args.ppath
    logger.info(f'Path to stored graphs: {str(graphs_path)}')

    return graphs_path, initial_proteins_path


def load_graphs_from_folder(folder_path: Path) -> List[Graph]:
    """
    "Load networkx graphs from a folder containing pickle files."

    :param folder_path: The path to the folder containing the pickle files
    :type folder_path: Path
    :return: A list of networkx graphs read from the folder.
    """
    graphs = []

    for file in folder_path.glob("*.pickle"):
        with open(file, "rb") as f:
            graph = pickle.load(f)
            if isinstance(graph, nx.Graph):
                graphs.append(graph)
    logger.info(f'Successfully loaded {len(graphs)} graphs')
    return graphs


def combine_graphs(graphs_list: List[Graph]) -> Graph:
    """
    Combine multiple networkx graphs into a single graph.

    :param graphs_list: A list of networkx graphs to be combined
    :type graphs_list: List[Graph]
    :return: A single networkx graph resulting from the combination of all graphs in the input list.
    """
    combined_graph = reduce(nx.compose, graphs_list)
    logger.info(f'Successfully combined {len(graphs_list)} graphs into one. '
                f'Total number of nodes is {len(combined_graph.nodes)}')
    return combined_graph


def filter_by_distance(graph: nx.Graph,
                       initial_proteins: List[str],
                       min_dist_treshold: int = 2) -> nx.Graph:
    """
    It removes all vertices that are more than N steps away from all special vertices

    :param graph: The input graph
    :type graph: nx.Graph
    :param initial_proteins: A list of initial proteins node names
    :type initial_proteins: List[str]
    :param min_dist_treshold: The distance threshold. Vertices farther than N steps from all special vertices will be removed, defaults to 5
    :type min_dist_treshold: int (optional)
    :return: A filtered graph with vertices more than N steps away from all special vertices removed.
    """
    vertices_to_remove = set()

    for protein in graph.nodes():
        if protein in initial_proteins:
            continue

        min_distance = float('inf')

        for initial_protein in initial_proteins:
            if nx.has_path(graph, protein, initial_protein):
                distance = nx.shortest_path_length(graph, protein, initial_protein)
                min_distance = min(min_distance, distance)
        if min_distance > min_dist_treshold:
            vertices_to_remove.add(protein)

    filtered_graph = graph.copy()
    filtered_graph.remove_nodes_from(vertices_to_remove)

    logger.info(f'Combined graph successfully filtered with minimal distance treshold={min_dist_treshold}. '
                f'Total number of nodes in filtered graph is {len(filtered_graph.nodes)}')

    logger.info(f'Filtered graph is connected: {nx.is_connected(filtered_graph)}')

    return filtered_graph


def read_proteins_file(protein_filepath: Path) -> List[str]:
    """
    The function reads a file containing proteins of interest, with the file path specified as a command line argument.
    It then returns a list of these proteins.

    :param protein_filepath: The path to the file to read
    :type protein_filepath: Path
    :return: A list of proteins.
    """
    with open(protein_filepath, 'r') as protein_file:
        lines = protein_file.readlines()
        protein_list = [line.strip() for line in lines]
        logger.info(f'There are {len(protein_list)} proteins in {protein_filepath}. Specifically:')
        logger.info('\n' + '\n'.join(protein_list))
        return protein_list


def store_filtered_graph(filtered_graph: Graph) -> None:
    # Store graph object
    string_db_base_path = Path('./Temporary_files/string_db/')
    graph_obj_file_name = 'filtered_graph.pickle'
    graph_obj_file_path = string_db_base_path.joinpath(graph_obj_file_name)

    with open(graph_obj_file_path, 'wb') as f:
        pickle.dump(filtered_graph, f)
    logger.info(f'Filtered graph object successfully stored as pickle at {str(graph_obj_file_path)}')

    # Save graph as svg figure
    graph_figure_file_names = 'filtered_graph.svg'
    graph_figure_file_path = string_db_base_path.joinpath(graph_figure_file_names)

    plt.figure(figsize=(10, 10))
    nx.draw(filtered_graph, with_labels=True)
    plt.savefig(graph_figure_file_path, bbox_inches='tight')
    logger.info(f'Filtered graph successfully saved as svg figure at {str(graph_figure_file_path)}')


def main() -> None:
    graphs_path, initial_proteins_path = parse_args()
    graphs_list = load_graphs_from_folder(graphs_path)
    initial_proteins_list = read_proteins_file(initial_proteins_path)

    combined_graph = combine_graphs(graphs_list)
    filtered_graph = filter_by_distance(combined_graph, initial_proteins_list)
    store_filtered_graph(filtered_graph)

if __name__ == '__main__':
    main()
