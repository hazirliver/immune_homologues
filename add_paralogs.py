import argparse
import os
import pickle
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from loguru import logger
from networkx.classes.graph import Graph

def parse_args() -> Tuple[Path, Path, Path]:
    parser = argparse.ArgumentParser(description="Read new paralogs' edges and add them to combined graph")

    parser.add_argument('graph', type=Path, help='Path to the pickle file containing the combined graph')
    parser.add_argument('edges', type=Path, help="Path to the .tsv file containing paralogs' edges")
    parser.add_argument('ppath', type=Path, help='Path to the file containing the initial proteins of interest')

    args = parser.parse_args()
    graph_filepath = args.graph
    edges_filepath = args.edges
    initial_proteins_path = args.ppath

    if not graph_filepath.is_file():
        raise FileNotFoundError(f'File not found: {str(graph_filepath)}')
    if not edges_filepath.is_file():
        raise FileNotFoundError(f'File not found: {str(edges_filepath)}')

    if os.path.getsize(graph_filepath) == 0:
        raise ValueError(f'File is empty: {str(graph_filepath)}')
    if os.path.getsize(edges_filepath) == 0:
        raise ValueError(f'File is empty: {str(edges_filepath)}')

    graph_filename = str(graph_filepath).split('/')[-1]
    logger.info(f'File {graph_filename} was successfully read from path {str(graph_filepath)}')

    edges_filename = str(edges_filepath).split('/')[-1]
    logger.info(f'File {edges_filename} was successfully read from path {str(edges_filepath)}')

    return Path(graph_filepath), Path(edges_filepath), Path(initial_proteins_path)


def load_graph(graph_filepath: Path) -> Graph:
    with open(graph_filepath, "rb") as f:
        return pickle.load(f)


def load_edges(edges_filepath: Path) -> pd.DataFrame:
    return pd.read_csv(edges_filepath, sep='\t', index_col=False, header=0)


def make_paralog_graph(edges: pd.DataFrame,
                       source_colname: str = 'external_gene_name',
                       target_colname: str = 'hsapiens_paralog_associated_gene_name'):
    return nx.from_pandas_edgelist(edges, source_colname, target_colname)

def combine_graphs(graph_list: List[nx.Graph]) -> nx.Graph:
    combined_graph = nx.Graph()
    for graph in graph_list:
        for node in graph.nodes():
            if node not in combined_graph:
                combined_graph.add_node(node)

        for edge in graph.edges():
            src, dest = edge
            if not combined_graph.has_edge(src, dest):
                combined_graph.add_edge(src, dest, **graph.edges[src, dest])

    logger.info(f'Total number of nodes in initial graph is {len(combined_graph.nodes)}')
    return combined_graph

def filter_by_distance(graph: nx.Graph,
                       initial_proteins: List[str],
                       min_dist_treshold: int = 1) -> nx.Graph:
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


def store_filtered_graph(filtered_graph: Graph,
                         initial_proteins_list: List[str]) -> None:
    # Store graph object
    string_db_base_path = Path('./Temporary_files/string_db/graph_with_paralogs/')
    graph_obj_file_name = 'paralogs_graph.pickle'
    graph_obj_file_path = string_db_base_path.joinpath(graph_obj_file_name)

    if not string_db_base_path.exists():
        os.makedirs(string_db_base_path)

    with open(graph_obj_file_path, 'wb') as f:
        pickle.dump(filtered_graph, f)
    logger.info(f'Graph with paralogs object successfully stored as pickle at {str(graph_obj_file_path)}')

    # Save graph as svg figure
    graph_figure_file_names = 'paralogs_graph.svg'
    graph_figure_file_path = string_db_base_path.joinpath(graph_figure_file_names)

    plt.figure(figsize=(10, 10))
    node_colors = [
        'red' if node in initial_proteins_list else 'skyblue'
        for node in filtered_graph.nodes()
    ]
    nx.draw(filtered_graph, with_labels=True, node_color=node_colors, pos=nx.spring_layout(filtered_graph,
                                                                                           k=1/10))

    plt.savefig(graph_figure_file_path, bbox_inches='tight')
    logger.info(f'Graph with paralogs graph successfully saved as svg figure at {str(graph_figure_file_path)}')

def main() -> None:
    graph_filepath, edges_filepath, initial_proteins_path = parse_args()
    combined_graph = load_graph(graph_filepath)
    edges = load_edges(edges_filepath)
    print(edges)
    initial_proteins = read_proteins_file(initial_proteins_path)
    paralog_graph = make_paralog_graph(edges)
    graph_with_paralogs = combine_graphs([combined_graph, paralog_graph])

    graph_with_paralogs_filtred = filter_by_distance(graph_with_paralogs, initial_proteins, 1)
    store_filtered_graph(graph_with_paralogs_filtred, initial_proteins)


if __name__ == '__main__':
    main()