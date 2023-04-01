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


def read_dataframe_from_file(filepath: Path) -> pd.DataFrame:
    """
    Reads a tab-separated file at the specified filepath and returns a pandas DataFrame object.

    :param filepath: The path to the file to be read.
    :type filepath: Path object
    :return: A pandas DataFrame object containing the data from the file.
    """
    nodes_df = pd.read_csv(filepath, sep='\t')
    logger.info(f'There are {nodes_df.shape[0]} edges in graph')
    return nodes_df


def parse_args() -> Tuple[Path, str]:
    parser = argparse.ArgumentParser(description='Read dataframe from a file')
    parser.add_argument('filepath', type=Path, help='Path to the file containing the dataframe')
    args = parser.parse_args()
    filepath = args.filepath

    if not filepath.is_file():
        raise FileNotFoundError(f'File not found: {str(filepath)}')

    if os.path.getsize(filepath) == 0:
        raise ValueError(f'File is empty: {str(filepath)}')

    filename = str(filepath).split('/')[-1]  # extract the filename from the filepath
    logger.info(f'File {filename} was successfully read from path {str(filepath)}')

    return filepath, filename.split('.')[0]


def make_graph(df: pd.DataFrame,
               left_node_colname: str = 'preferredName_A',
               right_node_colname: str = 'preferredName_B') -> Graph:
    """
    It takes a dataframe with two columns, one for the left node and one for the right node, and returns a graph

    :param df: the dataframe containing the edge list
    :type df: pd.DataFrame
    :param left_node_colname: the name of the column in the dataframe that contains the left node, defaults to #node1
    :type left_node_colname: str (optional)
    :param right_node_colname: str = 'node2', defaults to #ode2
    :type right_node_colname: str (optional)
    :return: A graph object
    """
    if df.empty:
        raise ValueError("Dataframe cannot be empty")
    elif left_node_colname not in df.columns:
        raise ValueError(f"Column '{left_node_colname}' not found in dataframe")
    elif right_node_colname not in df.columns:
        raise ValueError(f"Column '{right_node_colname}' not found in dataframe")
    else:
        return nx.from_pandas_edgelist(df, left_node_colname, right_node_colname)


def get_paths(graph: Graph,
              source_node: str,
              target_node: str,
              node_dist_threshold: int = 5) -> List[List[str]]:
    """
    Given a graph, a source node, and a target node,
    return all paths between the source and target nodes that are at most `node_dist_threshold` nodes long

    :param graph: Graph
    :type graph: Graph
    :param source_node: the node from which we want to start the path
    :type source_node: str
    :param target_node: the node we want to get to
    :type target_node: str
    :param node_dist_threshold: The maximum number of nodes between the source and target nodes, defaults to 5
    :type node_dist_threshold: int (optional)
    :return: A list of paths between the source and target nodes.
    """
    if not isinstance(node_dist_threshold, int):
        raise ValueError("node_dist_threshold must be an integer")
    elif not 2 < node_dist_threshold < 10:
        raise ValueError("node_dist_threshold must be in the range 2 < node_dist_threshold < 10")

    logger.info(f'All paths in the graph between {source_node} and {target_node} '
                f'whose length does not exceed {node_dist_threshold} will be found.')

    logger.info('Checking existence at least one path...')
    if nx.has_path(graph, source=source_node, target=target_node):
        logger.info(f'There some path(s) between {source_node} and {target_node}. \n'
                    f'Start searching for all paths shorter than {node_dist_threshold}...')
        selected_paths = list(nx.all_simple_paths(graph,
                                                  source=source_node,
                                                  target=target_node,
                                                  cutoff=node_dist_threshold))
    else:
        logger.error(f'There are no paths between {source_node} and {target_node}. Stop running.')
        sys.exit(1)

    logger.info(f'All paths were successfully found. Total: {len(selected_paths)} paths')

    return selected_paths


def get_unique_internal_nodes(paths_list: List[List[str]]) -> List[str]:
    """
    It takes a list of paths, and returns a list of the unique internal nodes in those paths

    :param paths_list: a list of lists of strings, where each string is a node name
    :type paths_list: List[List[str]]`
    :return: a list of strings, where each string is a unique internal node name
    """
    internal_nodes_list = []
    [internal_nodes_list.extend(path[1:-1]) for path in paths_list]

    unique_internal_nodes = list(set(internal_nodes_list))
    logger.info(f'There are {len(unique_internal_nodes)} unique internal nodes')

    return unique_internal_nodes


def get_inner_subgraph(graph: Graph,
                       protein: str,
                       radius: int = 4,
                       undirected=True) -> Graph:
    subgraph = nx.ego_graph(graph, protein, radius=radius, undirected=undirected)
    logger.info('Subgraph was successfully created')
    return subgraph


def save_graph(graph: Graph,
               protein: str) -> None:
    string_db_base_path = Path('./Temporary_files/string_db/')
    if not string_db_base_path.exists():
        logger.info("The path './Temporary_files/string_db/' is used for saving results but does not exist. "
                    "It has been created.")
        string_db_base_path.mkdir(parents=True)

    # Store graph object
    graph_obj_file_name = f'{protein}.pickle'
    graph_obj_file_path = string_db_base_path.joinpath(graph_obj_file_name)

    with open(graph_obj_file_path, 'wb') as f:
        pickle.dump(graph, f)
    logger.info(f'Graph object successfully stored as pickle at {str(graph_obj_file_path)}')

    # Save graph as svg figure
    graph_figure_file_names = f'{protein}.svg'
    graph_figure_file_path = string_db_base_path.joinpath(graph_figure_file_names)

    plt.figure(figsize=(10, 10))
    node_colors = [
        'red' if node == protein else 'skyblue'
        for node in graph.nodes()
    ]
    nx.draw(graph, with_labels=True, node_color=node_colors)
    plt.savefig(graph_figure_file_path, bbox_inches='tight')
    logger.info(f'Graph object successfully saved as svg figure at {str(graph_figure_file_path)}')


def main() -> None:
    filepath, protein = parse_args()

    nodes_df = read_dataframe_from_file(filepath)
    graph = make_graph(nodes_df)
    inner_subgraph = get_inner_subgraph(graph, protein)
    save_graph(inner_subgraph, protein)


if __name__ == '__main__':
    main()
