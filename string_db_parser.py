import argparse
import os
import pickle
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


def parse_args() -> Tuple[Path, str, str]:
    """
    It takes a filepath as an argument, extracts the filename from the filepath, and then splits the filename into two
    parts, the source node and the target node
    :return: A tuple of the filepath, source_node, and target_node
    """

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

    if '_' not in filename:
        raise ValueError(f'Invalid filename format: {filename}. Must contain an underscore.')
    elif len(filename.split('_')) != 2:
        raise ValueError(f'Invalid filename format: {filename}. '
                         f'Must be in the format of "SourceProteinName_TargetProteinName".')

    source_node, target_node = filename.split('_')
    target_node = target_node.split('.')[0]
    logger.info(f'File {filename} contains {source_node} as source protein '
                f'and {target_node} as target protein in graph')

    return filepath, source_node, target_node


def make_graph(df: pd.DataFrame,
               left_node_colname: str = '#node1',
               right_node_colname: str = 'node2') -> Graph:
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

    # node_dist_threshold=5 => total length of path = 6 (4 internal nodes + source and target nodes)
    logger.info(f'All paths in the graph between {source_node} and {target_node} '
                f'whose length does not exceed {node_dist_threshold} will be found.')
    path_list_len = node_dist_threshold + 1

    logger.info('Start searching for all paths...')
    all_paths_generator = nx.all_simple_paths(graph,
                                              source=source_node,
                                              target=target_node)

    selected_paths = [path_list for path_list in all_paths_generator if len(path_list) <= path_list_len]
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
                       unique_internal_nodes: List[str],
                       source_node: str,
                       target_node: str) -> Graph:
    """
    Given a graph, a list of unique internal nodes, and a source and target node,
    return a subgraph of the original graph that contains only the unique internal nodes and the source and target nodes

    :param graph: the graph we're working with
    :type graph: Graph
    :param unique_internal_nodes: a list of nodes that are internal to the subgraph
    :type unique_internal_nodes: List[str]
    :param source_node: the node from which the subgraph is to be extracted
    :type source_node: str
    :param target_node: the node that we want to reach
    :type target_node: str
    """
    subgraph = graph.subgraph(unique_internal_nodes + [source_node, target_node])
    logger.info('Subgraph was successfully created')
    return subgraph


def save_graph(graph: Graph,
               source_node: str,
               target_node: str) -> None:
    """
    "Save the graph to a file."

    The first line of the docstring is a one-sentence summary of the object's purpose. If you have no idea what your
    function does, you won't be able to write a good docstring

    :param graph: The graph to save
    :type graph: Graph
    :param source_node: The name of the source node
    :type source_node: str
    :param target_node: The node that we want to reach
    :type target_node: str
    """
    string_db_base_path = Path('./Temporary_files/string_db/')
    if not string_db_base_path.exists():
        logger.info("The path './Temporary_files/string_db/' is used for saving results but does not exist. "
                    "It has been created.")
        string_db_base_path.mkdir(parents=True)

    # Store graph object
    graph_obj_file_name = f'{source_node}_{target_node}.pickle'
    graph_obj_file_path = string_db_base_path.joinpath(graph_obj_file_name)

    with open(graph_obj_file_path, 'wb') as f:
        pickle.dump(graph, f)
    logger.info(f'Graph object successfully stored as pickle at {str(graph_obj_file_path)}')

    # Save graph as svg figure
    graph_figure_file_names = f'{source_node}_{target_node}.svg'
    graph_figure_file_path = string_db_base_path.joinpath(graph_figure_file_names)

    plt.figure(figsize=(10, 10))
    nx.draw(graph, with_labels=True)
    plt.savefig(graph_figure_file_path, bbox_inches='tight')
    logger.info(f'Graph object successfully saved as svg figure at {str(graph_figure_file_path)}')


def main() -> None:
    filepath, source_node, target_node = parse_args()

    nodes_df = read_dataframe_from_file(filepath)
    graph = make_graph(nodes_df)
    paths_list = get_paths(graph, source_node, target_node, node_dist_threshold=5)
    unique_internal_nodes = get_unique_internal_nodes(paths_list)
    inner_subgraph = get_inner_subgraph(graph, unique_internal_nodes, source_node, target_node)
    save_graph(inner_subgraph, source_node, target_node)


if __name__ == '__main__':
    main()
