import pickle
from pathlib import Path

import networkx as nx
import numpy as np
import string

from loguru import logger
from matplotlib import pyplot as plt


def create_random_graph(fixed_node_name, N=20, k1=1, k2=5) -> nx.Graph:
    # Create an empty graph
    G = nx.Graph()

    # Add the fixed node
    G.add_node(fixed_node_name, pos=(0, 0))

    # Generate N random nodes with random names from the English alphabet
    for _ in range(N):
        node_name = ''.join(np.random.choice(list(string.ascii_lowercase), size=5))
        angle = np.random.uniform(0, 2 * np.pi)
        radius = np.random.uniform(k1, k2)

        x = radius * np.cos(angle)
        y = radius * np.sin(angle)

        G.add_node(node_name, pos=(x, y))
        G.add_edge(fixed_node_name, node_name)

    return G


def store_random_graph(graph: nx.Graph) -> None:
    # Store graph object
    string_db_base_path = Path('./Temporary_files/string_db/')
    graph_obj_file_name = 'random_graph.pickle'
    graph_obj_file_path = string_db_base_path.joinpath(graph_obj_file_name)

    with open(graph_obj_file_path, 'wb') as f:
        pickle.dump(graph, f)
    logger.info(f'Random graph object successfully stored as pickle at {str(graph_obj_file_path)}')

    # Save graph as svg figure
    graph_figure_file_names = 'filtered_graph.svg'
    graph_figure_file_path = string_db_base_path.joinpath(graph_figure_file_names)

    plt.figure(figsize=(10, 10))
    nx.draw(graph, with_labels=True)
    plt.savefig(graph_figure_file_path, bbox_inches='tight')
    logger.info(f'Random graph successfully saved as svg figure at {str(graph_figure_file_path)}')


if __name__ == '__main__':
    ptx4 = create_random_graph('PTX4')
    store_random_graph(ptx4)