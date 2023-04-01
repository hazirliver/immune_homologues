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


def parse_args() -> Tuple[Path, Path, Path, Path]:
    parser = argparse.ArgumentParser(description="Read new paralogs' edges and add them to combined graph")

    parser.add_argument('graph', type=Path, help='Path to the pickle file containing the combined graph with paralogs')
    parser.add_argument('orthologs', type=Path, help='Path to the file containing info about orthologs')
    parser.add_argument('taxids', type=Path, help='Path to the file taxids for holologs filtering')
    parser.add_argument('ppath', type=Path, help='Path to the file containing the initial proteins of interest')


    args = parser.parse_args()
    graph_filepath = args.graph
    orthologs_filepath = args.orthologs
    taxids_filepath = args.taxids
    initial_proteins_path = args.ppath

    if not graph_filepath.is_file():
        raise FileNotFoundError(f'File not found: {str(graph_filepath)}')
    if not orthologs_filepath.is_file():
        raise FileNotFoundError(f'File not found: {str(orthologs_filepath)}')

    if os.path.getsize(graph_filepath) == 0:
        raise ValueError(f'File is empty: {str(graph_filepath)}')
    if os.path.getsize(orthologs_filepath) == 0:
        raise ValueError(f'File is empty: {str(orthologs_filepath)}')

    graph_filename = str(graph_filepath).split('/')[-1]
    logger.info(f'File {graph_filename} was successfully read from path {str(graph_filepath)}')

    edges_filename = str(orthologs_filepath).split('/')[-1]
    logger.info(f'File {edges_filename} was successfully read from path {str(orthologs_filepath)}')

    return Path(graph_filepath), Path(orthologs_filepath), Path(taxids_filepath), Path(initial_proteins_path)


def load_graph(graph_filepath: Path) -> Graph:
    with open(graph_filepath, "rb") as f:
        return pickle.load(f)

def load_ortholog_df(orthologs_filepath: Path) -> pd.DataFrame:
    return pd.read_csv(orthologs_filepath, sep='\t', index_col=False, header=0)

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

def load_taxids(taxids_filepath: Path) -> List[int]:
    with open(taxids_filepath, 'r') as taxids_file:
        lines = taxids_file.readlines()
        return [int(line.strip()) for line in lines if line.strip() != '']


def filter_ortholods_df(orthologs_df: pd.DataFrame,
                        taxids_list: List[int]) -> pd.DataFrame:
    orthologs_df_hs = orthologs_df.loc[orthologs_df['#tax_id'] == 9606, :]
    return orthologs_df_hs.loc[orthologs_df['Other_tax_id'].isin(taxids_list), :]


def get_nodes(graph: Graph) -> List[str]:
    return list(graph.nodes)


def get_orthologs_edges(orthologs_df_filtred: pd.DataFrame,
                        graph_nodes: List[str]):
    return orthologs_df_filtred.loc[orthologs_df_filtred[orthologs_df_filtred['GeneID'].isin(graph_nodes)],
                                    ['GeneID', 'Other_GeneID']]


def create_orthologs_graph(orthologs_edges: pd.DataFrame,
                           source_colname: str = 'GeneID',
                           target_colname: str = 'Other_GeneID') -> Graph:
    return nx.from_pandas_edgelist(orthologs_edges, source_colname, target_colname)


def main() -> None:
    graph_filepath, orthologs_filepath, taxids_filepath, initial_proteins_path = parse_args()
    graph = load_graph(graph_filepath)
    orthologs_df = load_ortholog_df(orthologs_filepath)
    taxids_list = load_taxids(taxids_filepath)
    initial_proteins = read_proteins_file(initial_proteins_path)

    orthologs_df_filtred = filter_ortholods_df(orthologs_df, taxids_list)
    graph_nodes = get_nodes(graph)
    orthologs_edges = get_orthologs_edges(orthologs_df_filtred, graph_nodes)
    orthologs_graph = create_orthologs_graph(orthologs_edges)

if __name__ == '__main__':
    main()