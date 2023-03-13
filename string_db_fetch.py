import argparse
import os
import sys
from itertools import combinations
from pathlib import Path
from random import uniform
from time import sleep
from typing import Tuple, List, TypeAlias
from io import StringIO

import pandas as pd
import requests
from loguru import logger
from requests import Response

STRING_API_URL = 'https://version-11-5.string-db.org/api'
STRING_PARAMS = {
    'species': 9606,
    'required_score': 400,
    'network_type': 'functional',
    'add_nodes': 50,
    'show_query_node_labels': 1,
}

ProteinCombinations: TypeAlias = List[Tuple[str, ...]]


def parse_args() -> Path:
    """
    It reads the file with the proteins of interest, checks if it's not empty and contains at least 2 rows
    :return: The filepath of the file containing the proteins of interest
    """

    parser = argparse.ArgumentParser(description='Read initial file with proteins of interest')
    parser.add_argument('filepath', type=Path, help='The path to the file containing the proteins of interest')
    args = parser.parse_args()
    filepath = args.filepath

    if not filepath.is_file():
        raise FileNotFoundError(f'File not found: {str(filepath)}')

    if os.path.getsize(filepath) == 0:
        raise ValueError(f'File is empty: {str(filepath)}')
    with open(filepath, 'r') as file:
        if len(file.readlines()) < 2:
            raise ValueError('File contain only 1 row. Minimum number of proteins of interest is 2')

    filename = str(filepath).split('/')[-1]
    logger.info(f'File {filename} was successfully read from path {str(filepath)}')

    return filepath


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
        logger.info(f'There are {len(protein_list)} in {protein_filepath}. Specifically:')
        logger.info('\n'.join(protein_list))
        return protein_list


def get_all_N_combinations(filepath: Path,
                           N: int = 2) -> ProteinCombinations:
    """
    It reads the protein names from the file, and then returns all possible unique combinations of N proteins

    :param filepath: The path to the file containing the protein names
    :type filepath: Path
    :param N: The number of proteins in each combination, defaults to 2, defaults to 2
    :type N: int (optional)
    :return: A list of tuples of all possible unique combinations of N proteins.
    """
    protein_names = read_proteins_file(filepath)
    protein_N_combinations = list(combinations(protein_names, N))
    logger.info(f'For {len(protein_names)} initial proteins, '
                f'{len(protein_N_combinations)} combinations of {N} proteins are obtained.')
    return protein_N_combinations


def write_response_tsv(response: Response, paroteins_filepath: Path) -> None:
    """
    It takes a response object from the requests library and writes the content to a file

    :param response: Response
    :type response: Response
    :param paroteins_filepath: The path to the file where the data will be saved
    :type paroteins_filepath: Path
    """
    content = response.content.decode('utf-8')
    # noinspection PyTypeChecker
    df = pd.read_csv(StringIO(content), sep='\t')
    df.drop_duplicates(inplace=True)

    df.to_csv(paroteins_filepath, sep='\t', index=False)

def process_one_combination_string_db(protein_combination: Tuple[str, ...]) -> None:
    """
    It takes a tuple of proteins and writes the response from the String-DB API to a file

    :param protein_combination: Tuple[str, ...]
    :type protein_combination: Tuple[str, ...]
    """
    output_format = 'tsv'
    method = 'network'
    request_url = '/'.join([STRING_API_URL, output_format, method])

    params = {**STRING_PARAMS, "identifiers": "%0d".join(protein_combination)}

    try:
        sleep(uniform(0.5, 5))

        response = requests.post(request_url, data=params)

        response.raise_for_status()  # raise an exception for 4xx or 5xx status codes
        logger.info(f"Successfully received a response from the String-DB "
                    f"for proteins [{', '.join(protein_combination)}]")
    except requests.exceptions.RequestException as e:
        # Catch any request exceptions and handle them appropriately
        logger.error(f"Error making request for [{', '.join(protein_combination)}]: {e}")
        sys.exit(1)

    paroteins_filepath = Path(f"./Temporary_files/string_db/{'_'.join(protein_combination)}.tsv")
    write_response_tsv(response, paroteins_filepath)
    logger.info(f"String-DB results for [{', '.join(protein_combination)}] saved to {str(paroteins_filepath)}")


def fetch_all_string_db(protein_combinations: ProteinCombinations) -> None:
    """
    For each protein combination, fetch the interactions from STRING-DB and save them to a file

    :param protein_combinations: a list of tuples of proteins
    :type protein_combinations: ProteinCombinations
    """
    number_of_combinations = len(protein_combinations)
    for id_, protein_combination in enumerate(protein_combinations, start=1):
        logger.info(f"[{id_}/{number_of_combinations}]. Start processing [{', '.join(protein_combination)}]")
        process_one_combination_string_db(protein_combination)
        logger.info(f"[{id_}/{number_of_combinations}]. Successfully processed")


def main():
    logger.info('Start fetching String-DB')
    protein_of_interest_filepath = parse_args()
    protein_combinations = get_all_N_combinations(protein_of_interest_filepath, N=2)
    fetch_all_string_db(protein_combinations)
    logger.info('All results successfully fetched.')


if __name__ == '__main__':
    main()
