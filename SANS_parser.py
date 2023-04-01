import random
import time

import pandas as pd
from Bio import SeqIO
import requests
from bs4 import BeautifulSoup
from loguru import logger

PROTS_ORDER_1 = ['ZC3HAV1.faa', 'EXOSC5', 'EXOSC7', 'EXOSC3',
                 'TRIM25', 'DDX58', 'JAK2', 'CAMK2A', 'CAMK2B',
                 'CAMK2D', 'CAMK2G', 'ARAF', 'RAF1']

# Set the path to your multifasta file
MULTIFASTA_PATH = "./Input_data/string_protein_sequences.fa"


def sleep_random_time():
    time_to_sleep = random.uniform(5, 15)
    time.sleep(time_to_sleep)


def read_fasta_aa(prot_name: str) -> str:
    multifasta_obj = SeqIO.parse(MULTIFASTA_PATH, "fasta")
    for record in multifasta_obj:
        if prot_name in record.description:
            return str(record.seq)


def post_to_dataframe(post_content):
    soup = BeautifulSoup(post_content, 'html.parser')
    table = soup.find('table', {'border': True})
    if table is not None:
        return pd.read_html(str(table))[0]
    print("No <table> found in POST content.")
    return None


def sans_post_content(fasta_aa: str,
                      config: dict[str, str]) -> bytes:
    SANS_URL = "http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/sans.cgi"
    config['seq'] = fasta_aa
    try:
        sleep_random_time()
        response = requests.post(SANS_URL, data=config)
        return response.content
    except Exception as e:
        logger.error('Something went wrong while request to SANS')
        logger.error(f'Error: {e}')


def SANS(prot_list: list[str]) -> None:
    for i, prot_name in enumerate(prot_list):
        logger.info(f'----- Start processing {prot_name} -----')
        fasta_aa = read_fasta_aa(prot_name)
        logger.info(f'1. Successfully read fasta for {prot_name}. Protein length is {len(fasta_aa)}')
        request_content = sans_post_content(fasta_aa, config)
        logger.info(f'2. Successfully got response from SANS for {prot_name}')
        sans_df = post_to_dataframe(request_content)
        logger.info(f'3. Successfully converted response content to dataframe for {prot_name}')
        sans_df.to_csv(f'./Temporary_files/SANS/{prot_name}.tsv', sep='\t')
        logger.info(f'4. Successfully saved dataframe for {prot_name} to ./Temporary_files/')
        logger.info(f'[{i + 1} / {len(prot_list)}]}} ----- Protein {prot_name} successfully processed -----')


if __name__ == '__main__':
    config = {
        'mode': 'table',
        'db': 'uniprot',
        'H': '100000',
        'protocol': '2',
    }
    logger.info(f"Configuration: {config}")
    SANS(PROTS_ORDER_1)
