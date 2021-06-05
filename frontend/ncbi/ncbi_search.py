from Bio import Entrez
from Bio import SeqIO
import xmltodict
import json

user_email = ''
Entrez.email = user_email
db = 'protein'
handle = Entrez.einfo(db=db)

def set_global_email_ncbi_search(email_str):
    global user_email
    user_email = email_str
def get_last_updated(record):
    return record['DbInfo']['LastUpdate']

def searchByTerm(term):
    """
    search by term to get a list of accession numbers next to description
    :param term:
    :return: queryTranslation and accession numbers
    """
    #todo: move this email inside not global
    #Entrez.email = "nobutaka@gatech.edu"

    db="protein"
    eSearch = Entrez.esearch(db=db, term=term, idtype="acc")
    res = Entrez.read(eSearch)
    accessions_arr = res['IdList']
    query_translation = res['QueryTranslation']

    return accessions_arr, query_translation


def get_full_GB_info(accession):
    """
    Full info with journals, versions, source, organism, reference authors
    features (source, gene, CDS), origin=sequence
    :return: str of gb info
    """
    gb = Entrez.efetch(db='protein', id=str(accession), rettype='gb', retmode='text')
    gb_str = gb.read()
    return gb_str


def get_fasta_by_accession(accession, full_fasta=False):
    """
    SeqIO.read returns a fasta object vs fasta_handle.read() returns full fasta string
    :param accession: acccession as string
    :param full_fasta: if True return full fasta as string
    :return: sequence string
    """
    fasta_handle = Entrez.efetch("protein", id=accession, rettype='fasta', retmode='text' )
    rec_fasta = SeqIO.read(fasta_handle, "fasta")


    seq_str = str(rec_fasta.seq)
    fasta_file_str = ">"

    if full_fasta:
        # return description as str
        fasta_file_str += rec_fasta.description
        fasta_file_str += "\n"
        return fasta_file_str

    return seq_str

