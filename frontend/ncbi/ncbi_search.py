from Bio import Entrez
from Bio import SeqIO
import xmltodict

def get_last_updated(record):
    return record['DbInfo']['LastUpdate']

def searchByTerm(term):
    """
    search by term to get a list of accession numbers next to description
    :param term:
    :return: queryTranslation and accession numbers
    """
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
    if full_fasta:
        return fasta_handle.read()
    seq_str = str(rec_fasta.seq)
    return seq_str

