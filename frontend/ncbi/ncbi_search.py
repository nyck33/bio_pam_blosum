from Bio import Entrez
import xmltodict

Entrez.email = "kimn13@mytru.ca"
Entrez.tool = "getProteinFastas" #homologs later
db = 'protein'
paramEutils = {'usehistory': 'Y'}

handle = Entrez.einfo()
record= Entrez.read(handle)
record.keys()
