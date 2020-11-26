import base64
import datetime
import io
import pandas as pd
import json

import dash_bio as dashbio
import six.moves.urllib.request as urlreq
from six import PY3
import dash_html_components as html
import plotly.express as px
import dash
from dash import Dash, exceptions, no_update, callback_context
from dash_table import DataTable
#from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State  # Load Data
from Bio import pairwise2
from Bio import SeqIO
from Bio.Align import substitution_matrices

#import layouts
from frontend.needleman_layouts.intro_layout import needleman_intro
from frontend.needleman_layouts.entrez_layout import entrez_page, matrix_names_arr
from frontend.needleman_layouts.blast_layout import blast_layout
from frontend.needleman_layouts.visual_layout import plots_page

#import register_callbacks
from frontend.needleman_callbacks.entrez_callbacks import register_entrez_callbacks
from frontend.needleman_callbacks.blast_callbacks import register_blast_callbacks
from frontend.needleman_callbacks.visual_callbacks import register_visual_callbacks
#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
from Bio import Entrez
from Bio import SeqIO
# Entrez connection
Entrez.email = "nobutaka@gatech.edu"
# tool defaults to BioPython
#Entrez.tool = "getProteinFastas" #homologs later
db = 'protein'
#todo: learn to use this and adjust retmax
#paramEutils = {'usehistory': 'Y'}
handle = Entrez.einfo(db='protein')
# record is a dictionary record['DbInfo']['FieldList] shows
#record= Entrez.read(handle)
# see all available db's

#register stylesheet
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = Dash(external_stylesheets=external_stylesheets)
app.title = "Blast Test"

app.layout=html.Div([blast_layout])

register_blast_callbacks(app)

if __name__=="__main__":
    app.run_server(debug=True, port=8080) #, dev_tools_ui=False, dev_tools_props_check=False)
    #app.run_server