import dash_bio as dashbio
import six.moves.urllib.request as urlreq
from six import PY3
import dash_html_components as html
import plotly.express as px
import dash
from dash import Dash
from dash_table import DataTable
#from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State # Load Data
from Bio import SeqIO


def register_entrez_callbacks(app):
    # callbacks for entrez_page
    # callbacks for processing file uploadprocessing file upload
    @app.callback(Output('output-data-upload', 'children'),
                  [Input('upload-data', 'contents')],
                  [State('upload-data', 'filename')]
                  )
    def process_file_upload():
        pass

    # helpers for functions
    # for uplaod
    def parse_contents(contents, filename, date):
        content_typpe, content_string, = contents.split(',')

        decoded = base64.b64decode(content_string)
        try:
            if 'fasta' in filename:

