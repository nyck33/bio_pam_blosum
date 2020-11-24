#from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State # Load Data
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
import six.moves.urllib.request as urlreq
from six import PY3
from dash import no_update
import dash
import dash_bio as dashbio
import dash_html_components as html
import base64

#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
#import needleman
from backend.bio_needleman import Needleman
import datetime
import json
import os

from frontend.needleman_layouts.visual_layout import DATASETS

"""
When button pressed, get fastas and write a new file 
to folder sample_fastas, 
then include filepath as data in html output to Div element holding 
alignment chart

"""
#todo: cheating with global
file_path = ""


def register_visual_callbacks(app):
    """
    updates the event occurring on alignment chart
    :param app:
    :return:
    """

    @app.callback(
        Output('alignment-viewer-output', 'children'),
        [Input('alignment-chart', 'eventDatum')]
    )
    def update_output(data):
        if data is None:
            data = '{}'

        data = json.loads(data)
        if len(data.keys())==0:
            return "no data"

        return [
            html.Div(f'-{key}: {data[key]}') for key in data.keys()
        ]

    #todo: aligned-fasta-store updated on entrez-callbacks
    # Handle file upload/selection into data store
    @app.callback(
        Output('aligned-fasta-store2', 'data'),
        [Input('alignment-dropdown', 'value'),
         Input('alignment-file-upload', 'contents'),
         Input('alignment-file-upload', 'filename')]
    )
    def update_storage(dropdown, contents, filename):
        if (contents is not None) and ('fasta' in filename):
            content_type, content_string = contents.split(',')
            content = base64.b64decode(content_string).decode('UTF-8')
        else:
            content = DATASETS[dropdown]

        #jsonify
        #json_content = json.loads(content)
                    
        return content


    @app.callback(
        [Output('alignment-chart', 'data'),
        Output('desc-1', 'children'),
         Output('vs', 'children'),
        Output('desc-2', 'children')],
        [Input('aligned-fasta-store2', 'data')],
        [State('descrip-A-store', 'data'),
         State('descrip-B-store', 'data')]
    )
    def update_chart(input_data, descrip1_json, descrip2_json):
        if descrip1_json is not None:
            descrip1 = json.loads(descrip1_json)
        else:
            descrip1 = "sequence A"
        if descrip2_json is not None:
            descrip2 = json.loads(descrip2_json)
        else:
            descrip2 = "sequence B"
        return input_data, descrip1, "vs", descrip2





