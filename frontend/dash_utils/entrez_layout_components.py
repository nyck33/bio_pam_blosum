import base64
import datetime
import io

import dash_bio as dashbio
import six.moves.urllib.request as urlreq
from six import PY3
import dash_html_components as html
import plotly.express as px
import dash
from dash import Dash
#from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State # Load Data
from dash_table import DataTable

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_stylesheets = [dbc.themes.BOOTSTRAP]

upload = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Label(
                "1.  To access a sequence from NCBI using accession number, enter it here"
            ),
            dcc.Input(
                id="accession-input"
            ),
            html.Label(
                "2.  Upload a fasta format sequence file"
            ),
            dcc.Upload(
                id='upload-fasta',
                children=html.Div([
                    "Drag and Drop or ",
                    html.A('Select Fasta Files')
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight':'60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'testAlign': 'center',
                    'margin': '10px'
                },
                multiple=True
                ),
            html.Label(
                    "3.  Or copy and paste or type the sequence here"
            ),
            dcc.Textarea(
                id='fasta-text-area',
                placeholder="type your fasta sequence here",
                style={'width': '75%%', 'height': 300}
            ),
            html.Button(
                'Submit', id='submit-fasta-text', n_clicks=0
            )
        ], width=6, md=6, sm=12
        ),
        dbc.Col([
            html.H3(
                "NCBI Search Results"
            ),
            DataTable(
                id="ncbi-search-res-table"
            )
        ])
    ])
])

upload_output = html.Div(
    id="upload-output"
)

progress = html.Div(
    [
        #https://dash-bootstrap-components.opensource.faculty.ai/docs/components/progress/
    ]
)
######################################################################################
#helpers for functions
# for uplaod
def parse_contents(contents, filename, date):
    content_typpe, content_string, = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'fasta' in filename:
