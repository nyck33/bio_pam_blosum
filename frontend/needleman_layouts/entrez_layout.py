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

entrez_page = dbc.Container([
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
                id='upload-data',
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
                style={'width': '90%%', 'height': 300}
            ),
            html.Button(
                'Submit', id='submit-fasta-text', n_clicks=0
            )
        ], width=6, md=6, sm=12
        ),
        dbc.Col([
            html.H3("NCBI Search Results"),
            DataTable(
                id="ncbi-search-res-table"
            ),
            html.Br(),
            html.H3("Check upload contents"),
            html.Div(#show the uploaded file contents as string
                id="output-data-upload"
            ),
        ], width=6, md=6, sm=12)
    ]),
    dbc.Row([
        dbc.Col([
            html.Label(
                "Input Match Score",
                    ),
            html.Br(),
            dcc.Input(
                id="match-score-input"
            )
            ],width=4),
        dbc.Col([
            html.Label(
                'Input mismatch score'
            ),
            html.Br(),
            dcc.Input(
                id='mismatch-score-input'
            ),
        ], width=4),
        dbc.Col([
            html.Label(
                'Input gap penalty'
            ),
            html.Br(),
            dcc.Input(
                id='gap-penalty-input'
            ),
        ], width=4)
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(#https://dash-bootstrap-components.opensource.faculty.ai/docs/components/progress/
                id="progress-bar"
            )
        ], width=12)
    ])
])


######################################################################################
