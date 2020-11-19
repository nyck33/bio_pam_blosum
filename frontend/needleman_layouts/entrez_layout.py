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

# for dropdown
matrix_names_arr = ['PAM 30', 'PAM 70', 'PAM 250',
                    'BLOSUM 80', 'BLOSUM 62',
                     'BLOSUM 45', 'BLOSUM 50', 'BLOSUM 90']

entrez_page = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Label(
                "1.  To access a sequence from NCBI using accession number, enter it here"
            ),
            dcc.Input(
                id="accession-input",
                type="text",
                size="50"
            ),
            html.Br(),
            html.Label(
                "Enter accession number for second sequence here"
            ),
            dcc.Input(
                id="accession-input-two",
                type="text",
                size="50"
            ),
            html.Br(),
            html.Label(
                "Search by term (including booleans), click on the link for boolean rules."
            ),
            html.Br(),
            html.A("NCBI boolean rules",
                href="https://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/Entrez/complex_boolean.html",
                target="_blank"
                   ),
            html.Br(),
            html.Label(
                "Enter terms using boolean rules here:"
            ),
            dcc.Input(
                id="term-search-input",
                type="text",
                size="100"
            ),
            html.Br(),
            html.Hr(),
            html.Label("NCBI Search Results appear here"),
            html.H3("NCBI Search Results"),
            DataTable(
                id="ncbi-search-res-table"
            ),
            html.Br(),
            html.Label(
                "2.  Upload a single or multi-fasta file"
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
            html.Br(),
            html.Label(
                "Upload a second fasta file (if not multi)"
            ),
            dcc.Upload(
                id='upload-data-two',
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
            html.Br(),
            html.H3("Check upload contents sequence 1"),
            html.Div(#show the uploaded file contents as string
                id="output-data-upload"
            ),
            html.Br(),
            html.H3("Check upload contents (if) sequence 2"),
            html.Div(#show the uploaded file contents as string
                id="output-data-upload-two"
            ),
            html.Br(),
            html.Label(
                    "3.  Or copy and paste or type sequence 1 here"
            ),
            html.Br(),
            dcc.Textarea(
                id='fasta-text-area',
                placeholder="type your fasta sequence here",
                style={'width': '100%', 'height': 300}
            ),
            html.Br(),
            html.Label(
                    "3.  Or copy and paste or type sequence 2 here"
            ),
            html.Br(),
            dcc.Textarea(
                id='fasta-text-area-two',
                placeholder="type your fasta sequence here",
                style={'width': '100%', 'height': 300}
            ),
            html.Br(),

        ], width=12, md=12, sm=12
        ), #end col
    ]),
    html.Hr(),
    dbc.Row([
        dbc.Col([
            html.Label("Input Match Score"),
            html.Br(),
            dcc.Input(
                id="match-score-input"
            )
            ],width=4),
        dbc.Col([
            html.Label('Input mismatch score'),
            html.Br(),
            dcc.Input(
                id='mismatch-score-input'
            ),
        ], width=4),
        dbc.Col([
            html.Label('Input gap penalty, default=10'),
            html.Br(),
            dcc.Input(
                id='gap-penalty-input'
            ),
        ], width=4)
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id='matrix-dropdown',
                options = [{'label': matrix, 'value': matrix} for matrix in matrix_names_arr
                            ],
                value="BLOSUM 62"
            )
        ],width=12),
        dbc.Col([
            html.Br(),
            dbc.Button(
                'Submit', id='submit-fasta', n_clicks=0
            )
        ], width=12)
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
