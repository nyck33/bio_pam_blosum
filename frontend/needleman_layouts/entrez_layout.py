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
    #########################################################################
    # debug and store
    dbc.Row([
        dbc.Col([
            dcc.Store(id="sequence-1-store"),
            html.Div(
                "check store 1",
                id="check-store-1"
            ),
            dcc.Store(id="sequence-2-store"),
            html.Div(
                "check store 2",
                id="check-store-2"
            ),
            html.Div(
                "context checker",
                id='ctx-check'
            ),
            html.Div(  # receive preformatted text html.Pre, for uploads
                "other debug output",
                id='debug-out'
            )
        ], width=12)
    ]),
    #######################################################################
    #accession input
    dbc.Row([
        dbc.Col([
            html.H3(
                "1.  To access a sequence from NCBI using accession number, enter it here"
            ),
            html.Br(),
            html.Label("Accession 1:"),
            dcc.Input(
                id="accession-input",
                type="text",
                size="50"
            ),
            dbc.Button(
                'Search 1',
                id='btn-acc-1',
                color="primary",
                n_clicks=0
            ),
            html.Br(),
            html.Label("Accession 2:"),
            dcc.Input(
                id="accession-input-two",
                type="text",
                size="50"
            ),
            dbc.Button(
                'Search 2',
                id='btn-acc-2',
                color="secondary",
                n_clicks=0
            )

        ], width=12),
    ]),
    dbc.Row([ # results of accession search
        dbc.Col([
            html.P(
                "fasta res 1 here",
                id='accession-fasta-1',
            ),
        ], width=6),
        dbc.Col([
            html.P(
                "fasta res 2 here",
                id='accession-fasta-2'
            )
        ], width=6)
    ]),
    html.Hr(),
    ################################################################
    #NCBI search area
    dbc.Row([
        dbc.Col([
            html.H3(
                "2. Search by term (including booleans)"
            ),
            html.Br(),
            html.Label("click on the link for boolean rules."),
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
            dbc.Button(
                'Find Accession Numbers',
                id='btn-search-ncbi',
                color="primary",
                n_clicks=0
            ),

        ], width=12)
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("NCBI Search by Text Results appear here"),
            html.Div(
                # full fasta for accession num input
                id="ncbi-textsearch-res"
            ),
            html.Br(),
            html.P(
                "query translation output",
                id="query-translation-output"
            ),
        ], width=12)
    ]),
    html.Hr(),
    dbc.Row([
        dbc.Col([
            html.H3("3. Upload 2 fasta files or a multi-fasta file of 2 sequences")
        ], width=12),
        dbc.Col([
            html.Label("fasta file 1"),
            dcc.Upload(
                id='upload-data',
                children=html.Div([
                    "Drag and Drop or ",
                    html.A('Select Fasta Files')
                ]),
                style={
                    'width': '70%',
                    'height': '100px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'testAlign': 'center',
                    'margin': '10px'
                },
                multiple=False
            ),
        ], width=6),
        dbc.Col([
            html.Label("fasta file 2"),
            dcc.Upload(
                id='upload-data-two',
                children=html.Div([
                    "Drag and Drop or ",
                    html.A('Select Fasta Files')
                ]),
                style={
                    'width': '70%',
                    'height': '100px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'testAlign': 'center',
                    'margin': '10px'
                },
                multiple=False
            ),
        ], width=6)
    ]),
    dbc.Row([
        dbc.Col([
            html.P("Check upload contents sequence 1"),
            html.Div(  # show the uploaded file contents as string
                id="output-data-upload"
            ),
        ], width=12),
        dbc.Col([
            html.P("Check upload contents (if) sequence 2"),
            html.Div(  # show the uploaded file contents as string
                id="output-data-upload-two"
            ),
        ], width=12)
    ]),
    html.Hr(),
    dbc.Row([
        dbc.Col([
            html.H3("3. Manual Input of Sequences, no headers, sequences only"),
            html.P(
                "Copy and Paste or Type",
                id="text-instruction"
            )
        ], width=12),
        dbc.Col([
            html.Label("sequence 1"),
            html.Br(),
            dcc.Textarea(
                id='fasta-text-area',
                placeholder="type sequence 1 here",
                style={'width': '100%', 'height': 300}
            ),
            html.Br(),
            dbc.Button(
                'Confirm sequence 1',
                id='btn-text-1',
                color="primary",
                n_clicks=0
            )
        ], width=6),
        dbc.Col([
            html.Label("sequence 2"),
            html.Br(),
            dcc.Textarea(
                id='fasta-text-area-two',
                placeholder="type your fasta sequence here",
                style={'width': '100%', 'height': 300}
            ),
            html.Br(),
            dbc.Button(
                'Confirm sequence 2',
                id='btn-text-2',
                color="primary",
                n_clicks=0
            )
        ], width=6)
    ]),
    html.Hr(),
    #Let’s say match is 1-10 mismatch is -1 to -10. Gap penalty is -1 to -10 as well
    dbc.Row([
        dbc.Col([
            html.Label("match score, default=5"),
            html.Br(),
            dcc.Slider(
                id="match-slider",
                min=1,
                max=10,
                step=1,
                value=5,
            ),
            html.Div(id="match-slider-val")
        ], width=4),
        dbc.Col([
            html.Label('mismatch score, default=5'),
            html.Br(),
            dcc.Slider(
                id='mismatch-slider',
                min=-10,
                max=-1,
                step=1,
                value=-5,
            ),
            html.Div(id="mismatch-slider-val")
        ], width=4),
        dbc.Col([
            html.Label('gap penalty, default=-5'),
            html.Br(),
            dcc.Slider(
                id='gap-penalty-slider',
                min=-10,
                max=-1,
                step=1,
                value=-5,
            ),
            html.Div(id="gap-penalty-slider-val"),
        ], width=4)
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select Matrix"),
            dcc.Dropdown(
                id='matrix-dropdown',
                options = [{'label': matrix, 'value': matrix} for matrix in matrix_names_arr
                            ],
                value="BLOSUM 62"
            )
        ], width=4),
        dbc.Col([
            html.Br(),
            dbc.Button(
                'Run Needleman Wunsch',
                id='submit-fastas',
                color="success",
                n_clicks=0
            )
        ], width=4)
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(#https://dash-bootstrap-components.opensource.faculty.ai/docs/components/progress/
                id="progress-bar"
            )
        ], width=12)
    ]),
    html.Hr(),
    dbc.Row([
        dbc.Col([
            html.H3("Needleman Results")
        ], width=12),
        dbc.Col([
            html.H3("Score"),
            html.P(
                id="score"
            )
        ], width=12),
        dbc.Col([
            html.H3("Alignments"),
            html.P(
                id="alignments"
            )
        ], width=12)
    ])
])


######################################################################################
