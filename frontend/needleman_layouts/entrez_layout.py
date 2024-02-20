'''
todo: accept user email here and store as a session variable for each user
'''


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
from dash_extensions import Download

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_stylesheets = [dbc.themes.BOOTSTRAP]

# for dropdown
'''
matrix_names_arr = ['PAM10','PAM30', 'PAM70', 'PAM250',
                    'BLOSUM62', 'BLOSUM45', 'BLOSUM50', 'BLOSUM80','BLOSUM90']
'''
matrix_names_arr = ['PAM30', 'PAM70', 'PAM250',
                    'BLOSUM62', 'BLOSUM45', 'BLOSUM50', 'BLOSUM80','BLOSUM90']
entrez_page = dbc.Container([
    #########################################################################
    # debug and store
    dbc.Row([
        dbc.Col([
            dcc.Store(id="user-email-store"),
            html.Div(
                "check user-email",
                id='check-user-email'
            ),
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
                "Input your email to use for NCBI Entrez API calls.  Do not abuse the service."
            ),
            html.Br(),
            html.Label("Email:"),
            dcc.Input(
                id="email-input",
                type="text",
                size="60"
            ),
            dbc.Button(
                'Store Email',
                id='btn-email',
                color="primary",
                n_clicks=0
            ),
        ])
    ]),
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
                value="NP_001239370.1",
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
                value="NP_715630.3",
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
            html.Label(
                "NCBI Search by Text Results appear here",
                style ={
                    "display": "none"
                }
            ),
            dcc.Loading(
                id='ncbi-textsearch-loading',
                children=[
                    html.Div(
                        # full fasta for accession num input
                        id="ncbi-textsearch-res"
                    )],
                    type="dot"
            ),
            html.Br(),
            dcc.Loading(
                id='query-translation-loading',
                children=[
                    html.Div(
                        id="query-translation-output"
                    )],
                    type="dot"
            )
        ], width=12)
    ]),
    html.Hr(),
    dbc.Row([
        dbc.Col([
            html.H3("3. Upload 2 fasta files or a multi-fasta file of 2 sequences")
        ], width=12),
        dbc.Col([
            dcc.Upload(
                id='upload-data',
                children=html.Div(
                    "Drag and Drop or \n"
                    "Select Fasta File 1 from Files"
                ),
                style={
                    'width': '95%',
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
            dcc.Upload(
                id='upload-data-two',
                children=html.Div([
                    "Drag and Drop or "
                    "Select Fasta File 2 from Files"
                ]),
                style={
                    'width': '95%',
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
    dbc.Row([
        dbc.Col([
            dbc.Button(
                'Get Genbank Info Seq1',
                id='btn-get-gb',
                color="primary",
                n_clicks=0
            ),
            html.Br(),
            dcc.Loading(
                id='genbank-info-loading',
                children=[
                    html.Div(
                        id="gb-output"
                    )],
                    type="circle"
            ),
            html.Br(),
            dbc.Button(
                'Get Genbank Info Seq2',
                id='btn-get-gb2',
                color="primary",
                n_clicks=0
            ),
            html.Br(),
            dcc.Loading(
                id='genbank2-info-loading',
                children=[
                    html.Div(
                        id="gb2-output"
                    )],
                    type="circle"
            ),
        ])
    ]),
    html.Hr(),
    #Let’s say match is 1-10 mismatch is -1 to -10. Gap penalty is -1 to -10 as well
    dbc.Row([
        dbc.Col([
            html.P(
                "Better alignments are usually obtained by penalizing gaps: " 
                "higher costs for opening a gap and lower costs for extending "
                "an existing gap. For amino acid sequences match scores are usually " 
                "encoded in matrices like PAM or BLOSUM. Defaults are -10 and -0.5."
            )
        ])
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("gap-open penalty"),
            html.Br(),
            dcc.Slider(
                id="gap-open-slider",
                min=-20,
                max=-0,
                step=0.5,
                value=-10,
            ),
            html.Div(id="gap-open-slider-val")
        ], width=4),
        dbc.Col([
            html.Label('gap-extend penalty'),
            html.Br(),
            dcc.Slider(
                id='gap-extend-slider',
                min=-5,
                max=0,
                step=0.5,
                value=-0.5,
            ),
            html.Div(id="gap-extend-slider-val")
        ], width=4),

    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Select Matrix"),
            dcc.Dropdown(
                id='matrix-dropdown',
                options = [{'label': matrix, 'value': matrix} for matrix in matrix_names_arr
                            ],
                value="BLOSUM62"
            )
        ], width=4),
        dbc.Col([
            html.Br(),
            dbc.Button(
                'Run Needleman-Wunsch (global)',
                id='run-needleman',
                color="success",
                n_clicks=0
            )
        ], width=4),
        dbc.Col([
            html.Br(),
            dbc.Button(
                'Run Smith-Waterman (local)',
                id='run-waterman',
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
            html.P(
                 "Bio.pairwise2 will return up to 1000 alignments."
            )
        ], width=12),
        dbc.Col([
            # todo: moved from app page to here
            dcc.Store(
                id="aligned-A"
            ),
            dcc.Loading(
                id='needleman-loading',
                children=[
                    html.Div(
                        id="needleman-output"
                    )],
                    type="cube"
            )
        ], width=12),
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Store(
                id="aligned-B"
            ),
            dcc.Loading(
                id='waterman-loading',
                children=[
                    html.Div(
                        id="waterman-output"
                    )],
                type="cube"
            )
        ], width=12),
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Store(
                id="aligned-fasta-store"
            ),
            Download(
                id="download-aligned-fasta"
            ),
            dbc.Button(
                "Make aligned Fasta for chart",
                id='btn-align-fasta',
                color="primary",
                n_clicks=0
            ),
            html.Br(),
            html.Br(),
            dbc.Button(
                "Download Aligned Fasta",
                id="download-btn",
                color="secondary",
                n_clicks=0
            )
        ], width=12),
        dbc.Col([
            dcc.Loading(
            id='aligned-fasta-loading',
            children=[
                html.Div(
                    id="aligned-fasta-output"
                )],
                type="cube"
            )
        ],width=12),
        dbc.Col([
            html.P("Needle Water context check"),
            html.Div(
                id="needle-water-ctx"
            )
        ], width=12),
    ])
])


######################################################################################
