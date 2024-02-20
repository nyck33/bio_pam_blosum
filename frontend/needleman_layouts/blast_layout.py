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

blast_layout = dbc.Container([
        dbc.Row([
        dbc.Col([
            html.H1("Run Blastp on Unknown sequence"),
            dcc.Store(
                id="blastp-res-store",
                storage_type='session'
            ),

        ], width=12),
        dbc.Col([
            html.Label("Upload blastp result XML if exists"),
            dcc.Upload(
                id='upload-blastp-xml',
                children=html.Div(
                    "Drag and Drop or \n"
                    "Select blastp XML file"
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
        ]),
        dbc.Col([
            html.Label("Enter unknown sequence"),
            html.Br(),
            dcc.Textarea(
                id='blast-text-area',
                placeholder="type sequence here",
                style={'width': '80%', 'height': 300}
            ),
            html.Br(),
            dbc.Button(
                "save query string to storage 1",
                id='save-query-button',
                color='primary',
                n_clicks=0
            ),
            html.Br(),
            html.Label("Enter XML file name"),
            html.Br(),
            dcc.Input(
                id='blast-file-name',
                type='text',
                size="50",
                placeholder="ex. my_blast"
            )
        ], width=12),
    ]),
    dbc.Row([
        dbc.Col([
            html.Label("Enter E-value"),
            dcc.Input(
                id='evalue-input',
                type='text',
                size="30",
                placeholder="1.0e-04 (0.00010) is an OK upper bound"
            ),
            html.Br(),
            dbc.Button(
                'Run blastp',
                id='btn-blastp',
                color="primary",
                n_clicks=0
            )
        ], width=12),
    ]),
    dbc.Row([
        dbc.Col([
            html.Label(
                "blastp results here",
                style ={
                    "display": "block"
                }
            ),
            dcc.Loading(
                id="blastp-len-alignments-loading",
                children=[
                    html.Div(
                        id='blastp-len'
                    )],
                    type="dot"
            ),
            dcc.Loading(
                id='blastp-loading',
                children=[
                    html.Div(
                        # full fasta for accession num input
                        id="blastp-res"
                    )],
                    type="dot"
            ),
        ], width=12)
    ]),
    dbc.Row([
        dbc.Col([
            html.H3("Filter Results"),
        ], width=12),
        dbc.Col([
            html.Label("Enter start idx"),
            dcc.Input(
                id='blastp-res-start',
                type='number',
                min=0, step=1,
                placeholder="first results are 0 to 10"
            )
        ], width=12),
        dbc.Col([
            html.Label('Enter ending idx'),
            dcc.Input(
                id='blastp-res-end',
                type='number',
                min=1, step=1,
                placeholder="first results are 0 to 10"
            )
        ], width=12),
        dbc.Col([
            dbc.Button(
                'Filter results',
                id='btn-filter-blastp',
                color="primary",
                n_clicks=0
            )
        ], width=12)
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Loading(
                id="filter-res-loading",
                children=[
                    html.Div(
                        id='filter-res'
                    )],
                type="dot"
            ),
        ], width=12)
    ])
])