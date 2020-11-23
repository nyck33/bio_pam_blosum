import dash_bio
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
import os
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm

print(os.getcwd())
#todo: make a plots page
script_dir = os.path.dirname(__file__)
print(f'scriptdir: {script_dir}')
rel_path = 'human_v_mus.fasta'
abs_file_path = os.path.join(script_dir, rel_path)
print(f'abs_file_path: {abs_file_path}')

data = open(abs_file_path).read()


plots_page = dbc.Container([
    dbc.Row([
        dbc.Col([
            dcc.Store(
                id="aligned-fasta-store"
            ),
            html.Div(
                id="aligned-fasta-output"
            ),
            html.H3("Alignment Chart"),
            html.P(
                "Input: alignment FASTA file "
                "click link for example: "
            ),
            html.A(
                "Alignment FASTA example",
                href="https://github.com/plotly/dash-bio/blob/master/tests/dashbio_demos/dash-alignment-chart/data/sample.fasta",
                target="_blank"
            ),
            dbc.Button(
                'Update Alignment Chart',
                id='btn-align-chart',
                color="primary",
                n_clicks=0
            ),
            html.Div([
                dash_bio.AlignmentChart(
                    id='my-alignment-viewer',
                    data = data
                ),
                html.Div(
                    id="alignment-viewer-output"
                )
            ]),
            html.Div([
                dash_bio.AlignmentChart(
                    id='my-alignment-viewer-json',
                    data = data #default
                ),
                html.Div(
                    id="alignment-viewer-output-json"
                )
            ])

        ], width=12)
    ])
])

