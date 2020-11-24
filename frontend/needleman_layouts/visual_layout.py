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

'''
print(os.getcwd())
#todo: make a plots page
script_dir = os.path.dirname(__file__)
print(f'scriptdir: {script_dir}')
rel_path = 'human_v_mus.fasta'
abs_file_path = os.path.join(script_dir, rel_path)
print(f'abs_file_path: {abs_file_path}')

data = open(abs_file_path).read()
'''

DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Datasets
with open(os.path.join(DATAPATH, 'sample.fasta'), encoding='utf-8') as data_file:
    dataset1 = data_file.read()

with open(os.path.join(DATAPATH, 'p53.fasta'), encoding='utf-8') as data_file:
    dataset2 = data_file.read()

with open(os.path.join(DATAPATH, 'p53_clustalo.fasta'), encoding='utf-8') as data_file:
    dataset3 = data_file.read()

DATASETS = {
    'dataset1': dataset1,
    'dataset2': dataset2,
    'dataset3': dataset3,
}


plots_page = dbc.Container([
    dbc.Row([
        dbc.Col([
            dcc.Store(
                id="aligned-fasta-store"
            ),
            dcc.Store(
                id="aligned-fasta-store2"
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
            html.Br(),
            dcc.Dropdown(
                id='alignment-dropdown',
                options=[
                    {
                        'label': 'Sample.fasta',
                        'value': 'dataset1'
                    },
                    {
                        'label': 'P53.fasta naive',
                        'value': 'dataset2'
                    },
                    {
                        'label': 'P53.fasta aligned (ClustalW)',
                        'value': 'dataset3'
                    },
                ],
                value='dataset3',
            ),
            html.Br(),
            dcc.Upload(
                id='alignment-file-upload',
                className='control-upload',
                children=html.Div([
                    "Drag and drop multifasta files or select files."
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
            html.Br(),
            html.Div(
                id='desc-1'
            ),
            html.Div(
                id='vs'
            ),
            html.Div(
                id='desc-2'
            ),
            html.Br(),
            dcc.Loading(
                children=[
                    html.Div([
                        dash_bio.AlignmentChart(
                                id='alignment-chart',
                                height=725,
                                tilewidth=50,
                                data=dataset3
                        )
                    ])
                ]
            ),
            html.Div(
                id="alignment-viewer-output"
            )
        ], width=12)
    ])
])

