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


needleman_intro = dbc.Container([
    html.H2("Talk about Needleman here and show some stuff", id="title-intro"),
    dbc.Row([
        dbc.Col([
            html.Div(
                html.P("Objective is to use the two matrices to score\
                           similarities between 2 protein sequences "),
                id="needleman-explanation",
            ),
        ], width=4),
        dbc.Col([
            html.Div(
                "Diagrams here",
                id="needleman-diagrams"
            )
        ])
    ])
])