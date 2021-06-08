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
    html.H2("Needleman-Wunsch, Smith-Waterman and Entrez",
            id="title-intro",
            style={
                'text-align': 'center'
            }),
    dbc.Row([
        dbc.Col([
            html.Div([
                html.H3("Objective"),
                html.P("To score protein alignments using the matrices from the NCBI website "
                       " downloaded from their FTP server.  This demonstrates construction of a "
                       " simple pipeline whose performance matches or exceeds those of EMBOSS if the BioPython "
                       " methods are used."),
                html.Br(),
                html.P(
                    "This is a work-in-progress to show my knowledge of Python, web programming, CI/CD, containerization,"
                    " and test automation as well as basic bioinformatics."
                    " My project repo (with to-do list) can be accessed at:"
                ),
                html.A(
                    "Link to bio_pam_blosum",
                    href="https://github.com/nyck33/bio_pam_blosum",
                    target="_blank"
                )

            ]),
        ], width=12),
        html.Br(),
        html.Br(),

        dbc.Row([
            dbc.Col([
                html.P(
                    "By Nobutaka Kim using Plotly Dash"
                    " and BioPython.  Testing with Selenium Web Driver, CI/CD with CircleCI and Docker."
                ),
                html.P(
                    "email: nobutaka@gatech.edu"
                )

            ], width=12)
        ])
    ])
])