#import dash_bio
#print(dash_bio.__version__)
import base64
import datetime
import io
import pandas as pd
import json

import dash_bio as dashbio
import six.moves.urllib.request as urlreq
from six import PY3
import dash_html_components as html
import plotly.express as px
import dash
from dash import Dash, exceptions, no_update, callback_context
from dash_table import DataTable
#from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State  # Load Data
from Bio import pairwise2
from Bio import SeqIO
from Bio.Align import substitution_matrices

#import layouts
from frontend.needleman_layouts.intro_layout import needleman_intro
from frontend.needleman_layouts.entrez_layout import entrez_page, matrix_names_arr
from frontend.needleman_layouts.blast_layout import blast_layout
from frontend.needleman_layouts.visual_layout import plots_page

#import register_callbacks
from frontend.needleman_callbacks.entrez_callbacks import register_entrez_callbacks
from frontend.needleman_callbacks.blast_callbacks import register_blast_callbacks
from frontend.needleman_callbacks.visual_callbacks import register_visual_callbacks

#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
from Bio import Entrez
from Bio import SeqIO
# Entrez connection
#todo: putting all this in ncbi_search.py
#Entrez.email = ""
# tool defaults to BioPython
#Entrez.tool = "getProteinFastas" #homologs later
#db = 'protein'
#todo: learn to use this and adjust retmax
#paramEutils = {'usehistory': 'Y'}
#handle = Entrez.einfo(db='protein')
# record is a dictionary record['DbInfo']['FieldList] shows
#record= Entrez.read(handle)
# see all available db's

#todo: cheating with global
file_path = ""

#register stylesheet
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = Dash(external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)
app.title = "Needleman Wunsch and NCBI"
#app.config['suppress_callback_exceptions'] = True

#######################################################################
#register callbacks
register_entrez_callbacks(app)
register_visual_callbacks(app)
register_blast_callbacks(app)
###############################################################
SIDEBAR_STYLE={
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}
###############################################################
#layout components on every page
sidebar = html.Div(
    [
        html.H3("PAM BLOSUM + NCBI", className='display-4'),
        html.Hr(),
        html.P(
            "Choose page to display", className="lead"
        ),
        dbc.Nav(# todo: show res below params using callback and plots on plots
            [
                dbc.NavLink("Project Intro", href="/project-intro", id="intro-link"),
                dbc.NavLink("Entrez Search, Needleman-Wunsch, Smith-Waterman", href="/entrez-parameters", id="parameters-link"),
                dbc.NavLink("Blast", href="/blast", id="protein-blast"),
                dbc.NavLink("Alignment Chart", href="/plots", id="plots"),

            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)
##############################################################
# todo: no fixed content
content = dbc.Container([
    dbc.Row([
        dbc.Col([
            dcc.Store(
                id="blast-seq1-store"
            ),
            dcc.Store(
                id="blast-seq2-store"
            ),
            dcc.Store(
                id="accession-store-1"
            ),
            dcc.Store(
                id="accession-store-2"
            ),
            dcc.Store(
                id='descrip-A-store'
            ),
            dcc.Store(
                id="descrip-B-store"
            ),
            dcc.Store(
                id="aligned-A"
            ),
            dcc.Store(
                id="aligned-B"
            ),
            dcc.Store(
                id="aligned-fasta-store"
            ),
            html.Div(
                id="aligned-fasta-output",
                style={'display': 'none'}
            ),

            html.Div(
                id='aligned-A-output',
                style={'display': 'none'}
            ),
            html.Div(
                id='aligned-B-output',
                style={'display': 'none'}
            ),
            html.Div(
                id="page-content",
                style=CONTENT_STYLE
            )
        ], width=12)
    ])
])

app.layout=html.Div([dcc.Location(id="url"), sidebar, content])

##############################################################
##############################################################################
#callbacks for entire app
@app.callback(
    [Output('intro-link', "active"),
    Output('parameters-link', "active"),
    Output('protein-blast','active'),
    Output('plots', "active")],
    [Input("url", "pathname")],
)
def toggle_active_links(pathname):
    if pathname == "/" or pathname=="/project-intro":
        #Treat page intro as the homepage/index
        return True, False, False, False
    elif pathname == "/entrez-parameters":
        return False, True, False, False
    elif pathname == "/blast":
        return False, False, True, False
    elif pathname=="/plots":
        return False, False, False, True

@app.callback(Output("page-content", "children"),
              [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname in ["/", "/project-intro"]:
        return needleman_intro
    elif pathname == "/entrez-parameters":
        return entrez_page
    elif pathname=="/blast":
        return blast_layout
    elif pathname == "/plots":
        return plots_page
    return dbc.Jumbotron([
    html.H1("404: Not found", className='text-danger'),
    html.Hr(),
    html.P(f'The pathname {pathname} was not recognized...')
])

#####################################################################################

if __name__=="__main__":
    app.run_server(debug=True, port=8080) #, dev_tools_ui=False, dev_tools_props_check=False)
    #app.run_server