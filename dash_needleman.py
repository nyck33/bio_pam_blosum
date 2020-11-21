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

#import layouts
from frontend.needleman_layouts.entrez_layout import entrez_page #, register_entrez_callbacks
from frontend.needleman_layouts.intro_layout import needleman_intro
from frontend.needleman_layouts.visual_layout import plots_page

#import register_callbacks
from frontend.needleman_callbacks.entrez_callbacks import register_entrez_callbacks

#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
from Bio import Entrez
from Bio import SeqIO
# Entrez connection
Entrez.email = "nobutaka@gatech.edu"
# tool defaults to BioPython
#Entrez.tool = "getProteinFastas" #homologs later
db = 'protein'
#todo: learn to use this and adjust retmax
#paramEutils = {'usehistory': 'Y'}
handle = Entrez.einfo(db='protein')
# record is a dictionary record['DbInfo']['FieldList] shows
#record= Entrez.read(handle)
# see all available db's

#todo: backend import
from backend.pam import main, trace_back, build_matrics, load, parse_name, compare

#register stylesheet
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = Dash(external_stylesheets=external_stylesheets)
app.title = "Needleman Wunsch and NCBI"
#app.config['suppress_callback_exceptions'] = True

#######################################################################
#register callbacks
register_entrez_callbacks(app)
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
                dbc.NavLink("Project intro", href="/project-intro", id="intro-link"),
                dbc.NavLink("Entrez Search, Needleman-Wunsch, Smith-Waterman", href="/entrez-parameters", id="parameters-link"),
                dbc.NavLink("Dash Bio Demo", href="/plots", id="plots"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)
##############################################################
# todo: no fixed content
content = html.Div(id="page-content", style=CONTENT_STYLE)

app.layout=html.Div([dcc.Location(id="url"), sidebar, content])

##############################################################
##############################################################################
#callbacks for entire app
@app.callback(
    [Output('intro-link', "active"),
    Output('parameters-link', "active"),
    Output('plots', "active")],
    [Input("url", "pathname")],
)
def toggle_active_links(pathname):
    if pathname == "/" or pathname=="/project-intro":
        #Treat page intro as the homepage/index
        return True, False, False
    elif pathname == "/entrez-parameters":
        return False, True, False
    else:  #todo: need error pages when url pathname entry is not a match
        return False, False, True

@app.callback(Output("page-content", "children"),
              [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname in ["/", "/project-intro"]:
        return needleman_intro
    elif pathname == "/entrez-parameters":
        return entrez_page
    elif pathname == "/plots":
        return plots_page
    return dbc.Jumbotron([
    html.H1("404: Not found", className='text-danger'),
    html.Hr(),
    html.P(f'The pathname {pathname} was not recognized...')
])

#################################################################################
#register callbacks
#register_entrez_callbacks(app)
#register_intro_callbacks(app)
#####################################################################################
#entrez callbacks



#####################################################################################
if __name__=="__main__":
    app.run_server(debug=True, port=8080) #, dev_tools_ui=False, dev_tools_props_check=False)