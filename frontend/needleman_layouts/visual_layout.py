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


#todo: make a plots page
data = open('project/example_fasta_files/human_v_mus.fasta').read()
plots_page = html.Div([
    dashbio.AlignmentChart(
        id='testchart',
        data=data,
        colorscale='hydro',
        conservationcolorscale='blackbody',
        tilewidth=50
    ),
    html.Div(id='alignment-output')
])