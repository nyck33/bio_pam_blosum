import dash_bio
print(dash_bio.__version__)

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

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_stylesheets = [dbc.theme.BOOTSTRAP]
df = px.data.tips()# Build App
app = Dash(external_stylesheets=external_stylesheets)

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

sidebar = html.Div(
    [
        html.H2("")
    ]
)

app.layout = dbc.Container(
    dbc.Row(
        dbc.Col()
    ))

@app.callback(
    Output('alignment-output', 'children'),
    [Input('testchart', 'eventDatum')]
)
def update_output(value):
    if value is None:
        return 'No data'
    return str(value)

if __name__ == "__main__":
    app.run_server(port=8888)