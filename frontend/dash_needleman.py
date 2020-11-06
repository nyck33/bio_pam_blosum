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

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_stylesheets = [dbc.themes.BOOTSTRAP]
#df = px.data.tips()# Build App
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
        html.H3("PAM BLOSUM + NCBI API", className='display-4'),
        html.Hr(),
        html.P(
            "Choose page to display", className="lead"
        ),
        dbc.Nav(# todo: show res below params using callback and plots on plots
            [
                dbc.NavLink("Needleman Wunsch intro", href="/project-intro", id="intro-link"),
                dbc.NavLink("Entrez Gene and Parameters", href="/entrez-parameters", id="parameters-link"),
                dbc.NavLink("Plots of Results", href="/plots", id="plots"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)
fixed_content = html.Div([
    html.Div(
        "Write intro here", id="intro"
    ),
    html.Div(
        html.P("by Fei Xiang and Nobutaka Kim")
    )
])

needleman_intro = dbc.Container([
    html.H2("Talk about Needleman here and show some stuff", id="title-intro"),
    dbc.Row([
        dbc.Col([
            html.Div(
                "Some text here",
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

entrez_page = dbc.Container([
    html.H1("PAM/BLOSUM protein sequence scorer"),
    html.Div([]),
    html.P("Objective is to use the following matrices to score\
           similarities between 2 protein sequences "),

    dbc.Row([
        dbc.Col([
            html.H2("Enter Query Sequence"),
            dbc.FormGroup([
                dbc.Label("Enter FASTA sequence(s)")
            ])])
    ]
    )])
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



content = html.Div([
    html.Div(
        id="fixed-content", style=CONTENT_STYLE
    ),
    html.Div(id="page-content", style=CONTENT_STYLE)
])



app.layout=html.Div([dcc.Location(id="url"), sidebar, content])

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

@app.callback([Output("fixed-content", "children"),
               Output("page-content", "children")],
              [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname in ["/", "/project-intro"]:
        return fixed_content, needleman_intro
    elif pathname == "/entrez-parameters":
        return fixed_content, entrez_page
    elif pathname == "/plots":
        return fixed_content, plots_page
    return fixed_content, dbc.Jumbotron([
    html.H1("404: Not found", className='text-danger'),
    html.Hr(),
    html.P(f'The pathname {pathname} was not recognized...')
])

if __name__=="__main__":
    app.run_server(debug=True, port=8080)