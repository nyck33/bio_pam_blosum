import base64
import datetime
import io
import pandas as pd

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
from Bio import SeqIO


def register_entrez_callbacks(app):
    # callbacks for entrez_page
    # callbacks for processing file uploadprocessing file upload
    @app.callback(Output('output-data-upload', 'children'),
                  [Input('upload-data', 'contents')],
                  [State('upload-data', 'filename'),
                   State('upload-data', 'last_modified')])
    def process_upload(list_of_contents, list_of_names, list_of_dates):
        print(f"listContents: {type(list_of_contents)}\n{list_of_contents[0]}")
        print(f"listFilename: {type(list_of_names)}\n{list_of_names[0]}")
        print(f"listModified: {type(list_of_dates)}\n{list_of_dates[0]}")
        if list_of_contents is not None:
            children=[
                parse_contents(c,n,d) for c,n,d in
                zip(list_of_contents, list_of_names, list_of_dates)
            ]
            return children

    # helpers for functions
    # for uplaod
    def parse_contents(contents, filename, date):
        content_type, content_string, = contents.split(',')
        print(f'name:{type(filename)}\n{filename}\n')
        print(f'type:{type(content_type)}\n{content_type}\n')
        print(f'string:{type(content_string)}\n{content_string}\n')
        decoded = base64.b64decode(content_string)
        print(f'decoded:{decoded}\n')
        try:
            
            seq = SeqIO.read(decoded, "fasta")
            print(f'\n\nseq:{seq}')
        except Exception as e:
            print(e)
            return html.Div([
                'Error processing file upload'
            ])

        return html.Div([
            html.H5(filename),
            html.H6(datetime.datetime.fromtimestamp(date)),

            html.Div(
                f'{decoded}'
            ),
            html.Hr(),

            #debugging content
            html.Div('Raw Content'),
            html.Pre(contents[0:100] + '...', style={
                'whiteSpace': 'pre-wrap',
                'wordBreak': 'break-all'
            })

        ])

