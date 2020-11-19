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
from dash import Dash, no_update, exceptions, callback_context
from dash_table import DataTable
#from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State # Load Data
from Bio import SeqIO


def register_entrez_callbacks(app):
    @app.callback(Output('match-slider-val', 'children'),
                [Input('match-slider', 'value')])
    def show_slider_vals(value):
        return f'match score: {value}'

    @app.callback(Output('mismatch-slider-val', 'children'),
                  [Input('mismatch-slider', 'value')])
    def show_mismatch_slider_val(value):
        return f'mismatch score: {value}'

    @app.callback(Output('gap-penalty-slider-val', 'children'),
                  [Input('gap-penalty-slider', 'value')])
    def show_gap_pen_slider_val(value):
        return f'gap penalty: {value}'

    @app.callback([Output('score', 'children'),
                   Output('alignments', 'children')],
                  [Input('submit-fastas', 'n_clicks')],
                  [State('sequence-1-store', 'data'),
                  State('sequence-2-store', 'data')])

    def run_needleman(run_btn, seq1_json, seq2_json):
        pass



    # callbacks for entrez_page
    # callbacks for processing file upload processing file upload
    @app.callback([Output('sequence-1-store', 'data'),
                   Output('output-data-upload', 'children')],
                  [Input('upload-data', 'contents')],
                  [State('upload-data', 'filename'),
                   State('upload-data', 'last_modified')])
    def process_upload(list_of_contents, list_of_names, list_of_dates):
        if not list_of_contents and not list_of_names and not list_of_dates:
            return no_update
        print(f"listContents: {type(list_of_contents)}\n{list_of_contents[0]}")

        for i in range(len(list_of_contents)):
            print(i, type(list_of_contents[i]), list_of_contents[i])
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
        contents = contents.replace('\n', '')
        content_type, content_string, = contents.split(',')
        print(f'name:{type(filename)}\n{filename}\n')
        print(f'type:{type(content_type)}\n{content_type}\n')
        print(f'string:{type(content_string)}\n{content_string}\n')
        decoded = decode_file_content(content_string)
        print(f'decoded:\n{decoded}\n')
        #get substring from line 2
        decoded_arr = decoded.split()
        seq=""
        for line in decoded_arr:
            if len(line) >10 and line.upper() == line and line.isalpha():
                seq+=line

        assert seq.isalpha()
        print(seq)
        print(len(seq))

        return json.dumps(seq),\
                    html.Div([
                    html.H5(filename),
                    html.H6(datetime.datetime.fromtimestamp(date)),

                    html.Div([
                        html.P(decoded)
                    ]),
                    html.Hr(),

                    #debugging content
                    html.Div('Raw Content'),
                    html.Pre(contents[:] + '...', style={
                        'whiteSpace': 'pre-wrap',
                        'wordBreak': 'break-all'
                    })

                ])

    # annotations (param: param type) -> return type
    def decode_file_content(file: str) -> str:
        """
        Decode file from base64 to ascii
        Parameters
        :param file: str, required
        A base64 string
        :return:
        results: decoded str
        """
        return base64.b64decode(file.split(',')[-1].encode('ascii')).decode()

