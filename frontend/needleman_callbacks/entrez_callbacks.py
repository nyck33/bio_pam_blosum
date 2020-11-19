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

from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm


def register_entrez_callbacks(app):

    #################################################################
    #slider callbacks
    @app.callback(
        dash.dependencies.Output('slider-output-container', 'children'),
        [dash.dependencies.Input('my-slider', 'value')])
    def update_output(value):
        return 'You have selected "{}"'.format(value)


    @app.callback(
        Output('match-slider-val', 'children'),
        [Input('match-slider', 'value')])
    def show_slider_vals(match_value):
        val_string = f'match score: {match_value}'
        return val_string


    @app.callback(Output('mismatch-slider-val', 'children'),
                  [Input('mismatch-slider', 'value')])
    def show_mismatch_slider_val(mismatch_value):
        val_string = f'mismatch score: {mismatch_value}'
        return val_string


    @app.callback(Output('gap-penalty-slider-val', 'children'),
                  [Input('gap-penalty-slider', 'value')])
    def show_gap_pen_slider_val(gap_pen_value):
        val_string = f'gap penalty: {gap_pen_value}'
        return val_string

    ######################################################################
    # run needle callbacks
    @app.callback([Output('score', 'children'),
                   Output('alignments', 'children')],
                  [Input('submit-fastas', 'n_clicks')],
                  [State('sequence-1-store', 'data'),
                   State('sequence-2-store', 'data')])
    def run_needleman(run_btn, seq1_json, seq2_json):
        return no_update, no_update


    ######################################################################
    # dcc.Upload callbacks
    @app.callback([Output('sequence-1-store', 'data'),
                   Output('output-data-upload', 'children')],
                  [Input('upload-data', 'contents')],
                  [State('upload-data', 'filename'),
                   State('upload-data', 'last_modified')])
    def process_upload1(contents, names, dates):
        if not contents and not names and not dates:
            return no_update
        # print(f"listContents: {type(contents)}\n{len(contents)}")

        for i in range(len(contents)):
            print(i, type(contents[i]), contents[i])
        # print(f"listFilename: {type(names)}\n{len(names)}")
        # print(f"listModified: {type(dates)}\n{len(dates)}")
        if contents is not None:
            seq_json, children = parse_contents(contents, names, dates)

            return seq_json, children


    # helpers for functions
    # for upload
    def parse_contents(contents, filename, date):
        contents = contents.replace('\n', '')
        content_type, content_string, = contents.split(',')
        print(f'name:{type(filename)}\n{filename}\n')
        print(f'type:{type(content_type)}\n{content_type}\n')
        print(f'string:{type(content_string)}\n{content_string}\n')
        decoded = decode_file_content(content_string)
        print(f'decoded:\n{decoded}\n')
        # get substring from line 2
        decoded_arr = decoded.split()
        seq = ""
        for line in decoded_arr:
            if len(line) > 10 and line.upper() == line and line.isalpha():
                seq += line

        assert seq.isalpha()
        print(seq)
        print(len(seq))
        seq_json = json.dumps(seq)

        test_output = html.Div([
            html.H5(filename),
            html.H6(datetime.datetime.fromtimestamp(date)),

            html.Div([
                html.P(decoded)
            ]),
            html.Hr(),

            # debugging content
            html.Div('Raw Content'),
            html.Pre(contents[:] + '...', style={
                'whiteSpace': 'pre-wrap',
                'wordBreak': 'break-all'
            })

        ])

        return seq_json, test_output


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


    # callbacks for entrez_page
    # callbacks for processing file upload processing file upload
    @app.callback([Output('sequence-2-store', 'data'),
                   Output('output-data-upload-two', 'children')],
                  [Input('upload-data-two', 'contents')],
                  [State('upload-data-two', 'filename'),
                   State('upload-data-two', 'last_modified')])
    def process_upload2(contents, names, dates):
        if not contents and not names and not dates:
            return no_update
        # print(f"listContents: {type(contents)}\n{len(contents)}")

        for i in range(len(contents)):
            print(i, type(contents[i]), contents[i])
        # print(f"listFilename: {type(names)}\n{len(names)}")
        # print(f"listModified: {type(dates)}\n{len(dates)}")
        if contents is not None:
            seq_json, children = parse_contents(contents, names, dates)

            return seq_json, children


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
        # get substring from line 2
        decoded_arr = decoded.split()
        seq = ""
        for line in decoded_arr:
            if len(line) > 10 and line.upper() == line and line.isalpha():
                seq += line

        assert seq.isalpha()
        print(seq)
        print(len(seq))
        seq_json = json.dumps(seq)

        test_output = html.Div([
            html.H5(filename),
            html.H6(datetime.datetime.fromtimestamp(date)),

            html.Div([
                html.P(decoded)
            ]),
            html.Hr(),

            # debugging content
            html.Div('Raw Content'),
            html.Pre(contents[:] + '...', style={
                'whiteSpace': 'pre-wrap',
                'wordBreak': 'break-all'
            })

        ])

        return seq_json, test_output


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

    ##############################################################################
    #debug json stores
    @app.callback(
        Output('check-store-1', 'children'),
        [Input('sequence-1-store', 'data')]
    )
    def show_seq1_json(seq1_json):
        seq1_str = json.loads(seq1_json)
        return seq1_str

    @app.callback(
        Output('check-store-2', 'children'),
        [Input('sequence-2-store', 'data')]
    )
    def show_seq1_json(seq2_json):
        seq2_str = json.loads(seq2_json)
        return seq2_str

    ##################################################################
    # for accession num search on Entrez and output json to store and show in Div
    @app.callback(
        [Output('accession-fasta-1', 'children'),
         Output('sequence-1-store', 'data')],
        [Input('btn-acc-1', 'n_clicks')],
        [State('accession-input', 'value')]
    )
    def search_accession(n_clicks, acc):
        if n_clicks <= 0:
            return no_update, no_update

        seq_str = get_fasta_by_accession(acc)



