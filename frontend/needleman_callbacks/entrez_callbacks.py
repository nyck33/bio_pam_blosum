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

#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
from backend.pam import main, trace_back, build_matrics, load, parse_name, compare


def register_entrez_callbacks(app):
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

    #################################################################
    #slider callbacks

    @app.callback(
        Output('match-slider-val', 'children'),
        [Input('match-slider', 'value')])
    def show_match_slider_val(match_value):
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
    # Update the dcc.Store from whichever input/button combo is pressed
    # todo: parse input of textbox 1 and 2 for multifasta?
    @app.callback(
        [Output('sequence-1-store', 'data'),
         Output('sequence-2-store', 'data'),
         Output('ctx-check', 'children'),
         Output('debug-out', 'children')],
        [Input('btn-acc-1', 'n_clicks'), #for top accession input
        Input('btn-acc-2', 'n_clicks'),
        Input('upload-data', 'contents'), #for dcc.Upload()
        Input('upload-data-two', 'contents'),
         Input('btn-text-1', 'n_clicks'),
         Input('btn-text-2', 'n_clicks')],
        [State('accession-input', 'value'), #top accession, 2 step process
        State('accession-input-two', 'value'),
        State('upload-data', 'filename'),
        State('upload-data', 'last_modified'),
         State('upload-data-two', 'filename'),
         State('upload-data-two', 'last_modified'),
         State('fasta-text-area', 'value'),
         State('fasta-text-area-two', 'value')])
    def update_store(btn_acc, btn_acc2, upload_contents, upload_contents2,
                     btn_txt, btn_txt2, acc_input, acc_input2, upload_name, upload_date,
                     upload_name2, upload_date2, txt_input1, txt_input2):
        if btn_acc<=0 and btn_acc2<=0 and upload_contents is None and upload_contents2 is None\
                and btn_txt<=0 and btn_txt2<=0:
            return no_update, no_update, no_update, no_update

        ctx = callback_context
        trigger = ctx.triggered[0]['prop_id'].split('.')[0]
        ctx_msg = json.dumps({
            'states': ctx.states,
            'triggered': ctx.triggered,
            'inputs': ctx.inputs
        }, indent=2)
        #accession input
        #store 1, store2, ctx, debug-out
        if trigger == "btn-acc-1":
            # calls search_fasta_by_accession()
            seq_str, seq_json = search_accession(acc_input)
            return seq_json, no_update, ctx_msg, seq_str
        elif trigger == "btn-acc-2":
            seq_str, seq_json = search_accession(acc_input2)
            return no_update, seq_json, ctx_msg, seq_str

        #upload inputs
        elif trigger == "upload-data":
            seq_json, html_output = process_upload(upload_contents,
                                                      upload_name,
                                                   upload_date)
            # make a string
            seq_str = json.loads(seq_json)
            return seq_json, no_update, ctx_msg, html_output
        elif trigger == "upload-data-two":
            seq_json, html_output = process_upload(upload_contents2,
                                                   upload_name2,
                                                   upload_date2)
            # make a string
            seq_str = json.loads(seq_json)
            return no_update, seq_json, ctx_msg, html_output
        #accepts seq part of fasta only, update store
        elif trigger == "btn-text-1":
            # store 1, store2, ctx, debug-out
            seq_json = json.dumps(txt_input1)
            return seq_json, no_update, ctx_msg, no_update
        elif trigger == "btn-text-2":
            seq_json = json.dumps(txt_input2)
            return no_update, seq_json, ctx_msg, no_update

    ##############################################
    # helpers
    def process_upload(contents, names, dates):
        if not contents and not names and not dates:
            return no_update
        # print(f"listContents: {type(contents)}\n{len(contents)}")

        for i in range(len(contents)):
            print(i, type(contents[i]), contents[i])
        # print(f"listFilename: {type(names)}\n{len(names)}")
        # print(f"listModified: {type(dates)}\n{len(dates)}")
        if contents is not None:
            seq_json, Pre_output = parse_contents(contents, names, dates)

            return seq_json, Pre_output


    # helpers for process upload()
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

        Pre_output = html.Div([
            html.H5(filename),
            html.H6(datetime.datetime.fromtimestamp(date)),

            html.Div([
                html.P(decoded)
            ]),
            html.Hr(),

            # debugging content
            html.Div('Raw Content'),
            html.Pre(contents[:],
                    style={
                        'whiteSpace': 'pre-wrap',
                        'wordBreak': 'break-all'
                    })

        ])

        return seq_json, Pre_output


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

    #helper 2
    def search_accession(acc1):
        #make calls to ncbi and get seq
        seq_str = get_fasta_by_accession(acc1)
        # store json
        seq_json = json.dumps(seq_str)
        # output for debug and store
        return seq_str, seq_json

    ##################################################################
    # for accession num search on Entrez and output json to store and show in Div
    @app.callback(
        [Output('ncbi-textsearch-res', 'children'),
        Output('query-translation-output', 'children')],
        Input('btn-search-ncbi', 'n_clicks'),
        State('term-search-input', 'value')
    )
    def get_accessions_arr(btn_search, searchterm):
        if btn_search<=0:
            return no_update, no_update
        accessions_arr, query_translation = searchByTerm(searchterm)

        acc_arr_Pre = html.Pre(accessions_arr[:],
                    style={
                        'whiteSpace': 'pre-wrap',
                        'wordBreak': 'break-all'
                    })

        return acc_arr_Pre, query_translation



################################################################################################################################
# run needle callbacks
    #might need:
    #State('match-slider', 'value') State('mismatch-slider', 'value)
    @app.callback([Output('matrix-output', 'children'),
                   Output('alignments-output', 'children')],
                  [Input('submit-fastas', 'n_clicks')],
                  [State('sequence-1-store', 'data'),
                   State('sequence-2-store', 'data'),
                   State('gap-penalty-slider', 'value'),
                   State('matrix-dropdown', 'value')])
    def run_needleman(run_btn, seq1_json, seq2_json, gap, matrix):
        if run_btn <=0:
            return no_update, no_update

        seq1 = json.loads(seq1_json)
        seq2 = json.loads(seq2_json)

        S1, S2, M = main(seq1, seq2, matrix, gap)

        align_dict = dict(zip(S1, S2))

        alignments_Pre = html.Pre(S1, S2, align_dict,
                             style={
                                 'whiteSpace': 'pre-wrap',
                                 'wordBreak': 'break-all'
                             }
                            )

        return M, alignments_Pre
'''
# run needle engine and output as check
        # todo: call main(param)
        S!, S2, F = main(params)
'''