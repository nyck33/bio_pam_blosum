import base64
import datetime
import io
import pandas as pd
import json
import textwrap

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
from Bio import pairwise2
from Bio.Align import substitution_matrices

#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
from backend.pam import main, trace_back, build_matrics, load, parse_name, compare
#import needleman
from backend.bio_needleman import Needleman

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
        Output('gap-open-slider-val', 'children'),
        [Input('gap-open-slider', 'value')])
    def show_gap_open_slider_val(gap_open_val):
        val_string = f'match score: {gap_open_val}'
        return val_string


    @app.callback(Output('gap-extend-slider-val', 'children'),
                  [Input('gap-extend-slider', 'value')])
    def show_gap_extend_slider_val(gap_extend_val):
        val_string = f'mismatch score: {gap_extend_val}'
        return val_string

    ######################################################################
    # Update the dcc.Store from whichever input/button combo is pressed
    # todo: parse input of textbox 1 and 2 for multifasta?
    @app.callback(
        [Output('sequence-1-store', 'data'),
         Output('sequence-2-store', 'data'),
         Output('accession-store-1', 'data'),
         Output('accession-store-2', 'data'),
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
            return no_update, no_update, no_update, no_update, no_update, no_update

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
            #jsonify acc_input
            acc_json = json.dumps(acc_input)
            return seq_json, no_update, acc_json, no_update, \
                    ctx_msg, seq_str
        elif trigger == "btn-acc-2":
            seq_str, seq_json = search_accession(acc_input2)
            # jsonify acc_input
            acc2_json = json.dumps(acc_input2)
            return no_update, seq_json, no_update, acc2_json, \
                    ctx_msg, seq_str

        #upload inputs
        elif trigger == "upload-data":
            seq_json, html_output, acc_num = process_upload(upload_contents,
                                                      upload_name,
                                                   upload_date)
            # make a string
            seq_str = json.loads(seq_json)
            # jsonify acc_num
            acc_json = json.dumps(acc_num)
            return seq_json, no_update, acc_json, no_update, \
                    ctx_msg, html_output
        elif trigger == "upload-data-two":
            seq_json, html_output, acc_num = process_upload(upload_contents2,
                                                   upload_name2,
                                                   upload_date2)
            # make a string
            seq_str = json.loads(seq_json)
            # jsonify acc_num
            acc_json = json.dumps(acc_num)

            return no_update, seq_json, no_update, acc_json,\
                   ctx_msg, html_output
        #accepts seq part of fasta only, update store
        elif trigger == "btn-text-1":
            # store 1, store2, ctx, debug-out
            seq_json = json.dumps(txt_input1)
            return seq_json, no_update, no_update, no_update, \
                   ctx_msg, no_update
        elif trigger == "btn-text-2":
            seq_json = json.dumps(txt_input2)
            return no_update, seq_json, no_update, no_update,\
                   ctx_msg, no_update

    ##############################################
    # helpers
    def process_upload(contents, names, dates):
        if not contents and not names and not dates:
            return no_update
        # print(f"listContents: {type(contents)}\n{len(contents)}")


        if contents is not None:
            seq_json, Pre_output, acc_num = parse_contents(contents, names, dates)

            return seq_json, Pre_output, acc_num


    # helpers for process upload()
    # for upload
    def parse_contents(contents, filename, date):
        contents = contents.replace('\n', '')
        content_type, content_string, = contents.split(',')
        #print(f'name:{type(filename)}\n{filename}\n')
        #print(f'type:{type(content_type)}\n{content_type}\n')
        #print(f'string:{type(content_string)}\n{content_string}\n')
        decoded = decode_file_content(content_string)
        #print(f'decoded:\n{decoded}\n')
        # get substring from line 2
        decoded_arr = decoded.split()
        seq = ""
        acc_num = ""
        for line in decoded_arr:
            if line[0] == ">":
                acc_num = line.split()[0]
                acc_num = acc_num.replace(">", "")
                print(f'acc_num: {acc_num}')
            if len(line) > 10 and line.upper() == line and line.isalpha():
                seq += line

        assert seq.isalpha()
        #print(seq)
        #print(len(seq))
        seq_json = json.dumps(seq)

        Pre_output = html.Div([
            html.H5(filename),
            html.H6(datetime.datetime.fromtimestamp(date)),
            html.P(f"accession: {acc_num}"),
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

        return seq_json, Pre_output, acc_num


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
    def get_accessions_from_text(btn_search, searchterm):
        """

        :param btn_search:
        :param searchterm:
        :return: html Pre of array of accessions, str query translation
        """
        if btn_search<=0:
            return no_update, no_update
        accessions_arr, query_translation = searchByTerm(searchterm)

        #return A links, click for full GB info
        acc_arr_links = html.Div(
                    children=[
                    html.P(
                    html.A(accessions_arr[i],
                    n_clicks=0,
                    href=f"https://www.ncbi.nlm.nih.gov/search/all/?term={accessions_arr[i]}",
                    target="_blank",
                    id=f'link-{accessions_arr[i]}'))
                        for i in range(len(accessions_arr))])

        return acc_arr_links, query_translation

################################################################################################################################
# run needle and waterman
    @app.callback(
        [Output('needleman-output', 'children'),
         Output('needleman-output-2', 'children'),
       Output('waterman-output', 'children'),
       Output('needle-water-ctx', 'children')],
      [Input('run-needleman', 'n_clicks'),
       Input('run-waterman', 'n_clicks')],
      [State('sequence-1-store', 'data'),
       State('sequence-2-store', 'data'),
       State('gap-open-slider', 'value'),
       State('gap-extend-slider', 'value'),
       State('matrix-dropdown', 'value')])
    def run_needleman(btn_needle, btn_waterman, seq1_json, seq2_json,
                        gap_open, gap_extend, mat_name):
        if btn_needle <=0 and btn_waterman <=0:
            return no_update, no_update, no_update

        ctx = callback_context
        trigger = ctx.triggered[0]['prop_id'].split('.')[0]
        ctx_msg = json.dumps({
            'states': ctx.states,
            'triggered': ctx.triggered,
            'inputs': ctx.inputs
        }, indent=2)

        seq1 = json.loads(seq1_json)
        seq2 = json.loads(seq2_json)

        # instantiate class
        needle = Needleman(seq1, seq2, mat_name, gap_open, gap_extend)
        needle.load_matrix()
        align_str_arr = []
        if trigger=='run-needleman':
            needle.align_global()
            # get list of named tuples
            alignments = needle.alignments
            alignments_html = format_output(alignments)
            new_html = reformat_output(alignments)
            return alignments_html, no_update, ctx_msg

        if trigger=="run-waterman":
            needle.align_local()
            alignments = needle.alignments
            alignments_html = format_output(alignments)
            return no_update, alignments_html, ctx_msg

    def format_output(alignments):
        """

        :param alignments:
        :return: alignments_html
        """
        # format output results in string
        align_str_arr = []
        for i in range(len(alignments)):
            align_str = pairwise2.format_alignment(*alignments[i])
            align_str_arr.append(align_str)
        # reformat

        alignments_html = html.Div([
            html.P(
                align_str_arr[i],
                style={
                    'word-wrap': 'break-word'
                }
            ) for i in range(len(align_str_arr))
        ])

        return alignments_html

    def reformat_output(alignments_arr):
        """
        Just return one
        :param alignments_arr:
        :return:
        """
        alignment = alignments_arr[0]
        seqA = str(alignment.seqA)
        seqB = str(alignment.seqB)

        # break into 50 characters into array
        lenA = len(seqA)
        lenB = len(seqB)
        '''
        arrA = split seqA by 50
        arrB
        alignments_html = html.Div([
            html.P(
                f'{arrA[i]}\n{arrB[i]}',
                style={
                    'word-wrap': 'break-word'
                }
            ) for i in range(len(align_str_arr))
        ])
        '''
        pass