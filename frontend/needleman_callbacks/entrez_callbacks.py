import base64
import datetime
import io
import pandas as pd
import json
import textwrap
import os

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
from Bio import pairwise2, Align
from Bio.Align import substitution_matrices
# todo: download temp aligned fasta file to client
from dash_extensions import Download
import time
#import ncbi search class
from frontend.ncbi.ncbi_search import (
                                set_global_email_ncbi_search,
                                get_last_updated, get_fasta_by_accession,
                                get_full_GB_info, searchByTerm)
#import needleman
from backend.bio_needleman import (Needleman, matrix_load)

# todo: worker for Heroku
#from Bio import pairwise2

from rq import Queue, Retry
from worker import conn
from .aligner_utils import (local_align, global_align2)#, get_protein_alignment, make_single_seq)

def register_entrez_callbacks(app):
    #debug json stores
    #todo: use email debug fn to set global
    @app.callback(
        Output('check-user-email', 'children'),
        [Input('user-email-store', 'data')]
    )
    def show_user_email(usr_email_json):
        usr_email_str = json.loads(usr_email_json)
        #todo: this sets on ncbi_search.py
        set_global_email_ncbi_search(usr_email_str)
        return usr_email_str

    @app.callback(
        Output('check-store-1', 'children'),
        [Input('sequence-1-store', 'data')]
    )
    def show_seq1_json(seq1_json):
        seq1_str = json.loads(seq1_json)
        return seq1_str[:200]

    @app.callback(
        Output('check-store-2', 'children'),
        [Input('sequence-2-store', 'data')]
    )
    def show_seq1_json(seq2_json):
        seq2_str = json.loads(seq2_json)
        return seq2_str[:200]

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
    @app.callback(
        Output('user-email-store', 'data'),
        [Input('btn-email', 'n_clicks')],
        [State('email-input','value')])
    def store_email(btn_email, email_input):
        if btn_email <= 0:
            return no_update
        email_json = json.dumps(email_input)
        return email_json

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
        if contents is not None:
            seq_json, Pre_output, acc_num = parse_contents(contents, names, dates)

            return seq_json, Pre_output, acc_num


    # helpers for process upload()
    # for upload
    def parse_contents(contents, filename, date):
        contents = contents.replace('\n', '')
        content_type, content_string, = contents.split(',')

        decoded = decode_file_content(content_string)
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
        seq_json = json.dumps(seq)

        Pre_output = html.Div([
            html.H5(filename),
            html.H6(datetime.datetime.fromtimestamp(date)),
            html.P(f"accession: {acc_num}"),
            html.Div([
                html.P(decoded[:200])
            ]),
            html.Hr(),

            # debugging content
            html.Div('Raw Content'),
            html.Pre(contents[0:200] + '...',
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
        #todo: might need "UTF-8" in decode()
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
#########################################################
    #get GenBank info
    @app.callback(
        [Output('gb-output', 'children'),
         Output('gb2-output', 'children')],
        [Input('btn-get-gb', 'n_clicks'),
         Input('btn-get-gb2', 'n_clicks')],
        [State('accession-store-1', 'data'),
         State('accession-store-2', 'data')]
    )
    def get_gb_info(btn1, btn2, acc1_json, acc2_json):
        if btn1<=0 and btn2<=0:
            return no_update, no_update
        ctx = callback_context
        trigger = ctx.triggered[0]['prop_id'].split('.')[0]
        ctx_msg = json.dumps({
            'states': ctx.states,
            'triggered': ctx.triggered,
            'inputs': ctx.inputs
        }, indent=2)

        if trigger=="btn-get-gb":
            acc1 = json.loads(acc1_json)
            gb1_str = get_full_GB_info(acc1)

            gb_html = html.P(
                gb1_str,
                style={
                    'word-wrap': 'break-word'
                }
            )
            return gb_html, no_update
        elif trigger=="btn-get-gb2":
            acc2 = json.loads(acc2_json)
            gb2_str = get_full_GB_info(acc2)

            gb2_html = html.P(
                gb2_str,

            )
            return no_update, gb2_html


################################################################################################################################
# run needle and waterman
    @app.callback(
        [Output('aligned-A', 'data'),
         Output('aligned-B', 'data'),
        Output('needleman-output', 'children'),
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
            return no_update, no_update, no_update, no_update,\
                        no_update

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
        #needle = Needleman(seq1, seq2, mat_name, gap_open, gap_extend)
        #needle.load_matrix()
        #todo: call functions outside of class for heroku
        # make queue for RQ
        q = Queue(connection=conn)

        matrix = matrix_load(mat_name)
        align_str_arr = []
        gap_open = -10
        gap_extend = -0.5
        #aligner_obj = pairwise2.align()
        if trigger=='run-needleman':
            #needle.align_global()
            # get list of named tuples
            #alignments = needle.alignments
            #alignments = global_align_biop(seq1, seq2, matrix)
            '''
            pairwise2.align.globalds(seq1, seq2, matrix,
                                     gap_open, gap_extend,
                                     penalize_end_gaps=False,
                                     one_alignment_only=True)
            '''
            job = q.enqueue(global_align2, args=(seq1, seq2, matrix))

            count = 0
            while True:
                if job.result is not None or count > 9999:
                    break
                time.sleep(1)
                count+=1

            alignment = job.result
            # want to separate string into each protein and connecting lines
            seqA, connector, seqB = get_protein_alignment(alignment)
            #aligned_seq1 = make_single_seq(seq1, connector)
            #aligned_seq2 = make_single_seq(seq2, connector)

            # get the sequences from the first alignment and store
            # for chart
            #aligned_seq1 = alignment[0].seqA
            #aligned_seq2 = alignment[0].seqB
            #jsonify
            seqA_str = str(seqA)
            seqB_str = str(seqB)
            aligned1_json = json.dumps(str(seqA_str))
            aligned2_json = json.dumps(str(seqB_str))

            #alignment_html = format_output(alignment, "Needleman-Wunsch")

            alignments_html = format_for_aligner(alignment, "Needleman-Wunsch")
            return aligned1_json, aligned2_json, alignments_html, no_update,  ctx_msg

        if trigger=="run-waterman":
            #needle.align_local()
            #alignments = needle.alignments

            job = q.enqueue(local_align, args=(seq1, seq2, matrix))

            count = 0
            while True:
                if job.result is not None or count > 9999:
                    break
                time.sleep(1)
                count += 1
                #print(f'job.get_id(): {job.get_id()}, '
                 #     f'job.result:{job.result}')


            alignment = job.result
            # want to separate string into each protein and connecting lines
            seqA, connector, seqB = get_protein_alignment(alignment)
            #aligned_seq1 = make_single_seq(seq1, connector)
            #aligned_seq2 = make_single_seq(seq2, connector)

            # get the sequences from the first alignment and store
            # for chart
            #aligned_seq1 = alignment[0].seqA
            #aligned_seq2 = alignment[0].seqB
            #jsonify
            seqA_str = str(seqA)
            seqB_str = str(seqB)
            aligned1_json = json.dumps(str(seqA_str))
            aligned2_json = json.dumps(str(seqB_str))

            #alignment_html = format_output(alignment, "Needleman-Wunsch")

            alignments_html = format_for_aligner(alignment, "Smith-Waterman")

            #alignments_html = format_output(alignments, "Smith-Waterman")
            return aligned1_json, aligned2_json, no_update, alignments_html, ctx_msg




    def get_protein_alignment(alignment):
        """

        :param alignments: iterator of alignment objects
        :param idx: choose alignments[idx]
        :return: seqA and seqB alignments to make aligned FASTA
        """
        target = str(alignment)
        seq_arr = target.split()
        seqA, connectors, seqB = seq_arr

        return seqA, connectors, seqB

    #def format_for_aligner(seqA, connector, seqB, algo_name):
    def format_for_aligner(alignment, algo_name):
        """
        make an html div
        :param alignments:
        :param algo_name:
        :return:
        """
        # find max count

        alignments_html = \
            html.Div([
                html.H3(f"{algo_name} Alignments"),
                html.Div([
                    html.P(
                        str(alignment),
                        style={
                            'word-wrap': 'break-word'
                        }
                    )
                ])
            ])

        return alignments_html

        '''
        max_count = len(connector) // 50
        alignments_html = \
            html.Div([
                html.P(
                    seqA[count*50: (count*50)+50]
                    connector[count*50: (count*50)+50]
                    seqB[count * 50: (count*50) + 50]
                ) for count in range(0, max_count+1)
            ])
        '''

    def format_output(alignments, algo_name):
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
        alignment_type = ""
        if algo_name == "Needleman-Wunsch":
            alignment_type="Global"
        else:
            alignment_type="Local"
        alignments_html = \
            html.Div([
                html.H3(f"{algo_name} {alignment_type} Alignments"),
                html.Div([
                    html.P(
                        align_str_arr[i],
                        style={
                            'word-wrap': 'break-word'
                        }
                    ) for i in range(len(align_str_arr))
                ])
            ])


        return alignments_html

    @app.callback(
        [Output('aligned-A-output', 'children'),
        Output('aligned-B-output', 'children')],
        [Input('aligned-A', 'data'),
         Input('aligned-B', 'data')]
    )
    def show_aligned_seqs(a_json, b_json):
        a_str = json.loads(a_json)
        b_str = json.loads(b_json)

        alignA_html = html.Div([
            html.P(
                f'A aligned\n{a_str}',
                style={
                    'word-wrap': 'break-word'
                }
            )
        ])

        alignB_html = html.Div([
            html.P(
                f'B aligned\n{b_str}',
                style={
                    'word-wrap': 'break-word'
                }
            )
        ])

        return alignA_html, alignB_html

    @app.callback(
        #[Output("aligned-fasta-store", 'data'),
         #[Output("aligned-fasta-output", 'children'),
        [Output("aligned-fasta-store",'data'),
         Output('descrip-A-store', 'data'),
         Output('descrip-B-store', 'data')],
        [Input('btn-align-fasta', 'n_clicks')],
        [State("aligned-A", 'data'),
         State("aligned-B", "data")],
        [State('accession-store-1', 'data'),
         State('accession-store-2', 'data')]
    )
    def make_aligned_fasta(n_clicks, align1_json, align2_json,
                           acc1_json, acc2_json):
        if n_clicks <=0:
            return no_update, no_update, no_update
        acc1 = json.loads(acc1_json)
        acc2 = json.loads(acc2_json)
        #print(f'acc: {acc1}, {acc2}')
        # load aligned string sequences
        seqA = json.loads(align1_json)
        seqB = json.loads(align2_json)
        #print(f'seqs from json:\n{seqA}\n{seqB}')

        fasta1_desc = get_fasta_by_accession(acc1, full_fasta=True)
        fasta2_desc = get_fasta_by_accession(acc2, full_fasta=True)

        #jsonify descriptions
        fasta1_desc_json = json.dumps(str(fasta1_desc))
        fasta2_desc_json = json.dumps(str(fasta2_desc))

        #write a string of both fastas
        fasta_str = ""
        fasta_str += fasta1_desc
        char=0
        for char in range(len(seqA)):
            if char % 50 == 0 and char > 0:
                fasta_str += "\n"
            else:
                fasta_str += seqA[char]
        #terminate fasta 1
        fasta_str += "\n"
        #fasta_str += ">"
        char=0
        fasta_str += fasta2_desc
        for char in range(len(seqB)):
            if char % 50 == 0 and char > 0:
                fasta_str += "\n"
            else:
                fasta_str += seqB[char]

        # try writing to file
        #todo: problem here with absolute path,
        # get __file__ path for this file, go to pardir then change dir
        '''
        #full_path = "/home/nobu/Desktop/BioInformatics/bio_pam_blosum/frontend/needleman_layouts/data"
        script_dir = os.path.dirname(__file__)
        #https://stackoverflow.com/a/2860321/9481613
        par_dir = os.path.dirname(script_dir)
        data_dir = 'needleman_layouts/data'
        data_dir = os.path.join(par_dir, data_dir)
        #todo: name this file first line of fasta with date
        #rel_path = "aligned.fasta"
        
        #abs_file_path = os.path.join(full_path, rel_path)
        abs_file_path = os.path.join(data_dir, rel_path)
        with open(abs_file_path, "w") as outfile:
            outfile.write(fasta_str)
        '''
        # return (content=fasta_str, filename=rel_path)
        #jsonify
        aligned_fasta_json = json.dumps(fasta_str)

        #return aligned_fasta_json, fasta_str, fasta1_desc_json, \
        return aligned_fasta_json,\
               fasta1_desc_json, fasta2_desc_json

    @app.callback(
        Output("aligned-fasta-output", "children"),
        [Input("aligned-fasta-store", "data")]
    )
    def output_aligned_fasta(aligned_fasta_json):
        aligned_fasta_str = json.loads(aligned_fasta_json)
        # get output in form of html.Div
        alignment_fasta_div = format_final_output(aligned_fasta_str)
        return alignment_fasta_div

    def format_final_output(aligned_fasta_str):
        """
        :param alignments:
        :return: alignments_html
        """
        aligned_fasta_arr = aligned_fasta_str.splitlines()
        alignments_html = \
            html.Div([
                html.H3(f"Aligned Fasta"),
                html.Div([
                    html.P(
                        aligned_fasta_arr[i],
                        style={
                            'word-wrap': 'break-word'
                        }
                    ) for i in range(len(aligned_fasta_arr))
                ])
            ])
        return alignments_html

    @app.callback(
        Output("download-aligned-fasta", "data"),
        [Input("download-btn", "n_clicks")],
        [State("aligned-fasta-store",'data')]
    )
    def download_aligned_fasta(n_clicks, aligned_json):
        if n_clicks <= 0:
            return no_update
        elif aligned_json is None:
            return no_update

        fasta_str = json.loads(aligned_json)
        return dict(content=fasta_str, filename="aligned.fasta")
