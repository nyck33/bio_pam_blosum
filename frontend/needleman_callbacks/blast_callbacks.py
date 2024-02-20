import base64
import datetime
import io
import pandas as pd
import json
import textwrap
import os
from copy import deepcopy

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

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
#import needleman
from backend.bio_needleman import Needleman

def register_blast_callbacks(app):

    @app.callback(
        [Output('blastp-res', 'children'),
         Output('blastp-res-store', 'data'),
         Output('blastp-len', 'children'),
         Output('blast-text-area', 'value')],
        [Input('btn-blastp', 'n_clicks'),
         Input('upload-blastp-xml', 'contents')],
        [State('blast-text-area', 'value'),
         State('blast-file-name', 'value'),
         State('evalue-input', 'value'),
         State('upload-blastp-xml', 'filename'),
         State('upload-blastp-xml', 'last_modified')])
    def run_blastp(n_clicks, xml_contents, seq, filename, evalue_txt,
                   xml_filename, xml_date):
        if n_clicks <=0 and xml_contents ==None:
            return no_update, no_update, no_update, no_update
        ctx = callback_context
        trigger = ctx.triggered[0]['prop_id'].split('.')[0]
        ctx_msg = json.dumps({
            'states': ctx.states,
            'triggered': ctx.triggered,
            'inputs': ctx.inputs
        }, indent=2)
        # default
        file_name = "my_blast.xml"
        # convert e-value
        if evalue_txt is None:
            evalue_txt = "1.0e-04"
        e_value_thresh = float(evalue_txt)
        if trigger == "upload-blastp-xml":
            print(os.getcwd())
            file_name = xml_filename

        elif trigger == "btn-blastp":

            #take out whitespace
            seq_arr = seq.split()
            cleaned_seq = ""
            for segment in seq_arr:
                cleaned_seq+=segment

            result_handle = NCBIWWW.qblast('blastp', 'nr',
                                           cleaned_seq)

            if filename:
                file_name = filename + ".xml"

            with open(file_name, "w") as out_handle:
                out_handle.write(result_handle.read())

            result_handle.close()

        result_handle = open(file_name)

        blast_records = NCBIXML.read(result_handle)

        alignments_dict_arr = []
        alignment_dict = {}
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_value_thresh:
                    alignment_dict['desc'] = alignment.title
                    alignment_dict['length'] = alignment.length
                    alignment_dict['e-value'] = hsp.expect
                    alignment_dict['query_seq'] = hsp.query[0:80]
                    alignment_dict['query_match'] = hsp.match[0:300]
                    alignment_dict['subject_seq'] = hsp.sbjct[0:]
            alignments_dict_arr.append(deepcopy(alignment_dict))
            alignment_dict.clear()
        #query string
        alignment_dict = alignments_dict_arr[0]
        query_str = alignment_dict['query_seq']

        #len
        len_alignments_arr = len(alignments_dict_arr)
        len_str = html.H3(f"There are {len_alignments_arr} hits")
        # store the alignment dict_arr
        blastp_res_json = json.dumps(alignments_dict_arr)

        #output Div with button for choice
        blast_res_html = html.Div([
            html.Div([
                html.Div([
                    html.Div([
                        html.Div([
                            html.H5(f"{k}"),
                            html.Br(),
                            dcc.Textarea(
                                value=f'{v}',
                                id=f'{v}-input',
                                style={'width': '100%',
                                       'height': 150}
                            ),
                            html.Br(),
                            dbc.Button(
                                "save sequence",
                                id=f'{v}-btn'
                            ),
                        ]) if k=='subject_seq' \
                        else \
                        html.Div([
                            html.H4(
                                f'{v}'
                            ),
                            html.Br(),

                        ]) if k=="desc" \
                        else
                        html.Div([
                            html.P(
                                f'{k}:\n {v}',
                                style={
                                    'word-wrap': 'break-word'
                                }
                            ),
                            html.Br()
                        ]) for k, v in alignment_dict.items()
                    ]) for alignment_dict in alignments_dict_arr
                ]),
            ])
        ])

        return blast_res_html, blastp_res_json, len_str, query_str
    """
    get button clicked, iterate the dicts_arr, 
    use a count a variable to get the desc (main)
    
    """
    @app.callback(
        Output('filter-res', 'children'),
        [Input('btn-filter-blastp', 'n_clicks')],
        [State('blastp-res-start', 'value'),
         State('blastp-res-end', 'value'),
         State('blastp-res-store', 'data')]
    )
    def filter_results(n_clicks, start, end, blastp_json):
        if n_clicks<=0:
            return no_update
        blastp_res_dict_arr = json.loads(blastp_json)

        filtered_dict_arr = blastp_res_dict_arr[start:end]

        filtered_html = html.Div([
            html.Div([
                html.Div([
                    html.Div([
                        html.Div([
                            html.H5(f"{k}"),
                            html.Br(),
                            dcc.Textarea(
                                value=f'{v}',
                                id=f'{v}-input',
                                style={'width': '100%',
                                       'height': 150}
                            ),
                            html.Br(),
                            dbc.Button(
                                "save sequence",
                                id=f'{v}-btn'
                            ),
                        ]) if k=='subject_seq' \
                        else \
                        html.Div([
                            html.H4(
                                f'{v}'
                            ),
                            html.Br(),

                        ]) if k=="desc" \
                        else
                        html.Div([
                            html.P(
                                f'{k}:\n {v}',
                                style={
                                    'word-wrap': 'break-word'
                                }
                            ),
                            html.Br(0)
                        ]) for k, v in alignment_dict.items()
                    ]) for alignment_dict in filtered_dict_arr
                ]),
            ])
        ])

        return filtered_html
    '''
    @app.callback(
        [Output('blast-seq1-store', 'data')],
        [Input('save-query-button', 'n_clicks'),
         Input(f'{v}-btn')]
    )
    
    
    @app.callback(
        Output('blast-seq2-store', 'data'),
        [Input(f'{v}-btn') for json.loads('blstp-res-store')[0][]]
    )
    '''