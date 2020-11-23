#from jupyter_dash import JupyterDash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State # Load Data
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
import six.moves.urllib.request as urlreq
from six import PY3
from dash import no_update
import dash
import dash_bio as dashbio
import dash_html_components as html

#import ncbi search class
from frontend.ncbi.ncbi_search import get_last_updated, get_fasta_by_accession, get_full_GB_info, searchByTerm
#import needleman
from backend.bio_needleman import Needleman
import datetime
import json
import os

"""
When button pressed, get fastas and write a new file 
to folder sample_fastas, 
then include filepath as data in html output to Div element holding 
alignment chart

"""
#todo: cheating with global
file_path = ""


def register_visual_callbacks(app):
    """
    updates the event occurring on alignment chart
    :param app:
    :return:
    """

    @app.callback(
        Output('alignment-viewer-output', 'children'),
        [Input('my-alignment-viewer', 'eventDatum')]
    )
    def update_output(value):
        if value is None:
            return 'No data'
        return str(value)

    @app.callback(
        [Output('aligned-fasta', 'data')],
         [Input('btn-align-chart', 'n_clicks')],
        [State('accession-store-1', 'data'),
         State('accession-store-2', 'data'),
         State('aligned-A', 'data'),
         State('aligned-B', 'data')]
    )
    def update_aligned_fasta_store(n_clicks, acc1_json, acc2_json,
                                seqA_json, seqB_json):
        if n_clicks <=0:
            return no_update
        acc1 = json.loads(acc1_json)
        acc2 = json.loads(acc2_json)
        print(f'acc: {acc1}, {acc2}')
        # load aligned string sequences
        seqA = json.loads(seqA_json)
        seqB = json.loads(seqB_json)
        print(f'seqs from json:\n{seqA}\n{seqB}')

        fasta1_desc = get_fasta_by_accession(acc1, full_fasta=True)
        fasta2_desc = get_fasta_by_accession(acc2, full_fasta=True)

        #write a string of both fastas
        fasta_str = ">"
        fasta_str += fasta1_desc
        for char in range(len(seqA)):
            if char % 50 == 0 and char > 0:
                fasta_str += "\n"
        #terminate fasta 1
        fasta_str += "\n"
        fasta_str += ">"
        fasta_str += fasta2_desc
        for char in range(len(seqB)):
            if char % 50 == 0 and char > 0:
                fasta_str += "\n"

        # todo: need json?
        align_fasta_json = json.dumps(fasta_str)
        # try writing to file
        script_dir = os.path.dirname(__file__)
        rel_path = "temp.fasta"
        abs_file_path = os.path.join(script_dir, rel_path)

        with open(rel_path, "w") as outfile:
                outfile.write(fasta_str)
        outfile.close()

        return align_fasta_json

    """
    show the aligned fasta in dcc.Store()
    
    """
    @app.callback(
        Output('aligned-fasta-output', 'children'),
        [Input('aligned-fasta-store', 'data')]
    )
    def show_aligned_fasta(align_json):
        if align_json is None:
            return no_update
        align_str = json.loads(align_json)

        align_html = html.P(
            align_str,

        )
        return align_html

    """
    https://github.com/plotly/dash-bio/blob/master/tests/dashbio_demos/dash-alignment-chart/app.py
    seems to return json so try string, then json
    """
    @app.callback(
        [Output('my-alignment-viewer', 'data'),
         Output('my-alignment-viewer-json', 'data')],
        [Input('aligned-fasta-store', 'data')]
    )
    def update_chart_data(fasta_json):
        if not fasta_json:
            return no_update

        fasta_str = json.loads(fasta_json)
        return fasta_str, fasta_json
