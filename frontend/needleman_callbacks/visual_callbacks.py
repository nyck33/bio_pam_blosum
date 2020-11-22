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
from backend.pam import main, trace_back, build_matrics, load, parse_name, compare
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
def register_visual_callbacks(app):

    @app.callback(
        [Output('alignment-viewer-output', 'children')],
         [Input('btn-align-chart', 'n_clicks')],
        [State('accession-store-1', 'data'),
         State('accession-store-2', 'data')]
    )
    def output_alignment_chart(n_clicks, acc1_json, acc2_json):
        if n_clicks <=0:
            return no_update
        acc1 = json.loads(acc1_json)
        acc2 = json.loads(acc2_json)
        print(f'acc: {acc1}, {acc2}')

        fasta1 = get_fasta_by_accession(acc1, full_fasta=True)
        fasta2 = get_fasta_by_accession(acc2, full_fasta=True)

        # list of lines
        fasta1_lines = fasta1.splitlines()
        newline_arr = ["\n"]
        fasta2_lines = fasta2.splitlines()

        multifasta_lines = fasta1_lines + newline_arr + fasta2_lines

        #timestamp = datetime.datetime.now()
        #time = str(timestamp).replace(" ", "_").replace(":", "_").replace(".", "_")

        script_dir = os.path.dirname(__file__)
        rel_path = acc1 + \
                   "VS" + acc2 + ".fasta"
        abs_file_path = os.path.join(script_dir, rel_path)
        with open(rel_path, "a") as outfile:
            for line in multifasta_lines:
                outfile.write(line)

        output_html = html.Div([
            dashbio.AlignmentChart(
                #id='my-chart',
                data=abs_file_path,
                #colorscale='hydro',
                #conservationcolorscale='blackbody',
                tilewidth=50
            ),
            html.Div(id='alignment-viewer-output')
        ])

        return output_html



