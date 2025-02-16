# Dash App Prep: imports, initializing, dataset!

import os
import base64
import io
from dash import Dash, html, dcc, Input, Output, ctx, dash_table, no_update, callback
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
import pandas as pd
import plotly.figure_factory as ff
import numpy as np
import math
from Bio import AlignIO, Phylo, SeqIO
import dash_cytoscape as cyto
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

import time
import argparse
import logging
import multiprocessing
import subprocess
import pickle
import zipfile
import traceback
from itertools import groupby

RUNNING_IN_DOCKER = os.getenv('ALGAEORTHO_DOCKER', 'False').lower() == 'true'
NO_CLUSTALO = os.getenv('ALGAEORTHO_NOCLUSTALO', 'False').lower() == 'true'

SAVE_FIGURES = os.getenv('SAVE_FIGURES', 'True').lower() == 'true'
FIGURE_WIDTH = int(os.getenv('FIGURE_WIDTH', '3000'))
FIGURE_HEIGHT = int(os.getenv('FIGURE_HEIGHT', '2000'))

# if the app is running in a docker container, the path is different
# note: there is a leading dot in the local path
data_path = "./data"
if RUNNING_IN_DOCKER:
    data_path = "/data"

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.debug("Starting AlgaeOrtho App...")

# Original Input File will be the results from SonicParanoid
#

# globals
df_to_merge = pd.DataFrame()
df_to_merge_loaded = False
selected_dl_file = "./hs2.pim.txt"
selected_heatmap_file = data_path + "/heatmap.png"

def load_shapes():
    # Start incorporating backend:
    logging.debug("loading ortholog_groups.tsv...")
    orthodf = pd.read_csv("ortholog_groups.tsv", sep="\t")
    horizontal = orthodf
    shape_num = horizontal.shape[1]
    horizontal.columns[4:shape_num].values[0:] = ([entry[:-9] for entry in horizontal.columns[4:shape_num]])
    pd.concat([horizontal.iloc[:, 0:3], horizontal.iloc[:, 4]], axis=1)

    d = {}
    for i in range(4, shape_num):
        logging.debug("loading shapes... " + str(i+1) + "/" + str(shape_num))
        [entry[:-9] for entry in horizontal.columns[4:shape_num]]
        d[horizontal.keys()[i]] = pd.DataFrame(pd.concat([horizontal.iloc[:, 0:3], horizontal.iloc[:, i]], axis=1))
    return d

def load_dataframes(d):
    d2 = dict.fromkeys(d.keys())
    newdf = pd.DataFrame()
    count = 0
    for var in d.keys():
        logging.debug("loading dataframes... " + str(count+1) + "/" + str(len(d.keys())))
        frame = d[var]
        tester_df = frame.replace("*", np.NaN)
        df1 = (tester_df.assign(category=tester_df.iloc[:, 3].str.split(','))
            .explode('category')
            .reset_index(drop=True))
        df1 = pd.concat([df1, df1['category'].str.split('|', expand=True)], axis=1)
        df2 = df1.iloc[:, 0:8]
        df2.columns = ["group_id", "group_size", "sp_in_grp", "full_entry", "category", "jgi", "species", "proteinId"]
        newdf = pd.concat([newdf, df2], axis=0)
        d2[var] = df2
        count += 1
    return newdf

def load_df_to_merge():
    d = load_shapes()
    newdf = load_dataframes(d)
    concise_df = newdf[newdf['proteinId'].notna()]
    #df_to_merge = concise_df[["group_id", "species", "proteinId"]]
    return concise_df[["group_id", "species", "proteinId"]]

def generate_tree_from_alignment(alignment, root_tree=True):
    calculator = DistanceCalculator('identity', )
    #distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator, "nj")
    tree = constructor.build_tree(alignment)
    if root_tree:
        tree.root_at_midpoint()
    return tree

load_figure_template('COSMO')
layout = {'name': 'preset'}
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
stylesheet = [
    {
        'selector': '.nonterminal',
        'style': {
            'label': 'data(confidence)',
            'background-opacity': 0,
            "text-halign": "left",
            "text-valign": "top",
        }
    },
    {
        'selector': '.support',
        'style': {'background-opacity': 0}
    },
    {
        'selector': 'edge',
        'style': {
            "source-endpoint": "inside-to-node",
            "target-endpoint": "inside-to-node",
        }
    },
    {
        'selector': '.terminal',
        'style': {
            'label': 'data(name)',
            'width': 10,
            'height': 10,
            "text-valign": "center",
            "text-halign": "right",
            #'background-color': '#222222',
            'backgroundColor': '#222222'
        }
    }
]

with open("lyco_ochr.clu", "r") as aln:
    logging.debug("loading lyco_ochr.clu...")
    alignment = AlignIO.read(aln, "clustal")

tree = generate_tree_from_alignment(alignment)
newick_string = str(tree)

# calculator = DistanceCalculator('identity', )
# distance_matrix = calculator.get_distance(alignment)
# constructor = DistanceTreeConstructor(calculator)
# tree = constructor.build_tree(alignment)
# newick_string = str(tree)


def generate_elements(tree, xlen=30, ylen=30, grabbable=False):
    logging.debug("generate_elements - xlen:" + str(xlen) + 
                  " ylen:" + str(ylen) + 
                  " grabbable:" + str(grabbable))

    def get_col_positions(tree, column_width=80):
        taxa = tree.get_terminals()

        # Some constants for the drawing calculations
        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width = column_width - max_label_width - 1

        """Create a mapping of each clade to its column position."""
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
            # Potential drawing overflow due to rounding -- 1 char per tree layer
        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
        cols_per_branch_unit = ((drawing_width - fudge_margin) /
                                float(max(depths.values())))
        return dict((clade, int(blen * cols_per_branch_unit + 1.0))
                    for clade, blen in depths.items())

    def get_row_positions(tree):
        taxa = tree.get_terminals()
        positions = dict((taxon, 2 * idx) for idx, taxon in enumerate(taxa))

        def calc_row(clade):
            for subclade in clade:
                # logging.debug("TEST: subclade" + subclade)
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = ((positions[clade.clades[0]] +
                                 positions[clade.clades[-1]]) // 2)

        # logging.debug("TEST: tree.root" + tree.root)
        calc_row(tree.root)
        return positions

    def add_to_elements(clade, clade_id):
        children = clade.clades

        pos_x = col_positions[clade] * xlen
        pos_y = row_positions[clade] * ylen

        cy_source = {
            "data": {"id": clade_id},
            'position': {'x': pos_x, 'y': pos_y},
            'classes': 'nonterminal',
            'grabbable': grabbable
        }
        nodes.append(cy_source)

        if clade.is_terminal():
            cy_source['data']['name'] = clade.name
            cy_source['classes'] = 'terminal'

        for n, child in enumerate(children):
            # The "support" node is on the same column as the parent clade,
            # and on the same row as the child clade. It is used to create the
            # 90 degree angle between the parent and the children.
            # Edge config: parent -> support -> child

            support_id = clade_id + 's' + str(n)
            child_id = clade_id + 'c' + str(n)
            pos_y_child = row_positions[child] * ylen

            cy_support_node = {
                'data': {'id': support_id},
                'position': {'x': pos_x, 'y': pos_y_child},
                'grabbable': grabbable,
                'classes': 'support'
            }

            cy_support_edge = {
                'data': {
                    'source': clade_id,
                    'target': support_id,
                    'sourceCladeId': clade_id
                },
            }

            cy_edge = {
                'data': {
                    'source': support_id,
                    'target': child_id,
                    'length': clade.branch_length,
                    'sourceCladeId': clade_id
                },
            }

            if clade.confidence and clade.confidence.value:
                cy_source['data']['confidence'] = clade.confidence.value

            nodes.append(cy_support_node)
            edges.extend([cy_support_edge, cy_edge])

            add_to_elements(child, child_id)

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)

    nodes = []
    edges = []

    add_to_elements(tree.clade, 'r')

    return nodes, edges


temp_node, temp_edge = generate_elements(tree)
elements = temp_node + temp_edge

logging.debug("loading Dash app...")
app = Dash(external_stylesheets=[dbc.themes.COSMO],suppress_callback_exceptions=True)

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "24rem",
    "padding": "2rem 1rem",
    #"background-color": "#7393B3",
    "backgroundColor": "#7393B3",
}

# set clustalo check box to disabled if NO_CLUSTALO is set
clustalo_disabled = False
clustalo_value = 'run_clustalo'
clustalo_values = [clustalo_value]
if NO_CLUSTALO:
    clustalo_disabled = True
    clustalo_values = []

sidebar_divs = [
        html.H2("Select Input"),
        html.Hr(),
        html.H3(children='The following two buttons select pre-loaded data'),
        html.Button('Chlorophyta bZIP1 case study 1', id='button1', n_clicks=0),
        html.Button('Ochrophyta LCYB case study 2', id='button2', n_clicks=0),

        html.Div(id='output-data-upload'),
        html.Div(id='container-button-timestamp'),

        html.Hr(),

        html.H3(children='Or, upload your own!'),

        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Upload ',
                html.A('your protein query')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px',
                'cursor': 'pointer'
            },
            # Allow multiple files to be uploaded
            multiple=False
        ),
        html.Hr(),

        html.H3(children='Options'),
        
        dcc.Checklist(
            options=[
                {
                    'label': 'If checked, trees will be displayed rooted.', 
                    'value': 'root_tree', 
                    'disabled': False
                 }
            ],
            value=['root_tree'],
            id="root_tree_option",
            inline=True
        ),
        dcc.Checklist(
            options=[
                {
                    'label': 'If checked, Clustal Omega will run to visualize fasta on user uploads (can be slow). If unchecked, fasta will be produced without visualization.', 
                    'value': clustalo_value, 
                    'disabled': clustalo_disabled
                 }
            ],
            value=clustalo_values,
            id="run_clustalo_option",
            inline=True
        ),
        html.Hr(),

        html.Button("Download results", id="btn-download-txt", n_clicks=0),
        dcc.Download(id="download-txt"),
    ]

if SAVE_FIGURES:
    sidebar_divs = sidebar_divs + [
        html.Hr(),

        html.Button("Download hi-res heatmap image", id="btn-download-heatmap", n_clicks=0),
        dcc.Download(id="download-heatmap"),

        html.Button("Download hi-res tree image", id="btn-download-tree-png", n_clicks=0),
    ]

sidebar = html.Div(
    sidebar_divs,
    style=SIDEBAR_STYLE,
)

# left column width out of a 12 column grid
LEFT_COL_WIDTH = 4

app.layout = html.Div(children=[
    # dbc.Col([
        dbc.Row([
            dbc.Col(sidebar, width=LEFT_COL_WIDTH),
            dbc.Col(html.H1(children='AlgaeOrtho'),width=(12-LEFT_COL_WIDTH)),
            # dbc.Col(html.H3(children='Heatmap (scroll down for tree!)'))
        ]),
        dbc.Row([
            dbc.Col(width=LEFT_COL_WIDTH),
            dbc.Col(dcc.Loading(
                id="loading-1",
                type="default",
                children=html.Div(
                    id="loading-output-1",
                    style={'width': '100%', 'height': '5vh'}
                )
            ),width=(12-LEFT_COL_WIDTH))
        ]),
        dbc.Row([
            dbc.Col(width=LEFT_COL_WIDTH),
            dbc.Col([html.H3(children='Heatmap (scroll down for tree!)'),
                    dcc.Graph(id='graph', style={'width': '100%', 'height': '90vh'})],width=(12-LEFT_COL_WIDTH))
        ]),
        # dbc.Row([
        #     dbc.Col(width=LEFT_COL_WIDTH),
        #     dbc.Col([html.Button("Download Results", id="btn-download-txt", n_clicks=0),
        #              dcc.Download(id="download-txt")])
        # ]),

        dbc.Row([
            dbc.Col(width=LEFT_COL_WIDTH),
            # dbc.Col(html.H3(children='Tree'),),
            dbc.Col([html.H3(children='Tree'),
                cyto.Cytoscape(
                    id='cytoscape-usage-phylogeny',
                    elements=elements,
                    stylesheet=stylesheet,
                    layout=layout,
                    style={
                        'height': '95vh',
                        'width': '100%'
                    }
                )],width=(12-LEFT_COL_WIDTH))
        ]),

        dbc.Row([
            dbc.Col(width=LEFT_COL_WIDTH),
            # dbc.Col(html.H3(children='Tree Text (newick)')),
            dbc.Col([html.H3(children='Tree Text (newick)'),
                html.Div(id='textarea-example-output', style={'whiteSpace': 'break-spaces'})],width=(12-LEFT_COL_WIDTH))
        ])
    # ])
])


@callback(
    Output("cytoscape-usage-phylogeny", "generateImage"),
    [
        Input("btn-download-tree-png", "n_clicks"),
    ])
def download_cyto_image(n_clicks):
    logging.debug("download_cyto_image - n_clicks:" + str(n_clicks))
    # 'store': Stores the image data in 'imageData' !only jpg/png are supported
    # 'download'`: Downloads the image as a file with all data handling
    # 'both'`: Stores image data and downloads image as file.
    ftype = 'png'
    action = 'store'

    if ctx.triggered:
        if ctx.triggered_id != "tabs-image-export":
            action = "download"
            ftype = ctx.triggered_id.split("-")[-1]

    return {
        'type': ftype,
        'action': action
    }

@app.callback(
    [Output("download-heatmap", "data")],
    [Input("btn-download-heatmap", "n_clicks")],
    prevent_initial_call=True
)
def download_heatmap_image(n_clicks):
    logging.debug("download_heatmap_image - n_clicks:" + str(n_clicks))
    if os.path.isfile(selected_heatmap_file):
        return [dcc.send_file(selected_heatmap_file)]
    else:
        return

# upload data table info
@app.callback(
    [Output("download-txt", "data")],
    [Input("btn-download-txt", "n_clicks")],
    prevent_initial_call=True
)
def download_data(n_clicks):
    logging.debug("download_data - n_clicks:" + str(n_clicks))
    global selected_dl_file
    return [dcc.send_file(selected_dl_file)]

@app.callback(
    [Output('container-button-timestamp', 'children'),
     Output('graph', 'figure'),
     Output('cytoscape-usage-phylogeny', 'elements'),
     Output('textarea-example-output', 'children'),
     # https://stackoverflow.com/questions/62375102/problem-dropping-same-file-twice-in-a-row
     Output('upload-data', 'contents'),
     Output('upload-data', 'filename'),
     Output("loading-output-1", "children")],
    [Input('button1', 'n_clicks'),
     Input('button2', 'n_clicks'),
     Input('upload-data', 'contents'),
     Input('run_clustalo_option', 'value'),
     Input('root_tree_option', 'value')
     ],
    prevent_initial_call=True
)
def update_data(n_clicks1, n_clicks2, upload_contents, run_clustalo_option, root_tree_option):
    logging.debug("update_data - n_clicks1:" + str(n_clicks1) + 
                  " n_clicks2:" + str(n_clicks2) + 
                  " contents2:" + str(upload_contents),
                  " run_clustalo_option:" + str(run_clustalo_option),
                  " root_tree_option:" + str(root_tree_option))
    global df_to_merge_loaded
    global selected_dl_file

    logging.debug("df_to_merge loaded: " + str(df_to_merge_loaded))
    
    df1 = pd.read_csv("hs2.pim.txt", skiprows=1, header=None, delim_whitespace=True)
    df2 = pd.read_csv("lyco_ochr.pim.txt", skiprows=1, header=None, delim_whitespace=True)
    button_id = ctx.triggered[0]['prop_id'].split(".")[0]
    logging.debug("button_id: " + button_id)
    
    df = pd.DataFrame()
    root_tree = len(root_tree_option) > 0

    if button_id == "button1":
        df = df1
        selected_dl_file = "./hs2.pim.txt"
        alignfile = "hs2.clu"
        with open(alignfile, "r") as aln:
            alignment = AlignIO.read(aln, "clustal")
            tree = generate_tree_from_alignment(alignment, root_tree)
            has_tree = False
            out_file_treetxt = data_path + "/hs2.tree.txt"
            with open(out_file_treetxt, "w") as tree_file:
                tree_file.write(str(tree))
                has_tree = True
        # zip up the outputs
        out_file_zip = data_path + "/hs2.zip"
        with zipfile.ZipFile(out_file_zip, 'w') as outzip:
            outzip.write(alignfile)
            outzip.write(selected_dl_file)
            if has_tree:
                outzip.write(out_file_treetxt, "hs2.tree.txt")
            logging.debug("zipped output files: " + out_file_zip)
            selected_dl_file = out_file_zip
        msg = "bZIP1 PIM ready for download"
    
    elif button_id == "button2":
        df = df2
        selected_dl_file = "./lyco_ochr.pim.txt"
        alignfile = "lyco_ochr.clu"
        with open(alignfile, "r") as aln:
            alignment = AlignIO.read(aln, "clustal")
            tree = generate_tree_from_alignment(alignment, root_tree)
            has_tree = False
            out_file_treetxt = data_path + "/lyco_ochr.tree.txt"
            with open(out_file_treetxt, "w") as tree_file:
                tree_file.write(str(tree))
                has_tree = True
        # zip up the outputs
        out_file_zip = data_path + "/lyco_ochr.zip"
        with zipfile.ZipFile(out_file_zip, 'w') as outzip:
            outzip.write(alignfile)
            outzip.write(selected_dl_file)
            if has_tree:
                outzip.write(out_file_treetxt, "lyco_ochr.tree.txt")
            logging.debug("zipped output files: " + out_file_zip)
            selected_dl_file = out_file_zip
        msg = "LycoCyclase PIM ready for download"
    
    elif button_id == "run_clustalo_option":
        if len(run_clustalo_option) == 0:
            msg = "Clustalo option disabled"
        else:
            msg = "Clustalo option enabled"
    # elif button_id == "btn-download-txt":
    #     msg = "Download clicked"

    # create fig
    # Will have the user text input the location of the first one! And then it will look for the fasta files relatively to this, as well as

    else:
        if upload_contents is not None and df_to_merge_loaded:
            try:
                content_type, content_string = upload_contents.split(",")
                decoded = base64.b64decode(content_string)
                decoded = decoded.decode('utf-8').replace("\\n", "\n")
                fasta_sio = io.StringIO(decoded)
                #prot_file = (pd.read_csv(fasta_sio))
                fasta_sio.seek(0)

                id_dict = {
                    "jgi": [],
                    "species": [],
                    "proteinId": [],
                    "modelId": [],
                }
                seq_dict = {
                    "modelId": [],
                    "proteinSeq": [],
                }
                for record in SeqIO.parse(fasta_sio, "fasta"):
                    #logging.debug(record.id.split("|"))
                    #logging.debug(record.seq)
                    record_id_split = record.id.split("|", 3)
                    id_dict["jgi"].append(record_id_split[0])
                    id_dict["species"].append(record_id_split[1])
                    id_dict["proteinId"].append(record_id_split[2])
                    id_dict["modelId"].append(record_id_split[3])
                    seq_dict["modelId"].append(record_id_split[3])
                    seq_dict["proteinSeq"].append(record.seq)
                
                logging.debug("loaded IDs and sequences from fasta file")

                # prot_file_2 = prot_file.iloc[:, 0].str.split("\n", expand=True)
                # prot_file_1 = prot_file.iloc[:, 0].str.split("|", 3, expand=True)
                #print("pf1: ", prot_file_1.info())
                #print("pf2: ",prot_file_2.info())
                # prot_ce = prot_file_2.iloc[:, 0].str.split("|", 3, expand=True)
                #print("pce:", prot_ce.info())
                # prot_ce.columns = ["jgi", "species", "proteinId", "modelId"]

                # P_sequence = prot_file_1.iloc[:, 3].str.split("\n", n=1, expand=True)
                # P_sequence.columns = ["modelId", "proteinSeq"]

                prot_ce = pd.DataFrame(id_dict)
                prot_ce["proteinId"] = pd.to_numeric(prot_ce["proteinId"])
                
                P_sequence = pd.DataFrame(seq_dict)
                logging.debug("created prot_ce and P_sequence DFs")

                desc_seq_tet = P_sequence.merge(prot_ce, how="inner", on="modelId")
                desc_seq_tet['proteinId'] = desc_seq_tet['proteinId'].astype('str')
                logging.debug("merged desc_seq_tet")
                
                merged_to_fasta = desc_seq_tet.merge(df_to_merge, how="left", on=["species", "proteinId"])
                logging.debug("merged to_fasta")
                merged_all = merged_to_fasta.merge(df_to_merge, how="left", on=["group_id"])
                logging.debug("merged all")

                merged_all['group_id'] = merged_all['group_id'].astype('str')
                merged_all['proteinId_y'] = merged_all['proteinId_y'].astype('str')

                merged_all["label"] = merged_all[['jgi', 'species_y', 'proteinId_y', 'group_id']].apply(
                    lambda x: '|'.join(x[x.notnull()]), axis=1)
                logging.debug("created labels")

                to_fasta = merged_all[['label', 'proteinSeq']].copy()

                to_fasta['label'] = '>' + to_fasta['label'].astype(str)
                #logging.debug(df_to_merge.head())
                #logging.debug(to_fasta.head())
                to_fasta['proteinSeq'] = to_fasta['proteinSeq'].astype(str)

                #if there are two sequences with the same name, only the longer of the two remains.
                to_fasta.sort_values(by='proteinSeq', key=lambda x: x.str.len(),inplace = True)
                to_fasta.drop_duplicates(subset='label', keep = 'last',inplace = True)

                logging.debug("removed duplicates")

                # to_fasta['proteinSeq'] = to_fasta['proteinSeq'].str.replace('*', '', regex=True)
                #logging.debug(to_fasta.head())
                logging.debug("created to_fasta")


                # get a unique filename
                filename = "ortho_" + str(int(time.time()))

                # get the full path of each file
                in_file_fasta = data_path + "/" + filename + ".fasta"
                out_file_clu = data_path + "/" + filename + ".clu"
                out_file_pimtxt = data_path + "/" + filename + ".pim.txt"
                out_file_treetxt = data_path + "/" + filename + ".tree.txt"
                out_file_zip = data_path + "/" + filename + ".zip"
                
            
    
                with open(in_file_fasta, 'w') as f:
                    for index, row in to_fasta.iterrows():
                        #logging.debug(str(row))
                        #print(row["label"], f'\n', row["proteinSeq"], file=in_file_fasta)
                        try:
                            f.write(row["label"] + "\n" + row["proteinSeq"] + "\n")
                        except Exception as e:
                            #logging.error("error writing to_fasta: " + str(e))
                            pass
                logging.debug("wrote to_fasta to " + in_file_fasta)

                selected_dl_file = in_file_fasta

                # run clustalo ONLY if the option is selected
                # otherwise, just return the fasta file
                if len(run_clustalo_option) == 0:
                    msg = 'Merged fasta file ready for download'
                    long_msg = 'Merged fasta file created, clustalo not run. Merged fasta can be downloaded using the "Download Results" button.'
                    # msg_fig = ff.create_table([[msg],[long_msg]], height_constant=20)
                    # msg_fig.layout.annotations[0].font.size = 20
                    
                    # Create an empty matrix
                    z_data = [[0]]
                    # Create empty annotated heatmap figure
                    fig = ff.create_annotated_heatmap(z=z_data, colorscale='deep')
                    fig.update_layout(title=msg)
                    return msg, fig, [], long_msg, None, None, no_update

                else:
                    # figure out the number of theads for clustalo
                    threads = multiprocessing.cpu_count()
                    if threads < 1:
                        threads = 1
                    elif threads > 8:
                        threads = 8
                    logging.debug("clustalo threads: " + str(threads))

                    # run clustalo the first time to get the clu file
                    logging.debug("running clustalo with input: " + in_file_fasta)
                    subprocess.run(["clustalo", "-i", in_file_fasta, "-o", out_file_clu, "--outfmt", "clu", "--threads", str(threads), "-v", "-v", "-v"], 
                                text=True, check=True)
                    logging.debug("ran clustalo with output: " + out_file_clu)

                    # run clustalo the second time to get the distance matrix file
                    # from docs: if no alignment is desired but only distance calculation and tree construction, 
                    #            then --max-hmm-iterations=-1 will terminate the calculation before the alignment stage
                    logging.debug("running clustalo with input: " + out_file_clu)
                    subprocess.run(["clustalo", "-i", out_file_clu, "--percent-id", "--distmat-out=" + out_file_pimtxt, "--full", "--max-hmm-iterations=-1"],
                                text=True, check=True)
                    logging.debug("ran clustalo with output: " + out_file_pimtxt)

                    # read the clu/alignment file
                    alignfile = out_file_clu
                    has_tree = False
                    with open(alignfile, "r") as aln:
                        alignment = AlignIO.read(aln, "clustal")
                        tree = generate_tree_from_alignment(alignment, root_tree)
                        with open(out_file_treetxt, "w") as tree_file:
                            tree_file.write(str(tree))
                            has_tree = True

                    logging.debug("loaded alignment from: " + alignfile)

                    # read the pim.txt (distance matrix) file
                    df3 = pd.read_csv(out_file_pimtxt, skiprows=1, header=None, delim_whitespace=True)
                    df = df3
                    selected_dl_file = out_file_pimtxt
                    logging.debug("loaded distance matrix df: " + str(df))

                    # zip up the outputs
                    with zipfile.ZipFile(out_file_zip, 'w') as outzip:
                        outzip.write(in_file_fasta)
                        outzip.write(out_file_clu)
                        outzip.write(out_file_pimtxt)
                        if has_tree:
                            outzip.write(out_file_treetxt)
                        logging.debug("zipped output files: " + out_file_zip)
                        selected_dl_file = out_file_zip

                    msg = 'Full results ready for download'
                    long_msg = 'Merged fasta file created, clustalo run and alignment created. Results zip can be downloaded using the button.'

            except Exception as e:
                msg = 'An error occurred'
                long_msg = 'An error occurred. Please check your input file, refresh the page, and try again.'
                #logging.error(long_msg)
                logging.error(str(e))
                # print stack trace 
                traceback.print_exc()
                
                #return html.Div(id=msg)
                error_fig = ff.create_table([[msg],[long_msg]], height_constant=20)
                # Make text size larger
                error_fig.layout.annotations[0].font.size = 20
                return msg, error_fig, [], long_msg, None, None, no_update


    if not df.empty:
        df = df.T
        df = df.rename(columns=df.iloc[0]).drop(df.index[0])
        z = df.to_numpy().astype(np.float64)
        z_text = np.zeros(z.shape)  # np.round(z, 2)
        x_text = df.columns.tolist()
        y_text = df.columns.tolist()
        annotations = [['' for _ in range(len(z[0]))] for _ in range(len(z))]

        # take text off of figure
        fig = ff.create_annotated_heatmap(
            z, 
            x=x_text, 
            y=y_text, 
            annotation_text=annotations, 
            colorscale='deep'
        )
        fig.update_xaxes(side="bottom")
        fig['data'][0]['showscale'] = True

        # fig_time = str(int(time.time()))
        if SAVE_FIGURES:
            # if not os.path.exists(data_path + "/figures"):
            #     os.mkdir(data_path + "/figures")
            # fig_filename = data_path + "/figures/fig_heat_" + fig_time + ".png"
            fig.write_image(
                # file=fig_filename, 
                file=selected_heatmap_file,
                format="png",
                width=FIGURE_WIDTH, 
                height=FIGURE_HEIGHT)

        # Define elements, stylesheet and layout
        # download this from http://www.phyloxml.org/examples/apaf.xml
        temp_node, temp_edge = generate_elements(tree)
        elements = temp_node + temp_edge

        newick_string = str(tree)
        return msg, fig, elements, newick_string, None, None, no_update
    else:
        return no_update, no_update, no_update, no_update, None, None, no_update
    #    break


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AlgaeOrtho')
    
    # noserver mode will perform the pickle file logic only, it won't start the server
    # this is used mainly for the docker image. the pickle file will be created during
    # the docker image build process to make the docker run process faster
    parser.add_argument('--noserver', action='store_true')
    parser.set_defaults(noserver=False)
    args = parser.parse_args()
    logging.info("no server mode: " + str(args.noserver))

    DF_PICKLE_FILE = "df_to_merge.pickle"
    if os.path.exists(DF_PICKLE_FILE):
        # if the pickle file exists, load the dataframe from it
        # this is much faster than loading the dataframe from the original files
        df_to_merge = pickle.load( open( DF_PICKLE_FILE, "rb" ) )
        logging.debug("df_to_merge loaded from pickle")
        df_to_merge_loaded = True

    if not df_to_merge_loaded:
        # if the pickle file doesn't exist, load the dataframe and save them to a pickle file
        # this will speed up the app start next time
        df_to_merge = load_df_to_merge()
        df_to_merge_loaded = True
        pickle.dump( df_to_merge, open( DF_PICKLE_FILE, "wb" ) )

    # start the server if not in noserver mode
    if not args.noserver:
        envPort = int(os.getenv('SERVER_PORT', '8050'))
        envDebug = os.getenv('DASH_DEBUG_MODE', 'True').lower() == 'true'
        app.run_server(debug=envDebug, host='0.0.0.0', port=envPort)
