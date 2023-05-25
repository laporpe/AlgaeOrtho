# Dash App Prep: imports, initializing, dataset!

import base64
import io
import os
from dash import Dash, html, dcc, Input, Output, ctx, dash_table
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
import pandas as pd
import plotly.figure_factory as ff
import numpy as np
import math
from Bio import AlignIO, Phylo
import dash_cytoscape as cyto
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

import subprocess

# Original Input File will be the results from SonicParanoid
#

# Start incorporating backend:
orthodf = pd.read_csv("ortholog_groups.tsv", sep = "\t")

horizontal = orthodf
shape_num = horizontal.shape[1]

horizontal.columns[4:shape_num].values[0:] = ([entry[:-9] for entry in horizontal.columns[4:shape_num]])

pd.concat([horizontal.iloc[:,0:3],horizontal.iloc[:,4]], axis = 1)


d = {}
for i in range(4,shape_num):
    [entry[:-9] for entry in horizontal.columns[4:shape_num]]
    d[horizontal.keys()[i]] = pd.DataFrame(pd.concat([horizontal.iloc[:,0:3],horizontal.iloc[:,i]], axis = 1))

d2 = dict.fromkeys(d.keys())
newdf = pd.DataFrame()
count = 0
for var in d.keys():
    frame = d[var]
    tester_df = frame.replace("*",np.NaN)
    df1 = (tester_df.assign(category = tester_df.iloc[:,3].str.split(','))
         .explode('category')
         .reset_index(drop=True))
    df1 = pd.concat([df1,df1['category'].str.split('|', expand = True)], axis = 1)
    df2 = df1.iloc[:,0:8]
    df2.columns = ["group_id","group_size","sp_in_grp", "full_entry", "category", "jgi", "species", "proteinId"]
    newdf = pd.concat([newdf, df2], axis = 0)
    d2[var] = df2

newdf

concise_df = newdf[newdf['proteinId'].notna()]

df_to_merge = concise_df[["group_id","species","proteinId"]]

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
            'background-color': '#222222'
        }
    }
]

with open("lyco_ochr.clu", "r") as aln:
    alignment = AlignIO.read(aln, "clustal")

calculator = DistanceCalculator('identity', )
distance_matrix = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor(calculator)
tree = constructor.build_tree(alignment)

newick_string = str(tree)


def generate_elements(tree, xlen=30, ylen=30, grabbable=False):
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
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = ((positions[clade.clades[0]] +
                                 positions[clade.clades[-1]]) // 2)

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

app = Dash(external_stylesheets=[dbc.themes.COSMO])

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "24rem",
    "padding": "2rem 1rem",
    "background-color": "#7393B3",
}

sidebar = html.Div(
    [
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
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=False
        ),
        html.Hr(),

    ],
    style=SIDEBAR_STYLE,
)

app.layout = html.Div(children=[
    dbc.Row([
        dbc.Col(sidebar),
        dbc.Col(html.H1(children='AlgaeOrtho'))
    ]),
    dbc.Row([
        dbc.Col(sidebar),
        dbc.Col(dcc.Graph(id='graph', style={'width': '100%', 'height': '90vh'}))
    ]),
    # dbc.Row([
    #     dbc.Col(sidebar),
    #     dbc.Col([html.Button("Download Percent Identity Matrix", id="btn-download-txt"),
    #              dcc.Download(id="table")])
    # ]),
    dbc.Row([
        dbc.Col(sidebar),
        dbc.Col(cyto.Cytoscape(
            id='cytoscape-usage-phylogeny',
            elements=elements,
            stylesheet=stylesheet,
            layout=layout,
            style={
                'height': '95vh',
                'width': '100%'
            }
        ))
    ]),

    dbc.Row([
        dbc.Col(),
        dbc.Col(html.Div(id='textarea-example-output', style={'whiteSpace': 'break-spaces'}))
    ]),
    ])


# upload data table info


@app.callback(
    # [# Output('container-button-timestamp', 'children'),
    # Output('graph', 'figure')],
    [Output('container-button-timestamp', 'children'),
     # Output('output-data-upload', 'children'),
     Output('graph', 'figure'),
     #Output('table', 'data'),
     # Output("download-dataframe-csv", "data"),
     Output('cytoscape-usage-phylogeny', 'elements'),
     Output('textarea-example-output', 'children')],
    [Input('button1', 'n_clicks'),
     Input('button2', 'n_clicks'),
     #Input("btn-download-txt", "n_clicks3"),
     Input('upload-data', 'contents')
     # Input('textarea-example', 'value')
     ],
    prevent_initial_call=True
)

def update_data(n_clicks1, n_clicks2, contents2):
    df1 = pd.read_csv("hs2.pim.txt", skiprows=1, header=None, delim_whitespace=True)
    df2 = pd.read_csv("lyco_ochr.pim.txt", skiprows=1, header=None, delim_whitespace=True)
    button_id = ctx.triggered[0]['prop_id'].split(".")[0]
    msg = ""
    df = df1
    datar = df1
    alignfile = "hs2.clu"
    with open(alignfile, "r") as aln:
        alignment = AlignIO.read(aln, "clustal")

    if button_id == "button1":
        df = df1
        alignfile = "hs2.clu"
        with open(alignfile, "r") as aln:
            alignment = AlignIO.read(aln, "clustal")
        msg = "bZIP1 clicked"
    elif button_id == "button2":
        df = df2
        alignfile = "lyco_ochr.clu"
        with open(alignfile, "r") as aln:
            alignment = AlignIO.read(aln, "clustal")
        msg = "LycoCyclase clicked"
    # create fig
    # Will have the user text input the location of the first one! And then it will look for the fasta files relatively to this, as well as

    elif contents2 is not None:
        content_type, content_string = contents2.split(",")
        decoded = base64.b64decode(content_string)
        try:
            for_phylogeny = (pd.read_csv(io.StringIO(decoded.decode('utf-8'))))

            prot_file = for_phylogeny

            prot_file_2 = prot_file.iloc[:, 0].str.split("\n", expand=True)
            prot_file_1 = prot_file.iloc[:, 0].str.split("|", 3, expand=True)

            prot_ce = prot_file_2.iloc[:, 0].str.split("|", 3, expand=True)
            # print(prot_ce.head())
            prot_ce.columns = ["jgi", "species", "proteinId", "modelId"]
            prot_ce["proteinId"] = pd.to_numeric(prot_ce["proteinId"])
            P_sequence = prot_file_1.iloc[:, 3].str.split("\n", n=1, expand=True)
            P_sequence.columns = ["modelId", "proteinSeq"]

            desc_seq_tet = P_sequence.merge(prot_ce, how="inner", on="modelId")
            desc_seq_tet['proteinId'] = desc_seq_tet['proteinId'].astype('str')

            merged_to_fasta = desc_seq_tet.merge(df_to_merge, how="left", on=["species", "proteinId"])

            merged_all = merged_to_fasta.merge(df_to_merge, how="left", on=["group_id"])

            merged_all['group_id'] = merged_all['group_id'].astype('str')
            merged_all['proteinId_y'] = merged_all['proteinId_y'].astype('str')

            merged_all["label"] = merged_all[['jgi', 'species_y', 'proteinId_y', 'group_id']].apply(
                lambda x: '|'.join(x[x.notnull()]), axis=1)

            to_fasta = merged_all[['label', 'proteinSeq']].copy()

            to_fasta['label'] = '>' + to_fasta['label'].astype(str)
            to_fasta['proteinSeq'] = to_fasta['proteinSeq'].str.replace('*', '', regex=True)

            f = io.StringIO()
            for index, row in to_fasta.iterrows():
                print(row["label"], f'\n', row["proteinSeq"], file = "./ortho.fasta")

            in_file = "./ortho.fasta"
            out_file2 = "./ortho.pim.txt"
            out_file4 = "./ortho.clu"

            subprocess.run(["clustalo","-i", "$in_file","-o","out_file4","--outfmt","clu"], text=True)
            subprocess.run(["clustalo", "-i", "$out_file4", "--percent-id","--distmat-out=$out_file2",  "--full"], text=True)

            alignfile = "ortho.clu"

            with open(alignfile, "r") as aln:
                alignment = AlignIO.read(aln, "clustal")

            df3 = pd.read_csv("ortho.pim.txt", skiprows=1, header=None, delim_whitespace=True)
            df = df3


        except Exception as e:
            print(e)
            # handle the





    def generate_elements(tree, xlen=30, ylen=30, grabbable=False):
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
                    if subclade not in positions:
                        calc_row(subclade)
                positions[clade] = ((positions[clade.clades[0]] +
                                     positions[clade.clades[-1]]) // 2)

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

    df = df.T
    df = df.rename(columns=df.iloc[0]).drop(df.index[0])
    z = df.to_numpy().astype(np.float64)
    z_text = np.zeros(z.shape)  # np.round(z, 2)
    x_text = df.columns.tolist()
    y_text = df.columns.tolist()
    annotations = [['' for _ in range(len(z[0]))] for _ in range(len(z))]

    # take text off of figure
    fig = ff.create_annotated_heatmap(z, x=x_text, y=y_text, annotation_text=annotations, colorscale='deep')
    fig.update_xaxes(side="bottom")
    fig['data'][0]['showscale'] = True

    # Define elements, stylesheet and layout
    # download this from http://www.phyloxml.org/examples/apaf.xml
    calculator = DistanceCalculator('identity', )
    distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.build_tree(alignment)

    newick_string = str(tree)
    temp_node, temp_edge = generate_elements(tree)
    elements = temp_node + temp_edge

    return msg, fig, elements, newick_string


if __name__ == '__main__':
    envPort = int(os.getenv('SERVER_PORT', '8050'))
    envDebug = os.getenv('DASH_DEBUG_MODE', 'True').lower() == 'true'
    app.run_server(debug=envDebug, host='0.0.0.0', port=envPort)
