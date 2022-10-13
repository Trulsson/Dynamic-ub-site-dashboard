import pandas as pd
import dash
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
import math
from dash import html, dcc, dash_table, no_update
import dash_bootstrap_components as dbc
import numpy as np


# Helper function to transform -log10 p-values
def power_to_transform(x):
    return 10 ** -x


# and back to -log10
def neg_log10(x):
    return np.log10(1 / x)


# Import Ubisite DDA Data
df = pd.read_csv('UbiSiteDDA.csv')

df['gene+pos'] = df['Gene names'] + "-K" + df['Positions within proteins'].astype(str)
df_tak = df[df['TAK Difference (log2)'].notnull()]
df_mg = df[df['MG Difference (log2)'].notnull()]
df_pr = df[df['PR Difference (log2)'].notnull()]

# Columns in data
p_cols = ['TAK -Log p-value', 'MG -Log p-value', 'PR -Log p-value']
fc_cols = ['TAK Difference (log2)', 'MG Difference (log2)', 'PR Difference (log2)']
df_all = df.drop(['Proteins', 'Protein names',
                  'TAK Significant (FDR=0.05 S0=0.1)', 'TAK q-value',
                  'TAK all sites average',
                  'MG Significant (FDR=0.05 S0=0.1)', 'MG q-value',
                  'MG all sites average',
                  'PR Significant (FDR=0.05 S0=0.1)', 'PR q-value',
                  'PR all sites average', 'Ub Enzyme', 'Enzyme',
                  'Protein in His10Ub', 'PEP', 'Delta score',
                  'gene+pos'], axis=1).copy()

df_all.rename(columns={'TAK -Log p-value': 'TAK243 p', 'TAK Difference (log2)': 'TAK243 fc',
                       'PR -Log p-value': 'PR619 p', 'PR Difference (log2)': 'PR619 fc',
                       'MG -Log p-value': 'MG132 p', 'MG Difference (log2)': 'MG132 fc',
                       'Score for localization': 'Pos. Score', 'Sequence window': 'Sequence',
                       'Positions within proteins': 'Lysine', 'Gene names': 'Gene'
                       }, inplace=True)

# Transform p-values
df_tak_p = df_tak.apply(lambda x: power_to_transform(x) if x.name in p_cols else x)
df_tak_p.reset_index(inplace=True)
df_mg_p = df_mg.apply(lambda x: power_to_transform(x) if x.name in p_cols else x)
df_mg_p.reset_index(inplace=True)
df_pr_p = df_pr.apply(lambda x: power_to_transform(x) if x.name in p_cols else x)
df_pr_p.reset_index(inplace=True)

# Initialize dashboard
app = dash.Dash(
    external_stylesheets=[dbc.themes.BOOTSTRAP]
)

# Layout of app
app.layout = dbc.Container([
    # Title
    dbc.Row([
        dbc.Col([
            html.H2('UbiSite DDA Data',
                    className='text-center text-primary, mb-4'),
        ])
    ]),
    # Input field
    dbc.Row([
        dbc.Col([
            html.Div([
                dcc.Input(id="gene_input", type="text", placeholder="Enter Gene Name",
                          debounce=True, className='col-3'),
                dbc.Tooltip(
                    "Search for identified Ubiquitin sites in your Gene of interest in this dataset",
                    target="gene_input",
                    placement='top',
                )
            ])
        ]),
    ]),
    # Alert for Gene not found
    dbc.Row([
        dbc.Col([
            html.Div([
                dbc.Alert(
                    html.H6("Gene Name not found", className='text-center'),
                    className='col-3 position-relative"',
                    is_open=False,
                    fade=False,
                    color="danger",
                    dismissable=True,
                    id='err',
                )
            ])
        ]),
    ], style={"height": "40px"}),
    # Sliders and Volcano
    dbc.Row([
        # TAK
        dbc.Col(
            html.Div([
                html.H3(children='Conjugation Inhibition',
                        className='text-center text-primary, mb-4',
                        style={'margin-right': '2px', 'margin-left': '2px'},
                        id='tak-title'
                        ),
                dbc.Tooltip(
                    "Inhibition by TAK-243 for 3h",
                    target="tak-title",
                    placement='top',
                ),
                # P-value slider TAK
                html.H5('P-value'),
                dcc.Slider(
                    tooltip={'always_visible': True},
                    id='tak-volcanoplot-input_p',
                    min=0,
                    max=max(df_tak['TAK -Log p-value']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(0, math.ceil(max(df_tak['TAK -Log p-value'])))},
                    value=1.3
                ),
                # Fold Change slider TAK
                html.H5('Fold Change'),
                dcc.RangeSlider(
                    tooltip={'always_visible': True},
                    id='tak-volcanoplot-input_fc',
                    min=min(df_tak_p['TAK Difference (log2)']),
                    max=max(df_tak_p['TAK Difference (log2)']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(math.floor(min(df_tak_p['TAK Difference (log2)'])),
                                                               math.ceil(max(df_tak_p['TAK Difference (log2)'])))},
                    value=[-2, 2]
                ),
                # Volcano TAK
                dcc.Graph(id='ubisite-tak-volcanoplot'),
            ]),
            md=12, lg=4
        ),
        # MG
        dbc.Col(
            html.Div([
                html.H3(children='Proteasome Inhibition',
                        className='text-center text-primary, mb-4',
                        id='mg-title'),
                dbc.Tooltip(
                    "Inhibition by MG-132 for 3h",
                    target="mg-title",
                    placement='top',
                ),
                # P-value slider MG
                html.H5('P-value'),
                dcc.Slider(
                    tooltip={'always_visible': True},
                    id='mg-volcanoplot-input_p',
                    min=0,
                    max=max(df_mg['MG -Log p-value']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(0, math.ceil(max(df_mg['MG -Log p-value'])))},
                    value=1.3
                ),
                # Fold Change slider MG
                html.H5('Fold Change'),
                dcc.RangeSlider(
                    tooltip={'always_visible': True},
                    id='mg-volcanoplot-input_fc',
                    min=min(df_mg_p['MG Difference (log2)']),
                    max=max(df_mg_p['MG Difference (log2)']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(math.floor(min(df_mg_p['MG Difference (log2)'])),
                                                               math.ceil(max(df_mg_p['MG Difference (log2)'])))},
                    value=[-2, 2]
                ),
                # Volcano MG
                dcc.Graph(id='ubisite-mg-volcanoplot')
            ]),
            md=12, lg=4
        ),
        # PR
        dbc.Col(
            html.Div([
                html.H3(children='DUB Inhibition',
                        className='text-center text-primary, mb-4',
                        id='pr-title'),
                dbc.Tooltip(
                    "Inhibition by PR619 for 3h",
                    target="pr-title",
                    placement='top',
                ),
                # P-value slider PR
                html.H5('P-value'),
                dcc.Slider(
                    tooltip={'always_visible': True},
                    id='pr-volcanoplot-input_p',
                    min=0,
                    max=max(df_pr['PR -Log p-value']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(0, math.ceil(max(df_pr['PR -Log p-value'])))},
                    value=1.3
                ),
                # Fold Change slider PR
                html.H5('Fold Change'),
                dcc.RangeSlider(
                    tooltip={'always_visible': True},
                    id='pr-volcanoplot-input_fc',
                    min=min(df_pr_p['PR Difference (log2)']),
                    max=max(df_pr_p['PR Difference (log2)']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(math.floor(min(df_pr_p['PR Difference (log2)'])),
                                                               math.ceil(max(df_pr_p['PR Difference (log2)'])))},
                    value=[-2, 2]
                ),
                # Volcano PR
                dcc.Graph(id='ubisite-pr-volcanoplot')
            ]),
            md=12, lg=4
        )
    ]),
    # Data Table
    dbc.Row([
        dbc.Col([
            html.Div([
                html.Br(),
                dash_table.DataTable(
                    id='poi-table',
                    columns=[{"name": i, "id": i} for i in df_all.columns],
                    style_table={'height': '300px', 'overflowY': 'auto'},
                    page_action='none',
                    sort_action='custom',
                    fixed_rows={'headers': True},
                    sort_mode='single',
                    sort_by=[],
                    cell_selectable=True,
                    export_format='csv',
                    export_columns='all',
                    row_deletable=True,
                    tooltip_header={'Lysine': 'Ub modified position within protein',
                                    'Gene': 'Peptides that match to several proteins are separated by ";"',
                                    'TAK243 p': 'Ub conjugation inhibition p-value (-log10)',
                                    'TAK243 fc': 'Ub conjugation inhibition fold-change (log2)',
                                    'MG132 p': 'Proteasome inhibition p-value (-log10)',
                                    'MG132 fc': 'Proteasome inhibition fold-change (log2)',
                                    'PR619 p': 'Deubiquitinating enzyme inhibition p-value (-log10)',
                                    'PR619 fc': "Deubiquitinating enzyme inhibition fold-change (log2)",
                                    'Score': 'Peptide Score: Higher value means higher confidence',
                                    'Pos. Score': 'Position Score: Higher value means higher confidence in '
                                                  'localization of modification within peptides',
                                    'Sequence': 'Amino acid sequence surrounding modification site'},
                    # Overflow into ellipsis
                    style_cell={
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                        'minWidth': 95, 'maxWidth': 95, 'width': 95
                    },
                    tooltip_delay=0,
                    tooltip_duration=None
                )
            ])
        ])
    ])
])


@app.callback(
    Output('err', 'is_open'),
    Output('gene_input', 'value'),
    Input('gene_input', 'value'),
    State("err", "is_open")
)
def check_input(input_gene, is_open):
    if input_gene is None:
        return no_update, no_update
    print("Search for Gene: ", input_gene)
    poi = df_all[df_all['Gene'].str.contains(input_gene.upper()) == True]
    if len(poi.index) == 0:
        is_open = True
        return is_open, no_update
    is_open = False
    return is_open, input_gene


@app.callback(
    Output('ubisite-tak-volcanoplot', 'figure'),
    Input('tak-volcanoplot-input_p', 'value'),
    Input('tak-volcanoplot-input_fc', 'value'),
    Input('gene_input', 'value')
)
def update_volcanoplot_tak(p_val, fc, input_gene):
    fig = dashbio.VolcanoPlot(
        dataframe=df_tak_p,
        effect_size='TAK Difference (log2)',
        p='TAK -Log p-value',
        gene='gene+pos',
        point_size=8,
        xlabel='Log2 Difference',
        effect_size_line=fc,
        genomewideline_value=p_val,
        highlight_color='#119DFF',
        col='#2A3F5F',
        snp='gene+pos',
    )
    # Hover data information
    single_gene_names = []
    for i, string in enumerate(df_tak_p['gene+pos']):
        single_gene_names.append(string.split(';', 1)[0])
    customdata = np.stack((single_gene_names, df_tak_p['Score']), axis=-1)
    hovertemplate = '<b>Gene: %{customdata[0]}</b> <br>Fold Change: %{x} <br>p-value: %{y} <br>Score:%{customdata[1]}'
    fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)

    # Figure layout
    fig.update_layout(showlegend=False, title_text=None,
                      margin={'l': 10, 'r': 10, 'b': 30, 't': 30, 'pad': 4},
                      )

    # Update graph with Gene of interest
    if input_gene:
        p_vals = df_tak['TAK -Log p-value'].loc[df_tak['Gene names'].str.contains(input_gene.upper()) == True]
        folds = df_tak['TAK Difference (log2)'].loc[df_tak['Gene names'].str.contains(input_gene.upper()) == True]
        name = df_tak['gene+pos'].loc[df_tak['Gene names'].str.contains(input_gene.upper()) == True]
        texts = []
        for fold, p, name in zip(folds, p_vals, name):
            texts.append(fig.add_annotation(x=fold, y=p, text=name, bgcolor=(
                '#69b8f0' if (p > p_val) and any([(fold < fc[0]), (fold > fc[1])]) else 'white')))

    return fig


@app.callback(
    Output('ubisite-mg-volcanoplot', 'figure'),
    Input('mg-volcanoplot-input_p', 'value'),
    Input('mg-volcanoplot-input_fc', 'value'),
    Input('gene_input', 'value')
)
def update_volcanoplot_mg(p_val, fc, input_gene):
    fig = dashbio.VolcanoPlot(
        dataframe=df_mg_p,
        effect_size='MG Difference (log2)',
        p='MG -Log p-value',
        gene='gene+pos',
        point_size=8,
        xlabel='Log2 Difference',
        effect_size_line=fc,
        genomewideline_value=p_val,
        highlight_color='#119DFF',
        col='#2A3F5F',
        snp='Sequence window',
        annotation='Score'
    )
    # Hover data information
    single_gene_names = []
    for i, string in enumerate(df_mg_p['gene+pos']):
        single_gene_names.append(string.split(';', 1)[0])
    customdata = np.stack((single_gene_names, df_mg_p['Score']), axis=-1)
    hovertemplate = '<b>Gene: %{customdata[0]}</b> <br>Fold Change: %{x} <br>p-value: %{y} <br>Score:%{customdata[1]}'
    fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)

    # Figure layout
    fig.update_layout(showlegend=False, title_text=None, margin={'l': 10, 'r': 10, 'b': 30, 't': 30, 'pad': 4})

    # Update graph with Gene of interest
    if input_gene:
        p_vals = df_mg['MG -Log p-value'].loc[df_mg['Gene names'].str.contains(input_gene.upper()) == True]
        folds = df_mg['MG Difference (log2)'].loc[df_mg['Gene names'].str.contains(input_gene.upper()) == True]
        name = df_mg['gene+pos'].loc[df_mg['Gene names'].str.contains(input_gene.upper()) == True]
        texts = []
        for fold, p, name in zip(folds, p_vals, name):
            texts.append(fig.add_annotation(x=fold, y=p, text=name, bgcolor='white'))

    return fig


@app.callback(
    Output('ubisite-pr-volcanoplot', 'figure'),
    Input('pr-volcanoplot-input_p', 'value'),
    Input('pr-volcanoplot-input_fc', 'value'),
    Input('gene_input', 'value')
)
def update_volcanoplot_pr(p_val, fc, input_gene):
    fig = dashbio.VolcanoPlot(
        dataframe=df_pr_p,
        effect_size='PR Difference (log2)',
        p='PR -Log p-value',
        gene='gene+pos',
        point_size=8,
        xlabel='Log2 Difference',
        effect_size_line=fc,
        genomewideline_value=p_val,
        highlight_color='#119DFF',
        col='#2A3F5F',
        snp='Sequence window',
        annotation='Score'
    )
    # Hover data information
    single_gene_names = []
    for i, string in enumerate(df_pr_p['gene+pos']):
        single_gene_names.append(string.split(';', 1)[0])
    customdata = np.stack((single_gene_names, df_pr_p['Score']), axis=-1)
    hovertemplate = '<b>Gene: %{customdata[0]}</b> <br>Fold Change: %{x} <br>p-value: %{y} <br>Score:%{customdata[1]}'
    fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)

    # Figure layout
    fig.update_layout(showlegend=False, title_text=None, margin={'l': 10, 'r': 10, 'b': 30, 't': 30, 'pad': 4})

    # Update graph with Gene of interest
    if input_gene:
        p_vals = df_pr['PR -Log p-value'].loc[df_pr['Gene names'].str.contains(input_gene.upper()) == True]
        folds = df_pr['PR Difference (log2)'].loc[df_pr['Gene names'].str.contains(input_gene.upper()) == True]
        name = df_pr['gene+pos'].loc[df_pr['Gene names'].str.contains(input_gene.upper()) == True]
        texts = []
        for fold, p, name in zip(folds, p_vals, name):
            texts.append(fig.add_annotation(x=fold, y=p, text=name, bgcolor='white'))

    return fig


@app.callback(
    Output('poi-table', "data"),
    Input('poi-table', "sort_by"),
    Input('gene_input', 'value'))
def update_table(sort_by, input_gene):
    if len(sort_by) and input_gene:
        poi = df_all[df_all['Gene'].str.contains(input_gene.upper()) == True].sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )
    elif input_gene:
        # No sort is applied
        poi = df_all[df_all['Gene'].str.contains(input_gene.upper()) == True]

    else:
        return

    return poi.to_dict('records')


if __name__ == '__main__':
    app.run_server(debug=True)
