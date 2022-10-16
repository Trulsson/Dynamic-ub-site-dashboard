import pandas as pd
import dash
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
import math
from dash import html, dcc, dash_table, no_update
import dash_bootstrap_components as dbc


# Helper function to transform -log10 p-values
def power_to_transform(x):
    return 10 ** -x


# Import Ubisite DDA Data
df = pd.read_csv('UbiSiteDIA.csv')
df['gene+pos'] = df['Genes'] + "-K" + df['Positions within proteins'].astype(str)

# Columns in data
df_selected = df.drop(['Proteins', 'Protein names',
                       'TAK Significant (FDR=0.05 S0=0.1)', 'TAK q-value',
                       'TAK all sites average',
                       'MG Significant (FDR=0.05 S0=0.1)', 'MG q-value',
                       'MG all sites average',
                       'PR Significant (FDR=0.05 S0=0.1)', 'PR q-value',
                       'PR all sites average', 'Ub Enzyme', 'Enzyme',
                       'Protein in His10Ub', 'PEP', 'Delta score', ], axis=1).copy()

df_selected.rename(columns={'TAK -Log p-value': 'TAK243 p', 'TAK Difference (log2)': 'TAK243 fc',
                            'PR -Log p-value': 'PR619 p', 'PR Difference (log2)': 'PR619 fc',
                            'MG -Log p-value': 'MG132 p', 'MG Difference (log2)': 'MG132 fc',
                            'Score for localization': 'Pos. Score', 'Sequence window': 'Sequence',
                            'Positions within proteins': 'Lysine', 'Gene names': 'Gene'
                            }, inplace=True)

p_cols = ['TAK243 p', 'MG132 p', 'PR619 p']
fc_cols = ['TAK243 fc', 'MG132 fc', 'PR619 fc']

# Transform p-values
for col in p_cols:
    df_selected[col + '_x'] = df_selected[col].apply(lambda x: power_to_transform(x))
df_selected.reset_index(inplace=True)

df_tak, df_mg, df_pr = df_selected[df_selected['TAK243 fc'].notnull()], df_selected[df_selected['MG132 fc'].notnull()], \
                       df_selected[df_selected['PR619 fc'].notnull()]
df_tak.reset_index(inplace=True)
df_mg.reset_index(inplace=True)
df_pr.reset_index(inplace=True)

p_cols_x = ['TAK243 p_x', 'MG132 p_x', 'PR619 p_x']

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
                    className='text-center col-12'),
        ]),
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
                ),
                dbc.Alert(
                    html.H6("Gene Name found", className='text-center'),
                    className='col-3 position-relative"',
                    is_open=False,
                    fade=True,
                    color="success",
                    dismissable=True,
                    id='success',
                    duration=2000,
                )
            ])
        ]),
    ], style={"height": "40px"}),
    # Sliders and Volcano
    dbc.Row([
        # TAK
        dbc.Col(
            html.Div([
                html.H4(children='Conjugation Inhibition',
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
                    max=max(df_tak['TAK243 p']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(0, math.ceil(max(df_tak['TAK243 p'])))},
                    value=1.3
                ),
                # Fold Change slider TAK
                html.H5('Fold Change'),
                dcc.RangeSlider(
                    tooltip={'always_visible': True},
                    id='tak-volcanoplot-input_fc',
                    min=min(df_tak['TAK243 fc']),
                    max=max(df_tak['TAK243 fc']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(math.floor(min(df_tak['TAK243 fc'])),
                                                               math.ceil(max(df_tak['TAK243 fc'])))},
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
                html.H4(children='Proteasome Inhibition',
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
                    max=max(df_mg['MG132 p']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(0, math.ceil(max(df_mg['MG132 p'])))},
                    value=1.3
                ),
                # Fold Change slider MG
                html.H5('Fold Change'),
                dcc.RangeSlider(
                    tooltip={'always_visible': True},
                    id='mg-volcanoplot-input_fc',
                    min=min(df_mg['MG132 fc']),
                    max=max(df_mg['MG132 fc']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(math.floor(min(df_mg['MG132 fc'])),
                                                               math.ceil(max(df_mg['MG132 fc'])))},
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
                html.H4(children='DUB Inhibition',
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
                    max=max(df_pr['PR619 p']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(0, math.ceil(max(df_pr['PR619 p'])))},
                    value=1.3
                ),
                # Fold Change slider PR
                html.H5('Fold Change'),
                dcc.RangeSlider(
                    tooltip={'always_visible': True},
                    id='pr-volcanoplot-input_fc',
                    min=min(df_pr['PR619 fc']),
                    max=max(df_pr['PR619 fc']),
                    step=0.1,
                    marks={i: {'label': str(i)} for i in range(math.floor(min(df_pr['PR619 fc'])),
                                                               math.ceil(max(df_pr['PR619 fc'])))},
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
                    columns=[{"name": i, "id": i} for i in ['Lysine', 'Gene', 'TAK243 p', 'TAK243 fc', 'MG132 p',
                                                            'MG132 fc', 'PR619 p', 'PR619 fc', 'Score', 'Pos. Score',
                                                            'Sequence']],
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
                ),
            ])
        ])
    ]),
    # References
    dbc.Row([
        dbc.Col([
            html.A('Github repository',
                     href='https://github.com/Trulsson/Dynamic-ub-site-dashboard'),
        ]),
        dbc.Col([
            html.A('Trulsson et al. Nat Commun 13. (2022)',
                     href='https://www.nature.com/articles/s41467-022-30376-7'),
        ]),
    ]), html.Br(),
],className="m-auto")


@app.callback(
    Output('err', 'is_open'),
    Output('success', 'is_open'),
    Output('gene_input', 'value'),
    Input('gene_input', 'value'),
    State("err", "is_open")
)
def check_input(input_gene, is_open):
    # Makes sure there is an input and then checks if
    # the input Gene of interest is present in the data
    if input_gene is None:
        return no_update, no_update, no_update
    print("Search for Gene: ", input_gene)
    poi = df_selected[df_selected['Gene'].str.contains(input_gene.upper()) == True]
    if len(poi.index) == 0:
        is_open = True
        return is_open, no_update, no_update
    is_open = True
    return no_update, is_open, no_update


@app.callback(
    Output('ubisite-tak-volcanoplot', 'figure'),
    Output('ubisite-mg-volcanoplot', 'figure'),
    Output('ubisite-pr-volcanoplot', 'figure'),
    Input('tak-volcanoplot-input_p', 'value'),
    Input('tak-volcanoplot-input_fc', 'value'),
    Input('mg-volcanoplot-input_p', 'value'),
    Input('mg-volcanoplot-input_fc', 'value'),
    Input('pr-volcanoplot-input_p', 'value'),
    Input('pr-volcanoplot-input_fc', 'value'),
    Input('gene_input', 'value')
)
def update_volcanoplot_tak(p_val_tak, fc_tak, p_val_mg, fc_mg, p_val_pr, fc_pr, input_gene):
    p_val = [p_val_tak, p_val_mg, p_val_pr]
    fc = [fc_tak, fc_mg, fc_pr]
    # No NaN values in volcanoes
    dfs = [df_tak, df_mg, df_pr]
    figs = []
    for p, f, data, p_col_x, fc_col, p_col in zip(p_val, fc, dfs, p_cols_x, fc_cols, p_cols):
        fig = dashbio.VolcanoPlot(
            dataframe=data,
            effect_size=fc_col,
            p=p_col_x,
            gene='gene+pos',
            point_size=8,
            xlabel='Log2 Difference',
            effect_size_line=f,
            genomewideline_value=p,
            highlight_color='#119DFF',
            col='#2A3F5F',
            snp='Score',
        )
        # Hover data information, order of stack not the same as plot order. Not functional atm, default hover enabled.
        # single_gene_names = []
        # for string in :
        #     single_gene_names.append(string.split(';', 1)[0])
        # customdata = np.stack((data['gene+pos'], data['Score'], data[p_col_x], data[fc_col]), axis=-1)
        # hovertemplate = '<b>Gene: %{customdata[0]}</b> <br>Fold Change: %{customdata[3]} <br>p-value: %{customdata[2]} <br>Score:%{customdata[1]}'
        # fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)

        # Figure layout
        fig.update_layout(showlegend=False, title_text=None,
                          margin={'l': 10, 'r': 10, 'b': 30, 't': 30, 'pad': 4},
                          )

        # Update graph with Protein of interest
        if input_gene:
            p_vals_poi = data[p_col].loc[data['Gene'].str.contains(input_gene.upper()) == True]
            folds_poi = data[fc_col].loc[data['Gene'].str.contains(input_gene.upper()) == True]
            names = data['gene+pos'].loc[data['Gene'].str.contains(input_gene.upper()) == True]

            texts = []
            for fold_poi, p_poi, name in zip(folds_poi, p_vals_poi, names):
                texts.append(fig.add_annotation(x=fold_poi, y=p_poi, text=name, bgcolor=(
                    '#69b8f0' if (p_poi > p) and any([(fold_poi < f[0]), (fold_poi > f[1])]) else 'white')))
        figs.append(fig)
    return figs[0], figs[1], figs[2]


@app.callback(
    Output('poi-table', "data"),
    Input('poi-table', "sort_by"),
    Input('gene_input', 'value'))
def update_table(sort_by, input_gene):
    if len(sort_by) and input_gene:
        poi = df_selected[df_selected['Gene'].str.contains(input_gene.upper()) == True].sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )
    elif input_gene:
        # No sort is applied
        poi = df_selected[df_selected['Gene'].str.contains(input_gene.upper()) == True]

    else:
        return

    return poi.to_dict('records')


if __name__ == '__main__':
    app.run_server(debug=True)
