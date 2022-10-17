from dash import html, dcc
from dash.dependencies import Input, Output

# Connect to main app.py file
from app import app
from app import server

# Connect to your app pages
from apps import dda, dia


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div([
        dcc.Link('DDA data | ', href='/apps/dda'),
        dcc.Link('DIA data', href='/apps/dia'),
    ], className="col-3"),
    html.Div(id='page-content', children=[])
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/apps/dda':
        return dda.layout
    if pathname == '/apps/dia':
        return dia.layout
    else:
        return dda.layout


if __name__ == '__main__':
    app.run_server(debug=False)