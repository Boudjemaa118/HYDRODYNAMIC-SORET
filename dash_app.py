import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from db_retriever import DB_retriever
import pandas as pd

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
external_scripts = ['https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML']
 
app = dash.Dash(__name__, 
    external_scripts=external_scripts,
    external_stylesheets=external_stylesheets)

#TODO: add bounding box
#TODO: animation + contour plot + imshow
#TODO: math symbols in labels

db = DB_retriever()
df = db.extract_dataframe({"N": 50})
soret_values = df.soret.unique()
Le_values = df.Le.unique()
m_values = df.m.unique()

decription = '''
We consider the convection flow in infinite horizontal cylinder of square cross-section filled with the isotropic porous material saturated with the binary fluid. The upper and lower boundaries of the cylinder are maintained at different constant temperatures and uniform vertical temperature gradient is imposed on lateral boundaries. All the boundaries are supposed to be rigid and impermeable. The fluid is viscous and incompressible.
'''
equations_text = '''
Full set of equations describing the filtration Soret-driven convection consists of the Darcy-Boussinesq equation, the continuity equation, the energy equation and the concentration evolution equation. The problem for the dimensionless perturbations is written as follows 
$$ 0 = \Delta_2 \psi + \\left( Ra + a \\cos \\omega t  \\right) \\left( \\frac{\partial T}{\partial x}+\\frac{\partial C}{\partial x} \\right),$$ 
$$\\frac{\partial T}{\partial t} + \\frac{\partial \psi}{\partial x} + J(T,\psi) =\Delta_2 T,$$
$$m\\frac{\partial C}{\partial t} +\\varepsilon \\frac{\partial \psi}{\partial x}+J(C,\psi) = Le \\left(\Delta_2 C -\\varepsilon \Delta_2 T \\right). $$
'''
reference = '''
The [paper](https://www.sciencedirect.com/science/article/pii/S1631072111000337) deals with two-dimensional Soret-driven convection in a porous cavity with perfectly conducting boundaries heated from below. It is shown that thermodiffusion effect destroys degeneracy existing in the case of single-component fluid. The scenario of the convection onset is discussed. The boundaries of the diffusive state instability to the small-amplitude and finite-amplitude monotonous and oscillatory perturbations are determined.
'''
app.layout = html.Div(children=[
    html.H1(children='Porous Convection Dashboard', style={'textAlign': 'center'}),
    html.Div(children=[
        dcc.Checklist(
            id='Problem-formulation-Checklist',
            options=[{'label': '', 'value': 'on'}],
            value=['on'],
            style={'width': '2.5%', 'display': 'inline-block', 'transform': 'scale(2)'}
        ),
        html.H3(children='Problem formulation', style={'width': '33%', 'textAlign': 'left', 'display': 'inline-block'})
    ]),
    dcc.Markdown(id='Problem-formulation-textblock',children=decription, style={'display': 'block'}),
    html.Div(children=[
        dcc.Checklist(
            id='Govering-equations-Checklist',
            options=[{'label': '', 'value': 'on'}],
            value=['on'],
            style={'width': '2.5%', 'display': 'inline-block', 'transform': 'scale(2)'}
        ),
        html.H3(children='Govering equations', style={'width': '33%', 'textAlign': 'left', 'display': 'inline-block'})
    ]),
    html.Div(id='Govering-equations-textblock', children=equations_text),
    html.Div(children=[
        dcc.Checklist(
            id='Reference-Checklist',
            options=[{'label': '', 'value': 'on'}],
            value=['on'],
            style={'width': '2.5%', 'display': 'inline-block', 'transform': 'scale(2)'}
        ),
        html.H3(children='Reference', style={'width': '33%', 'textAlign': 'left', 'display': 'inline-block'})
    ]),    
    dcc.Markdown(id='Reference-textblock',children=reference),
    html.H3(children='Parametric study of double-diffusive convection in porous medium with Soret effect'),
    html.Div(children='Figure and table contain the numerical results for nonlinear problem. The dependence of max stream funcition, Nusselt number, asymmetr on managing parameters are shown.'),
    
    html.Div([
        html.Div([
            html.Label('Select \(Le\):'),
            dcc.Dropdown(
                id='Le-dropdown', 
                options=[{'label': Le, 'value': Le} for Le in Le_values], 
                value=Le_values[0],
                multi=True
            )
        ],
        style={'width': '33%', 'display': 'inline-block'}),

        html.Div([
            html.Label('Select \(\\varepsilon\):'),
            dcc.Dropdown(
                id='soret-dropdown',
                options=[{'label': soret, 'value': soret} for soret in soret_values],
                value=soret_values[0],
                multi=True
                )
        ],
        style={'width': '33%',  'float': 'center', 'display': 'inline-block'}),

        html.Div([
            html.Label('Select \(m\):'),
            dcc.Dropdown(
                id='m-dropdown',
                options=[{'label': m, 'value': m} for m in m_values],
                value=m_values[0],
                multi=True
                )
        ],
        style={'width': '33%', 'display': 'inline-block'}),
    ]),

    html.H4(children='Fig 1. Maximum stream function on Rayleight number'),
    dcc.Graph(id='Ra-ampl_psi-graphic'),
    html.H4(children='Table 1. Computational results'),
    html.Div(id='Ra-ampl_psi-table')
])

def create_filtered_dataframe(soret, Le, m):
    if type(m) != list:
        m = [m]
    if type(Le) != list:
        Le = [Le]
    if type(soret) != list:
        soret = [soret]
    return df[(df['soret'].isin(soret)) & (df['Le'].isin(Le)) & (df['m'].isin(m))]

@app.callback(
    Output('Ra-ampl_psi-graphic', 'figure'),
    [Input('soret-dropdown', 'value'),
     Input('Le-dropdown', 'value'),
     Input('m-dropdown', 'value')]
)
def generate_graph(soret, Le, m):
    data = []
    if type(m) != list:
        m = [m]
    if type(Le) != list:
        Le = [Le]
    if type(soret) != list:
        soret = [soret]

    for _m in m:
        for _Le in Le:
            for _soret in soret:
                dff = create_filtered_dataframe(_soret, _Le, _m)
                data.append({
                        'x': dff['Ra'],
                        'y': dff['ampl_psi'],
                        'name': 'eps = {}, Le = {}, m = {}'.format(_soret, _Le, _m)
                    })
    return {
        'data': data,
        'layout': {
            'xaxis': {
                'title': 'Ra',
                #'showline': True,
                'showticklabels': True,
                'zeroline': True,
            },
            'yaxis': {
                'title': 'ampl_psi',
                'showline': True,
                'showticklabels': True,
                'rangemode': 'tozero'
            },
            'legend': {'x': 0, 'y': 1},
            'hovermode': 'closest',
            'transition': {'duration': 500},
        }
    }

@app.callback(
    Output('Ra-ampl_psi-table', 'children'),
    [Input('soret-dropdown', 'value'),
     Input('Le-dropdown', 'value'),
     Input('m-dropdown', 'value')]
)
def generate_table(soret, Le, m, max_rows=100):
    dff = create_filtered_dataframe(soret, Le, m)
    columns = ['soret', 'Le', 'm', 'Ra', 'ampl_psi']
    columns2 = ['soret', '$Le$', 'm', 'Ra', 'ampl_psi']
    return html.Table([
        html.Thead(
            html.Tr([html.Th(children=col,style={'font-style': 'italic'}) for col in columns2])
        ),
        html.Tbody([
            html.Tr([
                html.Td(dff.iloc[i][col]) for col in columns
            ]) for i in range(min(len(dff), max_rows))
        ])
    ])

@app.callback(
   Output('Problem-formulation-textblock', 'style'),
   [Input('Problem-formulation-Checklist', 'value')])
def show_hide_element1(visibility_state):
    if visibility_state:
        return {'display': 'block'}
    return {'display': 'none'}

@app.callback(
   Output('Govering-equations-textblock', 'style'),
   [Input('Govering-equations-Checklist', 'value')])
def show_hide_element2(visibility_state):
    if visibility_state:
        return {'display': 'block'}
    return {'display': 'none'}

@app.callback(
   Output('Reference-textblock', 'style'),
   [Input('Reference-Checklist', 'value')])
def show_hide_element3(visibility_state):
    if visibility_state:
        return {'display': 'block'}
    return {'display': 'none'}

if __name__ == '__main__':
    app.run_server(debug=True)