import json
from typing import List, Tuple

import pandas as pd
import streamlit as st
import peptacular as pt
import plotly.graph_objects as go

from constants import *
from utils import COLOR_DICT, get_fragment_color

# Parse query parameters
params = st.query_params
query_peptide_sequence = params.get('sequence', DEFAULT_PEPTIDE)
query_min_charge = int(params.get('min_charge', DEFAULT_MIN_CHARGE))
query_max_charge = int(params.get('max_charge', DEFAULT_MAX_CHARGE))
query_mass_type = params.get('mass_type', DEFAULT_MASS_TYPE)
query_fragment_types = list(params.get('fragment_types', DEFAULT_FRAGMENT_TYPES).split(';'))

# Initialize session state for fragment types if it doesn't exist (page refresh)
if 'fragment_types' not in st.session_state:
    st.session_state.fragment_types = query_fragment_types

st.set_page_config(page_title="peptidefragmenter", page_icon=":bomb:", layout="wide")

# Sidebar: Peptide Fragmenter input
with st.sidebar:
    st.title('Peptide Fragmenter :bomb:')
    st.caption("""
    A simple peptide fragment ion claculator. Now ProForma 2.0 compliant!
    
    See the help tab for more information on supported ProForma features.
    """)

    st.caption('Made with [peptacular](https://pypi.org/project/peptacular/)')

    peptide_sequence = st.text_input('Peptide Sequence',
                                     value=query_peptide_sequence,
                                     max_chars=MAX_PEPTIDE_LENGTH,
                                     help='Peptide sequence to fragment. Include modifications in square brackets.').strip()

    annotation = pt.parse(peptide_sequence)

    apply_carbamidomethyl = st.checkbox('Use carbamidomethyl', value=False, help='Apply carbamidomethyl mod')

    if apply_carbamidomethyl is True:
        static_mod = '[Carbamidomethyl]@C'
        annotation.add_static_mods(static_mod)

    c1, c2 = st.columns(2)
    c1.caption(f'Residues: {len(annotation)}/{MAX_PEPTIDE_AA_COUNT}')

    # Check peptide AA count is within limits
    if len(annotation) > MAX_PEPTIDE_AA_COUNT:
        st.error(f'Peptide length cannot exceed {MAX_PEPTIDE_AA_COUNT} amino acids')
        st.stop()

    try:
        neutral_sequence_mass = pt.mass(annotation, monoisotopic=True, ion_type='p', charge=0)
    except Exception as e:
        st.error(f'Error calculating peptide mass: {e}')
        st.stop()

    c2.caption(f'Neutral Mass: {neutral_sequence_mass:.5f}')

    if annotation.contains_sequence_ambiguity() or annotation.contains_residue_ambiguity() or annotation.contains_mass_ambiguity():
        st.error('Sequence cannot contain ambiguity!')
        st.stop()

    if len(annotation) == 0:
        st.error('Peptide sequence cannot be empty')
        st.stop()

    c1, c2 = st.columns(2)
    min_charge = c1.number_input(label='Min Charge',
                                 min_value=MIN_PEPTIDE_CHARGE,
                                 max_value=MAX_PEPTIDE_CHARGE,
                                 value=query_min_charge,
                                 help='Minimum charge to fragment')
    max_charge = c2.number_input(label='Max Charge',
                                 min_value=MIN_PEPTIDE_CHARGE,
                                 max_value=MAX_PEPTIDE_CHARGE,
                                 value=query_max_charge,
                                 help='Maximum charge to fragment')

    # verify min charge is less or equal to than max charge
    if min_charge > max_charge:
        st.error('Min charge must be less than or equal to max charge')
        st.stop()

    mass_type = st.radio(label='Mass Type',
                         options=['monoisotopic', 'average'],
                         help='Mass type to use for fragment calculation',
                         index=0 if query_mass_type == 'monoisotopic' else 1,
                         horizontal=True)
    is_monoisotopic = mass_type == 'monoisotopic'

    st.caption('Terminal Ions')

    c1, c2, c3 = st.columns(3)
    a, b, c = c1.checkbox('a', value='a' in st.session_state.fragment_types), \
        c2.checkbox('b', value='b' in st.session_state.fragment_types), \
        c3.checkbox('c', value='c' in st.session_state.fragment_types)
    c1, c2, c3 = st.columns(3)
    x, y, z = c1.checkbox('x', value='x' in st.session_state.fragment_types), \
        c2.checkbox('y', value='y' in st.session_state.fragment_types), \
        c3.checkbox('z', value='z' in st.session_state.fragment_types)

    internal_ions = st.checkbox('Internal Ions', value=False)

    ax, ay, az = False, False, False
    bx, by, bz = False, False, False
    cx, cy, cz = False, False, False
    if internal_ions:
        c1, c2, c3 = st.columns(3)
        ax, ay, az = c1.checkbox('ax', value='ax' in st.session_state.fragment_types), \
            c2.checkbox('ay', value='ay' in st.session_state.fragment_types), \
            c3.checkbox('az', value='az' in st.session_state.fragment_types)
        c1, c2, c3 = st.columns(3)
        bx, by, bz = c1.checkbox('bx', value='bx' in st.session_state.fragment_types), \
            c2.checkbox('by', value='by' in st.session_state.fragment_types), \
            c3.checkbox('bz', value='bz' in st.session_state.fragment_types)
        c1, c2, c3 = st.columns(3)
        cx, cy, cz = c1.checkbox('cx', value='cx' in st.session_state.fragment_types), \
            c2.checkbox('cy', value='cy' in st.session_state.fragment_types), \
            c3.checkbox('cz', value='cz' in st.session_state.fragment_types)

    c1, c2 = st.columns(2)
    deselected_all = c1.button(label='Deselect All', use_container_width=True)
    select_all = c2.button(label='Select All', use_container_width=True)

    if deselected_all:
        st.session_state.fragment_types = []
        st.experimental_rerun()
    if select_all:
        st.session_state.fragment_types = list(set(st.session_state.fragment_types + ['a', 'b', 'c', 'x', 'y', 'z']))
        st.experimental_rerun()

    fragment_types = []

    if a:
        fragment_types.append('a')

    if b:
        fragment_types.append('b')

    if c:
        fragment_types.append('c')

    if x:
        fragment_types.append('x')

    if y:
        fragment_types.append('y')

    if z:
        fragment_types.append('z')

    if ax:
        fragment_types.append('ax')

    if ay:
        fragment_types.append('ay')

    if az:
        fragment_types.append('az')

    if bx:
        fragment_types.append('bx')

    if by:
        fragment_types.append('by')

    if bz:
        fragment_types.append('bz')

    if cx:
        fragment_types.append('bx')

    if cy:
        fragment_types.append('cy')

    if cz:
        fragment_types.append('cz')

    # update session state
    st.session_state.fragment_types = fragment_types


def generate_app_url(sequence: str, min_charge: int, max_charge: int, mass_type: str, fragment_types: List[str]):
    # Generate the app URL
    url = f'{BASE_URL}?sequence={sequence}&min_charge={min_charge}&max_charge={max_charge}&mass_type={mass_type}&fragment_types={";".join(fragment_types)}'
    return url


url = generate_app_url(annotation.serialize(), min_charge, max_charge, mass_type, fragment_types)

# set query params

st.write(f'##### :link: [Sharable URL]({url})')

t1, t3, t4 = st.tabs(['Results', 'Wiki', 'Help'])

@st.cache_data
def create_fragment_table(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool) -> Tuple[
    List, pd.DataFrame]:
    fragments = pt.fragment(sequence=sequence,
                         ion_types=ion_types,
                         charges=charges,
                         monoisotopic=monoisotopic)

    seq_len = pt.sequence_length(sequence)

    # convert list of dataclasses to list of dicts
    frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])

    frag_df['number'] = None
    # for forward ions (a,b,c) set number to frag.end
    frag_df['number'] = frag_df.apply(lambda row: row['end'] if row['ion_type'] in 'abc' else None, axis=1)

    # for reverse ions (x,y,z) set number to frag.start
    frag_df['number'] = frag_df.apply(lambda row: seq_len - row['start'] if row['ion_type'] in 'xyz' else row['number'],
                                      axis=1)

    return fragments, frag_df


# Get all fragment ions
fragments, frag_df = create_fragment_table(sequence=annotation.serialize(),
                                           ion_types=fragment_types,
                                           charges=list(range(min_charge, max_charge + 1)),
                                           monoisotopic=is_monoisotopic)

if len(frag_df) == 0:
    st.warning('No fragments found. Please check your input and try again.')
    st.stop()

# drop isotope, loss and parent_sequence columns
frag_df = frag_df.drop(columns=['isotope', 'loss', 'parent_sequence'])

# sort by ion_type, then start
frag_df.sort_values(by=['charge', 'ion_type', 'start'], inplace=True)

# drop any duplicates
frag_df.drop_duplicates(subset=['charge', 'ion_type', 'start', 'end'], inplace=True)

frag_df_downloaded = frag_df.to_csv(index=False)

traces = []
seen = set()
for idx, row in frag_df.iterrows():

    # Determine the Scatter object based on the condition
    if row['ion_type'] in {'a', 'b', 'c'}:
        scatter = go.Scatter(
            x=[row['mz'], row['mz']],
            y=[row['start'], row['end']],
            mode='lines',
            line=dict(color=get_fragment_color(row)),
            name=row['ion_type'],
            legendgroup=row['ion_type'],
            showlegend=row['ion_type'] not in seen

        )
    elif row['ion_type'] in {'x', 'y', 'z'}:
        scatter = go.Scatter(
            x=[row['mz'], row['mz']],
            y=[row['start'] + 1, row['end'] + 1],
            mode='lines',
            line=dict(color=get_fragment_color(row)),
            name=row['ion_type'],
            legendgroup=row['ion_type'],
            showlegend=row['ion_type'] not in seen
        )
    else:
        scatter = go.Scatter(
            x=[row['mz'], row['mz']],
            y=[row['start'] + 1, row['end'] + 1],
            mode='lines',
            line=dict(color=get_fragment_color(row)),
            name=row['ion_type'],
            legendgroup=row['ion_type'],
            showlegend=row['ion_type'] not in seen
        )

    seen.add(row['ion_type'])

    # Append the Scatter object to the traces list
    traces.append(scatter)

# Assuming traces is a list of go.Scatter objects
min_x = min(trace['x'][0] for trace in traces)  # Find the smallest x-value
max_x = max(trace['x'][1] for trace in traces)  # Find the largest x-value

# Expand the x-axis range a bit
padding = (max_x - min_x) * 0.01  # 1% padding on each side
x_range = [min_x - padding, max_x + padding]

# Create layout for the plot with updated x-axis range
layout = go.Layout(
    title={
        'text': 'Fragment Segments',
        'x': 0.5,
        'xanchor': 'center',
        'font': {'size': 20}
    },
    xaxis=dict(title='M/Z', range=x_range),
    yaxis=dict(title='Sequence'),
    showlegend=True
)

components = list(annotation.sequence)
ion_graph = go.Figure(data=traces, layout=layout)
ion_graph.update_yaxes(ticktext=[''] + components + [''],
                       tickvals=list(range(len(annotation) + 2)))

dfs = []
jsons = []

combined_data = {'AA': list(components)}
for charge in range(min_charge, max_charge + 1):
    data = {'AA': components}
    json_data = {'AA': components, 'info': {'charge': charge, 'sequence': annotation.serialize(), 'mass_type': mass_type}}
    for ion_type in sorted(fragment_types):

        if ion_type not in {'a', 'b', 'c', 'x', 'y', 'z'}:
            continue

        ion_df = frag_df[
            (frag_df['ion_type'] == ion_type) & (frag_df['charge'] == charge)]
        ion_df.sort_values(by=['number'], inplace=True)

        frags = ion_df['mz'].tolist()

        json_data[ion_type] = frags

        if ion_type in 'xyz':
            frags = frags[::-1]

        data[ion_type] = frags

        combined_data[ion_type + str(charge)] = frags

    jsons.append(json.dumps(json_data))

    # Displaying the table
    df = pd.DataFrame(data)
    df['+#'] = list(range(1, len(df) + 1))
    df['-#'] = list(range(1, len(df) + 1))[::-1]

    # reorder columns so that # is first # +1 is last and AA is in the middle
    combined_cols = df.columns.tolist()
    combined_cols.remove('+#')
    combined_cols.remove('-#')
    combined_cols.remove('AA')
    forward_cols = [col for col in combined_cols if 'a' in col or 'b' in col or 'c' in col]
    reverse_cols = [col for col in combined_cols if 'x' in col or 'y' in col or 'z' in col]

    # sort
    forward_cols.sort()
    reverse_cols.sort(reverse=True)

    new_cols = ['+#'] + forward_cols + ['AA'] + reverse_cols + ['-#']
    df = df[new_cols]

    # reorder columns so that # is first # +1 is last and AA is in the middle
    dfs.append(df)


# Function to generate color array for cell text based on ion types
def generate_text_color_array(df, color_dict):
    # Initialize a list of lists (each sublist corresponds to a column in df)
    color_array = []
    for col in df.columns:
        if col in color_dict:
            # Apply the specific color for the column
            color_array.append([color_dict[col]] * len(df))
        else:
            # Default color
            color_array.append(['black'] * len(df))
    return color_array


# Function to generate background color array for cells
def generate_bg_color_array(df):
    bg_color_array = []
    for col in df.columns:
        if df[col].dtype.kind in 'f':  # Assuming fragment ion masses are stored in float columns
            # Set background color to white for fragment ion mass columns
            bg_color_array.append(['white'] * len(df))
        else:
            # Default background color for other columns
            bg_color_array.append(['white'] * len(df))
    return bg_color_array


number_of_rows = len(dfs[0])  # Assuming dfs[0] is your DataFrame
cell_height = 25  # As set in the cells dictionary
header_height = 30  # An estimate, adjust based on your header size
margin = 80  # Adjust as needed for title, padding, etc.

# Calculate total figure height
total_height = (number_of_rows * cell_height) + header_height + margin

# Now, create Plotly tables for each DataFrame in dfs
figs = []
for df, charge in zip(dfs, range(min_charge, max_charge + 1)):
    # Format numeric values to 3 decimal places
    formatted_cell_values = []
    for col in df.columns:
        if df[col].dtype.kind in 'f':  # Check if column is float or integer
            formatted_cell_values.append(df[col].map('{:.5f}'.format).tolist())  # Format to 3 decimal places
        else:
            formatted_cell_values.append(df[col].tolist())

    header_values = list(df.columns)

    # capitalize the ion types:
    header_values = [col.upper() if col in COLOR_DICT else col for col in header_values]

    cell_values = formatted_cell_values
    text_colors = generate_text_color_array(df, COLOR_DICT)
    # Generate the background color array for cells
    bg_colors = generate_bg_color_array(df)

    widths = [100] * len(df.columns)
    widths[0] = 20
    widths[-1] = 20

    fig = go.Figure(data=[go.Table(
        header=dict(values=header_values,
                    fill_color='white',
                    align='center',
                    font=dict(color=text_colors, size=16),
                    height=header_height,
                    line=dict(color='black', width=0)),  # Header font color
        cells=dict(values=cell_values,
                   fill=dict(color=bg_colors),  # Example cell colors, adjust as needed
                   align='center',
                   font=dict(color=text_colors, size=13),
                   height=cell_height,
                   line=dict(color='black', width=0),  # Make header lines invisible
                   ),
        columnwidth=widths,
    )])

    fig.update_layout(
        title={
            'text': f'{annotation.serialize()} | +{charge} | {mass_type}',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 15}
        },
        height=total_height,
        margin=dict(l=80, r=80, t=70, b=10)
    )
    # Centered Text
    figs.append(fig)

combined_df = pd.DataFrame(combined_data)
# sort columns based on alphabetical order
combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

is_monoisotopic = mass_type == 'monoisotopic'
with t1:
    st.markdown(f'<h4><center>{annotation.serialize()}</center></h4>', unsafe_allow_html=True)

    # Mass Table
    # st.markdown('<h5><center>Peptide Mass</center></h5>', unsafe_allow_html=True)

    c1, c2 = st.columns(2)
    for c in range(min_charge, max_charge + 1):

        precursor_mass = pt.mass(sequence=annotation, charge=c, monoisotopic=is_monoisotopic, ion_type='p')
        if c == 0:
            c1.markdown(f'<h6><center>(M)  {precursor_mass:.5f}</center></h6>', unsafe_allow_html=True)
        else:
            c1.markdown(f'<h6><center>(M+{c}H)<sup>+{c}</sup>  {precursor_mass:.5f}</center></h6>',
                        unsafe_allow_html=True)

    for c in range(min_charge, max_charge + 1):
        comp, delta_mass = pt.comp_mass(sequence=annotation, ion_type='p', charge=c)

        chem_formula = ' '.join([f'{e}({comp[e]})' for e in comp])

        if c == 0:
            c2.markdown(f'<h6><center>(M) + ΔMass: {chem_formula} + {delta_mass:.4f} Da</center></h6>', unsafe_allow_html=True)
        else:
            # st.markdown(f'<h5><center>Neutral Composition + Delta Mass</center></h5>', unsafe_allow_html=True)
            c2.markdown(f'<h6><center>(M+{c}H) + ΔMass: {chem_formula} + {delta_mass:.4f} Da</cenrter></h6>', unsafe_allow_html=True)

    st.markdown('---')

    st.markdown('<h5><center>Fragment Ions</center></h5>', unsafe_allow_html=True)
    for c, (fig, j) in enumerate(zip(figs, jsons), min_charge):
        # Set figure layout to adjust height (and width if necessary)
        st.plotly_chart(fig, use_container_width=True)
        st.json(j, expanded=False)
        #st.markdown(fig.to_html(), unsafe_allow_html=True)

    # Create a Figure and add the traces
    st.markdown('---')

    st.plotly_chart(ion_graph, use_container_width=True)

    st.markdown('---')

    # center
    st.markdown('<h5><center>Fragment Data</center></h5>', unsafe_allow_html=True)

    # Display the table
    st.dataframe(frag_df, use_container_width=True)

    # download
    st.download_button(label="Download Fragment Ions", data=frag_df_downloaded, file_name='fragment_ions.csv',
                       mime='text/csv', use_container_width=True)

with t3:
    st.markdown(WIKI)

with t4:
    st.markdown(HELP)
