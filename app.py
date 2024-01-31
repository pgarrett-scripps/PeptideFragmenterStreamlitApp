from typing import List, Tuple

import pandas as pd
import peptacular.constants
import streamlit as st
from peptacular.fragment import build_fragments
from peptacular.sequence import strip_modifications, get_modifications
import plotly.graph_objects as go

from constants import *
from utils import COLOR_DICT, is_float, get_fragment_color

# Parse query parameters
params = st.query_params
query_peptide_sequence = params.get('sequence', DEFAULT_PEPTIDE)
query_min_charge = int(params.get('min_charge', DEFAULT_MIN_CHARGE))
query_max_charge = int(params.get('max_charge', DEFAULT_MAX_CHARGE))
query_mass_type = params.get('mass_type', DEFAULT_MASS_TYPE)
query_fragment_types = list(params.get('fragment_types', DEFAULT_FRAGMENT_TYPES))

st.set_page_config(page_title="peptidefragmenter", page_icon=":bomb:", layout="wide")

# Sidebar: Peptide Fragmenter input
with st.sidebar:
    st.title('Peptide Fragmenter :bomb:')
    st.markdown(
        """A simple peptide fragment ion claculator. Specify terminal PTMs with [] and internal PTMs with ()."""
    )

    st.markdown('Note that B, X, and Z residues have a mass of 0.0 Da.')

    peptide_sequence = st.text_input('Peptide Sequence',
                                     value=query_peptide_sequence,
                                     max_chars=MAX_PEPTIDE_LENGTH,
                                     help='Peptide sequence to fragment. Include modification masses in parentheses.')
    peptide_len = len(strip_modifications(peptide_sequence))
    st.caption(f'Residues: {peptide_len}/{MAX_PEPTIDE_AA_COUNT}')

    # Check peptide AA count is within limits
    if peptide_len > MAX_PEPTIDE_AA_COUNT:
        st.error(f'Peptide length cannot exceed {MAX_PEPTIDE_AA_COUNT} amino acids')
        st.stop()

    # Verify the input sequence is valid
    unmodified_sequence = strip_modifications(peptide_sequence)
    additional_aa = {'B', 'X', 'Z'}
    valid_aa = additional_aa.union(peptacular.constants.AMINO_ACIDS)
    if not all(valid_aa for aa in unmodified_sequence):
        st.error(f'Invalid amino acid(s) detected.')
        st.stop()

    # verify modifications are valid
    mod_dict = get_modifications(peptide_sequence)
    if not all(is_float(mod) for mod in mod_dict.values()):
        st.error('Invalid modification mass detected.')
        st.stop()

    if peptide_len == 0:
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

    fragment_types = st.multiselect(label='Fragment Types',
                                    options=['a', 'b', 'c', 'x', 'y', 'z'],
                                    default=query_fragment_types,
                                    help='Fragment types to calculate')


def generate_app_url(sequence: str, min_charge: int, max_charge: int, mass_type: str, fragment_types: List[str]):
    # Generate the app URL
    url = f'{BASE_URL}?sequence={sequence}&min_charge={min_charge}&max_charge={max_charge}&mass_type={mass_type}&fragment_types={"".join(fragment_types)}'
    return url


url = generate_app_url(peptide_sequence, min_charge, max_charge, mass_type, fragment_types)

st.write(f'##### [Analysis URL]({url}) (copy me and send to your friends!)')

t1, t3, t4 = st.tabs(['Results', 'Wiki', 'Help'])

@st.cache_data
def create_fragment_table(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool,
                          internal: bool) -> Tuple[List, pd.DataFrame]:
    fragments = build_fragments(sequence=sequence,
                                ion_types=ion_types,
                                charges=charges,
                                monoisotopic=monoisotopic,
                                internal=internal,
                                aa_masses={aa : 0.0 for aa in additional_aa},)

    # convert list of dataclasses to list of dicts
    frag_df = pd.DataFrame([fragment.__dict__ for fragment in fragments])
    frag_df['start'] = [fragment.start for fragment in fragments]
    frag_df['end'] = [fragment.end for fragment in fragments]
    frag_df['label'] = [fragment.label.replace('*', '+') for fragment in fragments]
    frag_df['mz'] = [fragment.mz for fragment in fragments]

    #  for all ends that are negative add seq_len
    frag_df.loc[frag_df['end'] < 0, 'end'] += len(unmodified_sequence)
    frag_df.loc[frag_df['start'] < 0, 'start'] += len(unmodified_sequence)

    # where end is None, set to seq_len
    frag_df.loc[frag_df['end'].isna(), 'end'] = len(unmodified_sequence)
    return fragments, frag_df


# Get all fragment ions
fragments, frag_df = create_fragment_table(sequence=peptide_sequence,
                                           ion_types=fragment_types,
                                           charges=list(range(min_charge, max_charge + 1)),
                                           monoisotopic=is_monoisotopic,
                                           internal=False)

frag_df_downloaded = frag_df.to_csv(index=False)


traces = []
seen = set()
for idx, row in frag_df[frag_df['internal'] == False].iterrows():

    # Determine the Scatter object based on the condition
    if row['ion_type'] in 'abc':
        scatter = go.Scatter(
            x=[row['mz'], row['mz']],
            y=[row['start'], row['end']],
            mode='lines',
            line=dict(color=get_fragment_color(row)),
            name=row['ion_type'],
            legendgroup=row['ion_type'],
            showlegend=row['ion_type'] not in seen

        )
    else:
        scatter = go.Scatter(
            x=[row['mz'], row['mz']],
            y=[row['start']+1, row['end']+1],
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
    title="Fragment Segments",
    xaxis=dict(title='M/Z', range=x_range),
    yaxis=dict(title='Sequence'),
    showlegend=True
)

# Create a Figure and add the traces
fig = go.Figure(data=traces, layout=layout)
fig.update_yaxes(ticktext=['N-Term']+list(unmodified_sequence)+['C-Term'], tickvals=list(range(len(unmodified_sequence)+2)))

dfs = []
combined_data = {'AA': list(unmodified_sequence)}
for charge in range(min_charge, max_charge + 1):
    data = {'AA': list(unmodified_sequence)}
    for ion_type in sorted(fragment_types):
        ion_df = frag_df[
            (frag_df['ion_type'] == ion_type) & (frag_df['charge'] == charge) & (frag_df['internal'] == False)]
        ion_df.sort_values(by=['number'], inplace=True)

        frags = ion_df['mz'].tolist()

        if ion_type in 'xyz':
            frags = frags[::-1]

        data[ion_type] = frags

        combined_data[ion_type + str(charge)] = frags

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

combined_df = pd.DataFrame(combined_data)
# sort columns based on alphabetical order
combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

styled_dfs = []

def highlight_cells(data):
    # Initialize empty DataFrame with same index and columns as original
    styled = pd.DataFrame('', index=data.index, columns=data.columns)

    # Iterate over cells and update `styled` based on cell position
    for row in data.index:
        for col in data.columns:
            if col == 'AA' or col == '+#' or col == '-#':
                styled.loc[
                    row, col] = f'background-color: gainsboro; color: black; text-align: center; font-weight: bold;'
                continue

            styled.loc[
                row, col] = f'color: {COLOR_DICT[col]}; text-align: center;'

    return styled

for df in dfs:
    styled_df = df.style.format(precision=4).apply(highlight_cells, axis=None)
    styled_dfs.append(styled_df)

with t1:
    for styled_df, charge in zip(styled_dfs, list(range(min_charge, max_charge + 1))):
        st.subheader(f'Charge {charge}')
        st.dataframe(styled_df, height=(len(dfs[0]) + 1) * 35 + 3, hide_index=True)

    st.plotly_chart(fig, use_container_width=True)

    frag_df.drop(columns=['parent_number', 'isotope', 'loss', 'aa_masses', 'parent_sequence', 'internal'], inplace=True)

    st.dataframe(frag_df, use_container_width=True)

with t3:
    st.markdown(WIKI)

with t4:
    st.markdown(HELP)
