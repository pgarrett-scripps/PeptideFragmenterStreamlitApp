from typing import List

import pandas as pd
import peptacular.constants
import streamlit as st
from peptacular.fragmenter import build_fragments
from peptacular.sequence import strip_modifications, parse_modified_sequence
from sortedcontainers import SortedDict
import plotly.graph_objects as go
import plotly.express as px

from constants import *
from utils import color_by_ion_type, COLOR_DICT, is_float, get_fragment_color


# Parse query parameters
params = st.experimental_get_query_params()
query_peptide_sequence = params.get('sequence', [DEFAULT_PEPTIDE])[0]
query_min_charge = int(params.get('min_charge', [DEFAULT_MIN_CHARGE])[0])
query_max_charge = int(params.get('max_charge', [DEFAULT_MAX_CHARGE])[0])
query_mass_type = params.get('mass_type', [DEFAULT_MASS_TYPE])[0]
query_fragment_types = list(params.get('fragment_types', [DEFAULT_FRAGMENT_TYPES])[0])

st.set_page_config(page_title="peptidefragmenter", page_icon=":bomb:")

# Sidebar: Peptide Fragmenter input
with st.sidebar:
    st.title('Peptide Fragmenter :bomb:')
    st.markdown(
        """This app takes an amino acid sequence and calculates the fragment ions for a given charge range. 
        Modifications should be provided in parentheses with the mass difference in Daltons.""")

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
    if not all(aa in peptacular.constants.AMINO_ACIDS for aa in unmodified_sequence):
        st.error(f'Invalid amino acid(s) detected.')
        st.stop()

    # verify modifications are valid
    mod_dict = parse_modified_sequence(peptide_sequence)
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
                         index=0 if query_mass_type == 'monoisotopic' else 1)
    is_monoisotopic = mass_type == 'monoisotopic'

    fragment_types = st.multiselect(label='Fragment Types',
                                    options=['a', 'b', 'c', 'x', 'y', 'z'],
                                    default=query_fragment_types,
                                    help='Fragment types to calculate')

    internal_fragments = st.checkbox(label='Internal Fragments',
                                     value=False,
                                     help='Include internal fragments')


t1, t2, t3, t4 = st.tabs(['Results', 'Spectra', 'Wiki', 'Help'])


@st.cache_data
def create_fragment_table(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool,
                    internal_fragments: bool) -> pd.DataFrame:

    fragments = build_fragments(sequence=sequence,
                                ion_types=ion_types,
                                charges=charges,
                                monoisotopic=monoisotopic,
                                internal_fragments=internal_fragments)

    # convert list of dataclasses to list of dicts
    frag_df = pd.DataFrame([fragment.__dict__ for fragment in fragments])
    frag_df['start'] = [fragment.start for fragment in fragments]
    frag_df['end'] = [fragment.end for fragment in fragments]
    frag_df['label'] = [fragment.label.replace('*', '+') for fragment in fragments]

    #  for all ends that are negative add seq_len
    frag_df.loc[frag_df['end'] < 0, 'end'] += len(unmodified_sequence)
    frag_df.loc[frag_df['start'] < 0, 'start'] += len(unmodified_sequence)

    # where end is None, set to seq_len
    frag_df.loc[frag_df['end'].isna(), 'end'] = len(unmodified_sequence)
    return frag_df

# Get all fragment ions
frag_df = create_fragment_table(sequence=peptide_sequence,
                            ion_types=fragment_types,
                            charges=list(range(min_charge, max_charge + 1)),
                            monoisotopic=is_monoisotopic,
                            internal_fragments=internal_fragments)

frag_df_downloaded = frag_df.to_csv(index=False)

# make a plotly plot that will graph the segments end -> start on the y-axis, and mass on the x-axis
traces = []
for idx, row in frag_df.iterrows():
    traces.append(
        go.Scatter(
            x=[row['mass'], row['mass']],
            y=[row['start'], row['end']],
            mode='lines',
            line=dict(color=get_fragment_color(row)),
            name=row['label'],
        )
    )

# Create layout for the plot
layout = go.Layout(
    title="Fragment Segments",
    xaxis=dict(title='Mass'),
    yaxis=dict(title='Sequence'),
    showlegend=False
)

# Create a Figure and add the traces
fig = go.Figure(data=traces, layout=layout)
fig.update_yaxes(ticktext=list(unmodified_sequence), tickvals=list(range(len(unmodified_sequence))))

data = {'AA': list(unmodified_sequence)}
for ion_type in sorted(fragment_types):
    for charge in range(min_charge, max_charge + 1):
        ion_df = frag_df[(frag_df['ion_type'] == ion_type) & (frag_df['charge'] == charge) & (frag_df['internal'] == False)]
        ion_df.sort_values(by=['number'], inplace=True)
        frags = ion_df['mass'].tolist()

        if ion_type in 'xyz':
            frags = frags[::-1]

        data[f'{"+"*charge}{ion_type}'] = frags


# Displaying the table
df = pd.DataFrame(data)
styled_df = df.style.apply(color_by_ion_type)

# CSS to inject contained in a string
hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
            """

# Inject CSS with Markdown
st.markdown(hide_table_row_index, unsafe_allow_html=True)

with t1:

    st.subheader('Fragment Ions')
    st.table(styled_df)
    st.plotly_chart(fig)

    with st.expander('Fragment Ion Data'):

        st.dataframe(frag_df, use_container_width=True)
        st.download_button(label='Download CSV', data=frag_df_downloaded, file_name='fragment_ions.csv',
                           use_container_width=True)

with t2:

    st.subheader('Input Spectra')
    st.caption('Add spectra to match fragment ions to. One per line. Format: {m/z} {intensity}')

    c1, c2 = st.columns(2)
    tolerance_type = c2.radio(label='Tolerance Type',
                              options=TOLERANCE_OPTIONS,
                              index=DEFAULT_TOLERANCE_TYPE_INDEX,
                              help='Offset type to add to spectra')

    tolerance = c1.number_input(label='Tolerance',
                                value=DEFAULT_TOLERANCE_TH if tolerance_type == 'Th' else DEFAULT_TOLERANCE_PPM,
                                step=TOLERANCE_STEP_TH if tolerance_type == 'Th' else TOLERANCE_STEP_PPM,
                                min_value=MIN_TOLERANCE_VALUE,
                                max_value=MAX_TOLERANCE_VALUE_TH if tolerance_type == 'Th' else MAX_TOLERANCE_VALUE_PPM,
                                help='Tolerance to use when matching fragment ions to spectra')

    min_intensity = st.number_input(label='Min Intensity',
                                    value=DEFAULT_MIN_INTENSITY,
                                    step=1.0,
                                    min_value=0.0)
    spectra = st.text_area(label='Spectra',
                           value=open('samplespectra.txt').read(),
                           help='Spectra to match fragment ions to. One per line. Format: {m/z} {intensity}\\n',
                           max_chars=30_000)

    if spectra:

        spectra_df = pd.DataFrame()
        spectra_df['mz'] = [float(i.split(' ')[0]) for i in spectra.split('\n')]
        spectra_df['intensity'] = [float(i.split(' ')[1]) for i in spectra.split('\n')]

        max_spectra_mz = spectra_df['mz'].max()

        ions = SortedDict()
        for _, row in frag_df.iterrows():

            if row['mass'] > max_spectra_mz:
                continue

            ions[row['mass']] = row

        # filter by min intensity
        spectra_df = spectra_df[spectra_df['intensity'] >= min_intensity]

        spectral_peak_matches = []
        for mz in spectra_df['mz']:
            mz_lower = (mz - tolerance) if tolerance_type == 'Th' else (mz - tolerance * mz / 1_000_000)
            mz_upper = (mz + tolerance) if tolerance_type == 'Th' else (mz + tolerance * mz / 1_000_000)

            # find the closest ion to the mz
            keys = [k for k in ions.irange(minimum=mz_lower, maximum=mz_upper)]
            if len(keys) > 0:

                key_error = {abs(mz - k): k for k in keys}
                best_key = key_error[min(key_error.keys())]

                closest_row = ions[best_key]

                # calculate the error
                error = (mz - keys[0]) if tolerance_type == 'Th' else (mz - best_key) * 1_000_000 / mz
                closest_row['error'] = error
                spectral_peak_matches.append(closest_row)
            else:
                spectral_peak_matches.append(None)

        # fragment ion info
        spectra_df['ion_type'] = [frag['ion_type'] if frag is not None else '' for frag in spectral_peak_matches]
        spectra_df['internal'] = [frag['internal'] if frag is not None else '' for frag in spectral_peak_matches]
        spectra_df['charge'] = [frag['charge'] if frag is not None else '' for frag in spectral_peak_matches]
        spectra_df['number'] = [frag['number'] if frag is not None else '' for frag in spectral_peak_matches]
        spectra_df['parent_number'] = [frag['parent_number'] if frag is not None else '' for frag in spectral_peak_matches]
        spectra_df['sequence'] = [frag['sequence'] if frag is not None else '' for frag in spectral_peak_matches]

        # mass error
        spectra_df['error'] = [frag['error'] if frag is not None else 0 for frag in spectral_peak_matches]
        spectra_df['abs_error'] = spectra_df['error'].abs()

        # for keep only the lowest abs_error for ion_type, charge, num
        spectra_df.sort_values(by='abs_error', inplace=True)

        # Find duplicates based on 'ion_type', 'charge', 'num', 'internal
        duplicates = spectra_df.duplicated(subset=['ion_type', 'charge', 'number', 'internal'], keep='first')

        spectra_df['ion_color_type'] = spectra_df['ion_type']
        spectra_df.loc[spectra_df['internal'] == True, 'ion_color_type'] = 'i'

        ion_labels = []
        for _, row in spectra_df.iterrows():
            charge_str = '+'*int(row['charge']) if row['charge'] else ''
            ion_labels.append(f"{charge_str}{row['ion_type']}{str(row['parent_number'])}")

        spectra_df['ion_label'] = ion_labels
        spectra_df.loc[spectra_df['internal'] == True, 'ion_label'] += 'i'


        COLOR_DICT.setdefault('', 'grey')
        fig = px.bar(spectra_df, x='mz', y='intensity', color='ion_color_type', hover_data=['charge', 'error', 'sequence'],
                     color_discrete_map=COLOR_DICT)
        fig.update_layout(title='Spectra Plot', xaxis_title='M/Z', yaxis_title='Intensity')

        for i, row in spectra_df.iterrows():
            if row['ion_type']:
                fig.add_annotation(
                    x=row['mz'],
                    y=row['intensity'],
                    text=row['ion_label'],
                    showarrow=False,
                    yshift=10,
                    font=dict(
                        size=13,
                        color=COLOR_DICT[row['ion_color_type']]
                    ),
                )

        st.plotly_chart(fig, use_container_width=True)

        with st.expander('Spectra Data'):
            spectra_df.sort_values(by=['mz'], inplace=True)
            st.dataframe(spectra_df)
            st.download_button(label='Download CSV', data=spectra_df.to_csv(index=False), file_name='spectra_results.csv',
                               use_container_width=True)

    with t3:
        st.markdown(WIKI)

    with t4:
        st.markdown(HELP)
