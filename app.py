from typing import List, Tuple

import pandas as pd
import peptacular.constants
import streamlit as st
from peptacular.fragment import build_fragments, Fragment
from peptacular.score import hyper_score, binomial_score, compute_fragment_matches, FragmentMatch
from peptacular.sequence import strip_modifications, get_modifications
import plotly.graph_objects as go
import plotly.express as px

from constants import *
from utils import color_by_ion_type, COLOR_DICT, is_float, get_fragment_color

# Parse query parameters
params = st.query_params
query_peptide_sequence = params.get('sequence', DEFAULT_PEPTIDE)
query_min_charge = int(params.get('min_charge', DEFAULT_MIN_CHARGE))
query_max_charge = int(params.get('max_charge', DEFAULT_MAX_CHARGE))
query_mass_type = params.get('mass_type', DEFAULT_MASS_TYPE)
query_fragment_types = list(params.get('fragment_types', DEFAULT_FRAGMENT_TYPES))

query_spectra = params.get('spectra', DEFAULT_SPECTRA)
query_spectra = '\n'.join([f'{pair.split(":")[0]} {pair.split(":")[1]}' for pair in query_spectra.split(';')])

st.set_page_config(page_title="peptidefragmenter", page_icon=":bomb:", layout="wide")

# Sidebar: Peptide Fragmenter input
with st.sidebar:
    st.title('Peptide Fragmenter :bomb:')
    st.markdown(
        """This app takes an amino acid sequence and calculates the fragment ions for a given charge range. 
        Modifications should be provided in parentheses with the mass difference in Daltons. Terminal modifications 
        use square brackets.""")

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
                         index=0 if query_mass_type == 'monoisotopic' else 1)
    is_monoisotopic = mass_type == 'monoisotopic'

    fragment_types = st.multiselect(label='Fragment Types',
                                    options=['a', 'b', 'c', 'x', 'y', 'z'],
                                    default=query_fragment_types,
                                    help='Fragment types to calculate')

    internal_fragments = st.checkbox(label='Internal Fragments',
                                     value=False,
                                     help='Include internal fragments')


def generate_app_url(sequence: str, min_charge: int, max_charge: int, mass_type: str, fragment_types: List[str]):
    # Generate the app URL
    url = f'{BASE_URL}?sequence={sequence}&min_charge={min_charge}&max_charge={max_charge}&mass_type={mass_type}&fragment_types={"".join(fragment_types)}'
    return url


url = generate_app_url(peptide_sequence, min_charge, max_charge, mass_type, fragment_types)

st.write(f'##### [Analysis URL]({url}) (copy me and send to your friends!)')

t1, t2, t3, t4 = st.tabs(['Results', 'Spectra', 'Wiki', 'Help'])


@st.cache_data
def create_fragment_table(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool,
                          internal: bool) -> Tuple[List, pd.DataFrame]:
    fragments = build_fragments(sequence=sequence,
                                ion_types=ion_types,
                                charges=charges,
                                monoisotopic=monoisotopic,
                                internal=internal)

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
                                           internal=internal_fragments)

frag_df_downloaded = frag_df.to_csv(index=False)

# make a plotly plot that will graph the segments end -> start on the y-axis, and mass on the x-axis
traces = []
for idx, row in frag_df[frag_df['internal'] == False].iterrows():
    traces.append(
        go.Scatter(
            x=[row['mz'], row['mz']],
            y=[row['start'], row['end']],
            mode='lines',
            line=dict(color=get_fragment_color(row)),
            name=row['label'],
        )
    )

# Create layout for the plot
layout = go.Layout(
    title="Fragment Segments",
    xaxis=dict(title='M/Z'),
    yaxis=dict(title='Sequence'),
    showlegend=False
)

# Create a Figure and add the traces
fig = go.Figure(data=traces, layout=layout)
fig.update_yaxes(ticktext=list(unmodified_sequence), tickvals=list(range(len(unmodified_sequence))))

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
    df['# (abc)'] = list(range(1, len(df) + 1))
    df['# (xyz)'] = list(range(1, len(df) + 1))[::-1]

    # reorder columns so that # is first # +1 is last and AA is in the middle
    df = df[['AA'] + ['# (abc)'] + [col for col in df.columns if col not in ['AA', '# (abc)', '# (xyz)']] + ['# (xyz)']]
    dfs.append(df)

combined_df = pd.DataFrame(combined_data)
# sort columns based on alphabetical order
combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

styled_dfs = []

for df in dfs:
    styled_df = df.style.apply(color_by_ion_type)

    # Set table styles with increased horizontal padding for more space between columns,
    # centered text, and no borders
    styles = [
        dict(selector="td", props=[("padding", "2px 2px"), ("text-align", "center"), ("border", "none")]),
        dict(selector="th", props=[("padding", "2px 2px"), ("text-align", "center"), ("border", "none")])
    ]
    styled_df = styled_df.set_table_styles(styles)
    styled_dfs.append(styled_df)

# CSS to inject contained in a string
hide_table_row_index_and_adjust_padding = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            td, th {padding: 0px}  /* Padding adjustment for all table cells */
            </style>
            """

# Inject CSS with Markdown
st.markdown(hide_table_row_index_and_adjust_padding, unsafe_allow_html=True)

with t1:
    st.header('Fragment Ions')

    for styled_df, charge in zip(styled_dfs, list(range(min_charge, max_charge + 1))):
        st.subheader(f'Charge {charge}')
        st.table(styled_df)

    st.plotly_chart(fig, use_container_width=True)

    with st.expander('Fragment Ion Data'):
        st.dataframe(frag_df, use_container_width=True)
        st.download_button(label='Download CSV', data=frag_df_downloaded, file_name='fragment_ions.csv',
                           use_container_width=True)

with t2:
    st.header('Input Spectra')
    st.caption('Add spectra to match fragment ions to. One per line. Format: {m/z} {intensity}')

    c1, c2 = st.columns(2)
    tolerance_type = c2.radio(label='Tolerance Type',
                              options=TOLERANCE_OPTIONS,
                              index=DEFAULT_TOLERANCE_TYPE_INDEX,
                              help='Offset type to add to spectra')

    tolerance = c1.number_input(label='Tolerance',
                                value=DEFAULT_TOLERANCE_TH if tolerance_type == 'th' else DEFAULT_TOLERANCE_PPM,
                                step=TOLERANCE_STEP_TH if tolerance_type == 'th' else TOLERANCE_STEP_PPM,
                                min_value=MIN_TOLERANCE_VALUE,
                                max_value=MAX_TOLERANCE_VALUE_TH if tolerance_type == 'th' else MAX_TOLERANCE_VALUE_PPM,
                                help='Tolerance to use when matching fragment ions to spectra')

    min_intensity = st.number_input(label='Min Intensity',
                                    value=DEFAULT_MIN_INTENSITY,
                                    step=1.0,
                                    min_value=0.0)
    spectra = st.text_area(label='Spectra',
                           value=query_spectra,
                           help='Spectra to match fragment ions to. One per line. Format: {m/z} {intensity}\\n',
                           max_chars=30_000)

    if spectra:

        mz_values, intensity_values = [], []

        for line in spectra.split('\n'):
            mz, intensity = line.split(' ')
            mz = float(mz)
            intensity = float(intensity)

            if intensity <= min_intensity:
                continue

            mz_values.append(mz)
            intensity_values.append(intensity)

        max_spectra_mz = max(mz_values)

        fragment_matches = compute_fragment_matches(fragments, mz_values, intensity_values, tolerance, tolerance_type)
        fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)
        fragment_matches = {fm.mz: fm for fm in fragment_matches}  # keep the best error for each fragment

        data = []

        for mz, i in zip(mz_values, intensity_values):
            fm = fragment_matches.get(mz, None)

            if fm:
                data.append(
                    {'sequence': fm.fragment.sequence, 'charge': fm.fragment.charge, 'ion_type': fm.fragment.ion_type,
                     'number': fm.fragment.number, 'internal': fm.fragment.internal,
                     'parent_number': fm.fragment.parent_number, 'monoisotopic': fm.fragment.monoisotopic, 'mz': mz,
                     'intensity': i, 'error': fm.error, 'abs_error': abs(fm.error)})

            else:
                data.append({'sequence': '', 'charge': 0, 'ion_type': '', 'number': 0, 'internal': False,
                             'parent_number': 0, 'monoisotopic': True, 'mz': mz,
                             'intensity': i, 'error': 0, 'abs_error': 0})

        spectra_df = pd.DataFrame(data)

        # for keep only the lowest abs_error for ion_type, charge, num
        spectra_df.sort_values(by='abs_error', inplace=True)

        spectra_df['ion_color_type'] = spectra_df['ion_type']
        spectra_df.loc[spectra_df['internal'] == True, 'ion_color_type'] = 'i'

        ion_labels = []
        for _, row in spectra_df.iterrows():

            try:
                charge_str = '+' * int(row['charge'])
                ion_type_str = row['ion_type']
                parent_number_str = str(int(row['parent_number']))
            except ValueError:
                charge_str = ''
                ion_type_str = ''
                parent_number_str = ''

            ion_labels.append(f"{charge_str}{ion_type_str}{parent_number_str}")

        spectra_df['ion_label'] = ion_labels
        spectra_df.loc[spectra_df['internal'] == True, 'ion_label'] += 'i'

        COLOR_DICT.setdefault('', 'grey')
        fig = px.bar(spectra_df, x='mz', y='intensity', color='ion_color_type',
                     hover_data=['charge', 'error', 'sequence'],
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

        spectra_df.sort_values(by='mz', inplace=True)

        st.caption('Score are under development and may not be accurate')
        hs = hyper_score(fragments, spectra_df['mz'].tolist(), spectra_df['intensity'].tolist(), tolerance,
                         tolerance_type)
        st.metric(f'Hyperscore', hs)
        bs = binomial_score(fragments, spectra_df['mz'].tolist(), spectra_df['intensity'].tolist(), tolerance,
                            tolerance_type)
        st.metric(f'Binomial Score', bs)


        def highlight_cells(data):
            # Initialize empty DataFrame with same index and columns as original
            styled = pd.DataFrame('', index=data.index, columns=data.columns)

            # Iterate over cells and update `styled` based on cell position
            for row in data.index:
                for col in data.columns:
                    if col == 'AA':
                        continue
                    label = '+' * int(col[1:]) + col[0] + str(row + 1)
                    if label in accepted_normal_ions:
                        styled.loc[row, col] = 'background-color: yellow'
                    elif label + 'i' in accepted_internal_ions:
                        styled.loc[row, col] = 'background-color: magenta'

            return styled


        matched_ions = spectra_df[spectra_df['ion_type'] != '']
        accepted_normal_ions = matched_ions[matched_ions['internal'] == False]['ion_label'].tolist()
        accepted_internal_ions = matched_ions[matched_ions['internal'] == True]['ion_label'].tolist()

        combined_df = combined_df.style.apply(highlight_cells, axis=None)
        st.table(combined_df)

        with st.expander('Spectra Data'):
            spectra_df.sort_values(by=['mz'], inplace=True)
            st.dataframe(spectra_df)
            st.download_button(label='Download CSV', data=spectra_df.to_csv(index=False),
                               file_name='spectra_results.csv',
                               use_container_width=True)

    with t3:
        st.markdown(WIKI)

    with t4:
        st.markdown(HELP)
