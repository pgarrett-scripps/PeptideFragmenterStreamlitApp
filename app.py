import pandas as pd
import peptacular.constants
import streamlit as st
from peptacular.sequence import fragment_series, strip_modifications, calculate_mz, parse_modified_sequence
import plotly.graph_objects as go

from utils import color_by_ion_type, COLOR_DICT

params = st.experimental_get_query_params()
PEPTIDE_SEQUENCE = params.get('sequence', ['(-3.14)PEP(123.456)TIDE'])[0]
MIN_CHARGE = int(params.get('min_charge', [1])[0])
MAX_CHARGE = int(params.get('max_charge', [2])[0])
MASS_TYPE = params.get('mass_type', ['monoisotopic'])[0]
FRAGMENT_TYPES = list(params.get('fragment_types', ['by'])[0])

st.set_page_config(page_title="peptidefragmenter", page_icon=":bomb:")

with st.sidebar:
    st.title('Peptide Fragmenter')
    st.markdown("""This app takes an amino acid sequence and calculates the fragment ions for a given charge range. Modifications should be provided in parentheses with the mass difference in Daltons.""")
    st.markdown("""For example, (-3.14)PEP(123.456)TIDE contians a -3.14 N-Term modifiction and a 123.456 modifcation on the second Proline (P).""")

    peptide_sequence = st.text_input('Peptide Sequence',
                                     value=PEPTIDE_SEQUENCE,
                                     help='Peptide sequence to fragment. Include modification masses in parentheses.')

    c1, c2 = st.columns(2)
    min_charge = c1.number_input('Min Charge', min_value=1, value=MIN_CHARGE, help='Minimum charge to fragment')
    max_charge = c2.number_input('Max Charge', min_value=1, value=MAX_CHARGE, help='Maximum charge to fragment')

    if min_charge > max_charge:
        st.error('Min charge must be less than or equal to max charge')
        st.stop()

    mass_type = st.radio('Mass Type',
                             ['monoisotopic', 'average'],
                             help='Mass type to use for fragment calculation',
                             index=0 if MASS_TYPE == 'monoisotopic' else 1)
    is_monoisotopic = mass_type == 'monoisotopic'

    fragment_types = st.multiselect('Fragment Types',
                                    ['a', 'b', 'c', 'x', 'y', 'z'],
                                    default=FRAGMENT_TYPES,
                                    help='Fragment types to calculate')

st.warning('This is a work in progress. Please report any issues or suggestions to pgarrett@scripps.edu.')

# check if anything is invalid
unmodified_sequence = strip_modifications(peptide_sequence)
for aa in unmodified_sequence:
    if aa not in peptacular.constants.AMINO_ACIDS:
        st.error(f'Invalid amino acid: {aa}')
        st.stop()

mod_dict = parse_modified_sequence(peptide_sequence)
for i, mod in mod_dict.items():
    # check if mod can be a float
    try:
        float(mod)
    except ValueError:
        st.error(f'Invalid modification mass: {mod}')
        st.stop()


st.subheader('Results')
data = {'AA': list(unmodified_sequence)}
for ion_type in fragment_types:
    for charge in range(min_charge, max_charge + 1):
        frags = list(fragment_series(sequence=peptide_sequence, ion_type=ion_type,
                                     charge=charge, monoisotopic=is_monoisotopic))
        if ion_type in 'abc':
            frags = frags[::-1]

        data[f'{ion_type}{charge}'] = frags

# Displaying the table
df = pd.DataFrame(data)
df = df.reindex(sorted(df.columns), axis=1)
df_downloaded = df.to_csv(index=False)
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

# Display a static table
st.table(styled_df)
st.download_button(label='Download CSV', data=df_downloaded, file_name='fragment_ions.csv', use_container_width=True)

peptide_mass_data = {}
for charge in range(min_charge, max_charge + 1):
    peptide_mass_data[f'+{charge}'] = [calculate_mz(peptide_sequence,
                                                    charge=charge,
                                                    monoisotopic=is_monoisotopic)]
mass_df = pd.DataFrame(peptide_mass_data)

# Plotting the graph
fig = go.Figure()

stripped_sequence = strip_modifications(peptide_sequence)
forward_seq = [stripped_sequence[:i] for i in range(1, len(stripped_sequence) + 1)]
reverse_seq = [stripped_sequence[i:] for i in range(len(stripped_sequence))]
y_axis = reverse_seq + [''] + forward_seq

df['forward'] = forward_seq
df['forward_ion_num'] = list(range(1, len(df) + 1))
df['reverse'] = reverse_seq
df['reverse_ion_num'] = list(range(1, len(df) + 1))[::-1]

for col in df.columns[1:-4]:  # skip the 'AA' column and the forward/reverse columns
    ion_type = col[0]  # extract ion type from column name
    charge = int(col[1:])  # extract charge from column name
    color = COLOR_DICT.get(ion_type, 'grey')  # get color or default to grey if not found

    if ion_type in 'abc':
        ion_nums = list(range(1, len(df) + 1))
        ion_labels = [f'{ion_type}{i}' for i in ion_nums]
        hover_text = df.apply(lambda row: f"Ion Type: {ion_type}<br>m/z: {row[col]}<br>Charge: {charge}<br>Sequence: "
                                          f"{row['forward']}<br>Ion Num: {row['forward_ion_num']}", axis=1)
        fig.add_trace(go.Bar(x=df[col], y=list(range(1, len(df) + 1)), orientation='v', marker_color=color, name=col,
                             hovertext=hover_text, hoverinfo='text', ))
    else:
        ion_nums = list(range(len(df), 0, -1))
        ion_labels = [f'{ion_type}{i}' for i in ion_nums]
        hover_text = df.apply(lambda row: f"Ion Type: {ion_type}<br>m/z: {row[col]}<br>Charge: {charge}<br>Sequence: "
                                          f"{row['reverse']}<br>Ion Num: {row['reverse_ion_num']}", axis=1)
        fig.add_trace(go.Bar(x=df[col], y=list(range(-len(df), 0)), orientation='v', marker_color=color, name=col,
                             hovertext=hover_text, hoverinfo='text', ))
#make bars wider
fig.update_traces(width=2)
fig.update_layout(barmode='group', yaxis_title="Fragment Ion Sequence", xaxis_title="M/Z", title='Fragment Ions')
fig.update_yaxes(ticktext=y_axis, tickvals=list(range(-len(df), len(df) + 1)))
#make taller
fig.update_layout(height=750)

st.plotly_chart(fig, use_container_width=True)

st.subheader('Peptide Mass/Charge Table')
st.table(mass_df)
