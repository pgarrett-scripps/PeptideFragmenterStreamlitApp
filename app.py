import pandas as pd
import streamlit as st
from peptacular.sequence import fragment_series, strip_modifications, calculate_mz
import plotly.graph_objects as go

COLOR_DICT = {'a': 'brown', 'b': 'blue', 'c': 'green', 'x': 'orange', 'y': 'red', 'z': 'purple'}


# Function to apply color according to ion type
def color_by_ion_type(col):
    ion_type = col.name[0]
    color = COLOR_DICT.get(ion_type, 'grey')  # get color or default to grey if not found
    return ['color: %s' % color] * len(col)


with st.sidebar:
    st.title('Fragment Ion Calculator')
    st.markdown("""This app takes an amino acid sequence and calculates the fragment ions for a given charge range. Modifications should be provided in parentheses with the mass in Daltons.""")
    st.markdown("""For example, (-3.14)PEP(123.456)TIDE contians a -3.14 N-Term modifaction and a 123.456 modifcation on the second Proline.""")

    peptide_sequence = st.text_input('Peptide Sequence',
                                     value='(-3.14)PEP(123.456)TIDE',
                                     help='Peptide sequence to fragment. Include modification masses in parentheses.')

    min_charge, max_charge = st.slider('Charge Range',
                                       min_value=1,
                                       max_value=10,
                                       value=(1, 2),
                                       help='Charge range to fragment')

    mass_type = st.radio('Mass Type',
                             ['monoisotopic', 'average'],
                             help='Mass type to use for fragment calculation')
    is_monoisotopic = mass_type == 'monoisotopic'

    fragment_types = st.multiselect('Fragment Types',
                                    ['a', 'b', 'c', 'x', 'y', 'z'],
                                    default=['b', 'y'],
                                    help='Fragment types to calculate')

st.warning('This is a work in progress. Please report any issues or suggestions to pgarrett@scripps.edu.')

st.subheader('Fragment Ions')
data = {'AA': list(strip_modifications(peptide_sequence))}
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

peptide_mass_data = {}
for charge in range(min_charge, max_charge + 1):
    peptide_mass_data[f'+{charge}'] = [calculate_mz(peptide_sequence, charge=charge,
                                                               monoisotopic=is_monoisotopic)]
mass_df = pd.DataFrame(peptide_mass_data)
st.subheader('Peptide Mass/Charge Table')
st.table(mass_df)

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

fig.update_layout(barmode='group', yaxis_title="Fragment Ion Sequence", xaxis_title="M/Z")
fig.update_yaxes(ticktext=y_axis, tickvals=list(range(-len(df), len(df) + 1)))
st.subheader('Fragment Ion Plot')
st.plotly_chart(fig, use_container_width=True)
