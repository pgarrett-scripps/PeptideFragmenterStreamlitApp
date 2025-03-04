from typing import List, Tuple, Optional, Dict

import pandas as pd
import streamlit as st
import peptacular as pt
import streamlit_permalink as stp

st.set_page_config(page_title="peptidefragmenter", page_icon=":bomb:", layout="wide")

# Default values
DEFAULT_PEPTIDE = '[Acetyl]-PEPTIDES[UNIMOD:21]'
DEFAULT_CHARGE = 2
DEFAULT_MASS_TYPE = 'monoisotopic'
DEFAULT_FRAGMENT_TYPES = {'a','b', 'x', 'y'}
DEFAULT_USE_MASS_BOUNDS = False
DEFAULT_MIN_MZ = 150.0
DEFAULT_MAX_MZ = 2000.0
DEFAULT_PRECISION = 5
DEFAULT_ROW_PADDING = 5
DEFAULT_COLUMN_PADDING = 15
DEFAULT_SHOW_BORDERS = False
DEFAULT_A_COLOR = '#8c564b'
DEFAULT_B_COLOR = '#1f77b4'
DEFAULT_C_COLOR = '#2ca02c'
DEFAULT_X_COLOR = '#ff7f0e'
DEFAULT_Y_COLOR = '#d62728'
DEFAULT_Z_COLOR = '#9467bd'

# if no query params update to default
if not st.query_params:
    st.query_params.clear()
    st.query_params['peptide'] = DEFAULT_PEPTIDE
    st.query_params['charge'] = DEFAULT_CHARGE
    st.query_params['mass_type'] = DEFAULT_MASS_TYPE
    for ft in DEFAULT_FRAGMENT_TYPES:
        st.query_params[ft] = True
    st.query_params['mass_bounds'] = DEFAULT_USE_MASS_BOUNDS
    st.query_params['min_mass'] = DEFAULT_MIN_MZ
    st.query_params['max_mass'] = DEFAULT_MAX_MZ
    st.query_params['decimal_places'] = DEFAULT_PRECISION
    st.query_params['row_padding'] = DEFAULT_ROW_PADDING
    st.query_params['column_padding'] = DEFAULT_COLUMN_PADDING
    st.query_params['show_borders'] = DEFAULT_SHOW_BORDERS
    st.rerun()

# Sidebar: Peptide Fragmenter input
with st.sidebar:
    st.title('Peptide Fragmenter :bomb:')
    st.caption("""
    A peptide fragment ion calculator ([ProForma 2.0 compliant](https://github.com/HUPO-PSI/ProForma/blob/master/SpecDocument/ProForma_v2_draft15_February2022.pdf)).""")

    st.caption('''**This pages URL automatically updates with your input, and can be shared with others.**''')

    peptide_help_msg = """
    **Peptide Sequence**: Enter the peptide sequence to fragment. Include modifications in square brackets.
    """

    st.subheader('Sequence', divider='grey')

    peptide_sequence = stp.text_input('Peptide',
                   value=DEFAULT_PEPTIDE,
                   max_chars=2000,
                   help=peptide_help_msg,
                   url_key='peptide')
    st.caption(
        '''Common mods: C[Carbamidomethyl], M[Oxidation], [Acetyl]-''')



    c1, c2 = st.columns(2)

    with c1:
        charge = stp.number_input('Charge',
                                 min_value=0,
                                 value=DEFAULT_CHARGE,
                                 help='Charge state of the peptide',
                                  url_key='charge')
    with c2:
        charge_aduct = stp.text_input('Aduct',
                                    value='+H',
                                    help='Adduct to use for charge state calculation',
                                    url_key='adduct', disabled=True)

    mass_type = stp.radio(label='Mass Type',
                         options=['monoisotopic', 'average'],
                         help='Mass type to use for fragment calculation',
                         index=['monoisotopic', 'average'].index(DEFAULT_MASS_TYPE),
                         horizontal=True,
                            url_key='mass_type')

    is_monoisotopic = mass_type == 'monoisotopic'


    try:
        annotation = pt.parse(peptide_sequence)
    except Exception as e:
        st.error(f'Error parsing peptide sequence: {e}')
        st.stop()

    # if contains chage state error
    if annotation.charge is not None:
        st.error('Peptide sequence cannot contain charge state!')
        st.stop()

    # Check peptide AA count is within limits
    if len(annotation) > 1000:
        st.error(f'Peptide length cannot exceed {1000} amino acids')
        st.stop()

    # if contains adduct error
    if annotation.charge_adducts is not None:
        st.error('Peptide sequence cannot contain adduct!')
        st.stop()

    if annotation.contains_sequence_ambiguity() or annotation.contains_residue_ambiguity() or annotation.contains_mass_ambiguity():
        st.error('Sequence cannot contain ambiguity!')
        st.stop()

    if len(annotation) == 0:
        st.error('Peptide sequence cannot be empty')
        st.stop()

    st.subheader('Fragment Ions', divider='grey')
    c1, c2, c3 = st.columns(3)

    with c1:
        a = stp.checkbox('a', value='a' in DEFAULT_FRAGMENT_TYPES, url_key='a')
        x = stp.checkbox('x', value='x' in DEFAULT_FRAGMENT_TYPES, url_key='x')

    with c2:
        b = stp.checkbox('b', value='b' in DEFAULT_FRAGMENT_TYPES, url_key='b')
        y = stp.checkbox('y', value='y' in DEFAULT_FRAGMENT_TYPES, url_key='y')

    with c3:
        c = stp.checkbox('c', value='c' in DEFAULT_FRAGMENT_TYPES, url_key='c')
        z = stp.checkbox('z', value='z' in DEFAULT_FRAGMENT_TYPES, url_key='z')

    fragment_types = [ft for ft, flag in zip('abcxyz', [a, b, c, x, y, z]) if flag]

    st.subheader('Additional Options', divider='grey')

    use_mass_bounds = stp.checkbox('Use Mass Bounds', value=DEFAULT_USE_MASS_BOUNDS, url_key='mass_bounds')

    min_mz, max_mz = None, None
    if use_mass_bounds:
        c1, c2 = st.columns(2)
        with c1:
            min_mz = stp.number_input('Min m/z', value=DEFAULT_MIN_MZ, url_key='min_mz', step=100.0)
        with c2:
            max_mz = stp.number_input('Max m/z', value=DEFAULT_MAX_MZ, url_key='max_mz', step=100.0)

    with st.expander('Format Options'):
        precision = stp.number_input('Decimal Places', value=DEFAULT_PRECISION, min_value=0, max_value=10,
                                         url_key='decimal_places')

        c1, c2 = st.columns(2)
        with c1:

            row_padding = stp.number_input('Row Padding (px)', value=DEFAULT_ROW_PADDING, min_value=0, max_value=100,
                                           url_key='row_padding')
        with c2:
            column_padding = stp.number_input('Column Padding (px)', value=DEFAULT_COLUMN_PADDING, min_value=0, max_value=100,
                                              url_key='column_padding')

        show_borders = stp.checkbox('Show Borders', value=DEFAULT_SHOW_BORDERS, url_key='show_borders')


        st.subheader('Fragment Colors', divider='grey')
        c1, c2, c3 = st.columns(3)
        with c1:
            a_color = stp.color_picker('a color', DEFAULT_A_COLOR)
            x_color = stp.color_picker('x color', DEFAULT_X_COLOR)

        with c2:
            b_color = stp.color_picker('b color', DEFAULT_B_COLOR)
            y_color = stp.color_picker('y color', DEFAULT_Y_COLOR)

        with c3:
            c_color = stp.color_picker('c color', DEFAULT_C_COLOR)
            z_color = stp.color_picker('z color', DEFAULT_Z_COLOR)


def create_fragment_table(sequence: str, ion_types: List[str], charges: List[int], monoisotopic: bool) -> Tuple[
    List, pd.DataFrame]:
    fragments = pt.fragment(sequence=sequence,
                         ion_types=ion_types,
                         charges=charges,
                         monoisotopic=monoisotopic)

    # convert list of dataclasses to list of dicts
    frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])

    frag_df['number'] = None
    # for forward ions (a,b,c) set number to frag.end
    frag_df['number'] = frag_df.apply(lambda row: row['end'] if row['ion_type'] in 'abc' else row['number'], axis=1)

    # for reverse ions (x,y,z) set number to frag.start
    frag_df['number'] = frag_df.apply(lambda row: row['start'] if row['ion_type'] in 'xyz' else row['number'], axis=1)

    return fragments, frag_df


def style_fragment_table(
    sequence: str,
    fragment_types: List[str],
    charge: int,
    is_monoisotopic: bool,
    color_map: Optional[Dict[str, str]] = None,
    show_borders: bool = True,
    aa_col: Optional[str] = "Seq",
    pos_col: Optional[str] = "#>",
    neg_col: Optional[str] = "<#",
    caption: Optional[str] = None,
    decimal_places: int = 4,
    row_padding: int = 4,
    column_padding: int = 10,
        min_mass: Optional[float] = None,
        max_mass: Optional[float] = None,
):
    # Define default colors
    default_colors = {
        "A": "#8B4513",  # Brown
        "B": "#1f77b4",  # Blue
        "C": "#2ca02c",  # Green
        "X": "#ff7f0e",  # Orange
        "Y": "#d62728",  # Red
        "Z": "#9467bd",  # Purple
    }

    # If the user provides a color map, update defaults
    if color_map:
        color_map = {k.upper(): v for k, v in color_map.items()}
        default_colors.update(color_map)

    # Generate fragment data
    fragments, frag_df = create_fragment_table(sequence, fragment_types, [charge], is_monoisotopic)

    if frag_df.empty:
        st.warning("No fragments found. Please check your input and try again.")
        st.stop()

    # Drop unnecessary columns
    frag_df = frag_df.drop(columns=['isotope', 'loss', 'parent_sequence'])
    frag_df.sort_values(by=['charge', 'ion_type', 'start'], inplace=True)
    frag_df.drop_duplicates(subset=['charge', 'ion_type', 'start', 'end'], inplace=True)

    components = pt.split(sequence)
    data = {aa_col: components} if aa_col else {}

    # Process fragment ions
    for ion_type in sorted(fragment_types):

        if ion_type not in 'abcxyz':
            continue

        ion_df = frag_df[(frag_df['ion_type'] == ion_type) & (frag_df['charge'] == charge)]
        ion_df.sort_values(by=['number'], inplace=True)
        frags = ion_df['mz'].tolist()
        if ion_type in 'XYZ':
            frags = frags[::-1]
        data[ion_type.upper()] = frags

    # Create DataFrame
    df = pd.DataFrame(data)

    forward_cols = [col for col in df.columns if 'A' == col or 'B' == col or 'C' == col]
    reverse_cols = [col for col in df.columns if 'X' == col or 'Y' == col or 'Z' == col]

    df = df[forward_cols + [aa_col] + reverse_cols]

    # Add positional columns
    if pos_col and forward_cols:
        df.insert(0, pos_col, list(range(1, len(df) + 1)))
    if neg_col and reverse_cols:
        df[neg_col] = list(range(len(df), 0, -1))
    border = "1px solid #999" if show_borders else "none"

    styles = [
        {'selector': 'table', 'props': [
            ('border-collapse', 'collapse'),
            ('border-spacing', '0'),
            ('border', border),
        ]},
        {'selector': 'th, td', 'props': [
            ('text-align', 'center'),
            ('padding', f'{row_padding}px {column_padding}px'),
            ('line-height', '1'),
            ('border', border),
        ]},
        {'selector': 'tr', 'props': [
            ('border', border),
        ]},
        {'selector': 'th', 'props': [
            ('background-color', '#ffffcc'),
            ('font-weight', 'bold'),
        ]},
    ]

    styled_df = df.style \
        .hide(axis='index') \
        .set_table_styles(styles)

    def highlight_columns(val, color):
        return f'color: {color}; font-weight: bold;' if val else ''

    if 'A' in df.columns:
        styled_df.applymap(lambda val: highlight_columns(val, default_colors['A']), subset=['A'])

    if 'B' in df.columns:
        styled_df.applymap(lambda val: highlight_columns(val, default_colors['B']), subset=['B'])

    if 'C' in df.columns:
        styled_df.applymap(lambda val: highlight_columns(val, default_colors['C']), subset=['C'])

    if 'X' in df.columns:
        styled_df.applymap(lambda val: highlight_columns(val, default_colors['X']), subset=['X'])

    if 'Y' in df.columns:
        styled_df.applymap(lambda val: highlight_columns(val, default_colors['Y']), subset=['Y'])

    if 'Z' in df.columns:
        styled_df.applymap(lambda val: highlight_columns(val, default_colors['Z']), subset=['Z'])

    # color background of any cell to #ffcccc with a value greater outside of mass bounds
    if min_mass and max_mass:
        styled_df.applymap(lambda val: 'background-color: #ffcccc' if val > max_mass or val < min_mass else '', subset=forward_cols + reverse_cols)

    if caption:
        styled_df.set_caption(f'{caption}')

    # Get indices of the max values in 'C' and 'X' columns
    max_index_C = df['C'].idxmax() if 'C' in df.columns else None
    max_index_X = df['X'].idxmax() if 'X' in df.columns else None

    # Define function to replace max values visually without changing the DataFrame
    def hide_max_value(val, col_name):
        if (col_name == 'C' and df.loc[max_index_C, 'C'] == val) or (
                col_name == 'X' and df.loc[max_index_X, 'X'] == val):
            return 'color: transparent;'  # Hide the value (alternative: 'color: white;' or 'visibility: hidden;')
        return ''

    # Apply styling to hide max values
    if 'C' in df.columns:
        styled_df.applymap(lambda val: hide_max_value(val, 'C'), subset=['C'])
    if 'X' in df.columns:
        styled_df.applymap(lambda val: hide_max_value(val, 'X'), subset=['X'])

    styled_df = styled_df.format(precision=decimal_places)

    return styled_df


frag_colors = {
    'a': a_color,
    'b': b_color,
    'c': c_color,
    'x': x_color,
    'y': y_color,
    'z': z_color
}


mass_type_abr = 'Monoisotopic' if is_monoisotopic else 'Avgerage'

style_df = style_fragment_table(
    sequence=annotation.serialize(),
    fragment_types=fragment_types,
    charge=charge,
    is_monoisotopic=is_monoisotopic,
    show_borders=show_borders,  # No outer borders
    decimal_places=precision,
    row_padding=row_padding,
    column_padding=column_padding,
    min_mass=min_mz,
    max_mass=max_mz,
    color_map=frag_colors,
    caption=f'<b>{annotation.serialize()}</b><br>Mass Type: {mass_type_abr}',
)

def center_table(val):
    html = style_df.to_html()
    # Update the column headers to include the superscript charge state
    for col in ["A", "B", "C", "X", "Y", "Z"]:
        html = html.replace(f'{col}</th>', f'{col}<sup>+{charge}</sup></th>')

    st.markdown(html, unsafe_allow_html=True)


st.subheader('Peptide Fragmenter Results')
#st.markdown(f'**Sequence:** {annotation.serialize()}')
#st.markdown(f'**Stripped Sequence:** {annotation.sequence}')

try:
    neutral_sequence_mass = pt.mass(annotation, monoisotopic=is_monoisotopic, ion_type='p', charge=0)
except Exception as e:
    st.error(f'Error calculating peptide mass: {e}')
    st.stop()

try:
    sequence_mz = pt.mz(annotation, monoisotopic=is_monoisotopic, ion_type='p', charge=charge)
except Exception as e:
    st.error(f'Error calculating peptide mass: {e}')
    st.stop()

center_table(style_df)

if use_mass_bounds:
    st.markdown(f'**Bounds:** {min_mz} - {max_mz} *m/z*')


st.caption('Made with [peptacular](https://pypi.org/project/peptacular/)')

