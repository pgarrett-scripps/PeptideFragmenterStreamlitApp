COLOR_DICT = {
    'a': '#8c564b',  # Darker Brown
    'b': '#1f77b4',  # Vivid Blue
    'c': '#2ca02c',  # Green
    'x': '#ff7f0e',  # Bright Orange
    'y': '#d62728',  # Red
    'z': '#9467bd',  # Purple
    'ax': '#e377c2', # Magenta
    'ay': '#f7b6d2', # Pink
    'az': '#c5b0d5', # Lavender
    'bx': '#1f77b4', # Darker Navy (Consider adjusting if too similar to 'b')
    'by': '#aec7e8', # Light Blue
    'bz': '#98df8a', # Light Green
    'cx': '#238b45', # Dark Green
    'cy': '#d62728', # Different shade of Red (Consider changing to avoid duplication)
    'cz': '#d62728'  # (Adjust to introduce distinction, e.g., teal)
}


# Function to apply color according to ion type
def color_by_ion_type(col):
    ion_type = col.name[-1]
    color = COLOR_DICT.get(ion_type, 'grey')  # get color or default to grey if not found
    return ['color: %s' % color] * len(col)


def get_fragment_color(row):
    color = COLOR_DICT.get(row['ion_type'], 'black')

    return color


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
