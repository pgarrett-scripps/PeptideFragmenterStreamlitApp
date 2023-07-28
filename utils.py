
COLOR_DICT = {'a': 'brown', 'b': 'blue', 'c': 'green', 'x': 'orange', 'y': 'red', 'z': 'purple', 'i': 'magenta'}


# Function to apply color according to ion type
def color_by_ion_type(col):
    ion_type = col.name[-1]
    color = COLOR_DICT.get(ion_type, 'grey')  # get color or default to grey if not found
    return ['color: %s' % color] * len(col)

def get_fragment_color(row):
    color = COLOR_DICT.get(row['ion_type'], 'black')

    if row['internal'] is True:
        color = 'magenta'

    return color

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False