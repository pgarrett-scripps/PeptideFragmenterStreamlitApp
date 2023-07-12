
COLOR_DICT = {'a': 'brown', 'b': 'blue', 'c': 'green', 'x': 'orange', 'y': 'red', 'z': 'purple'}


# Function to apply color according to ion type
def color_by_ion_type(col):
    ion_type = col.name[0]
    color = COLOR_DICT.get(ion_type, 'grey')  # get color or default to grey if not found
    return ['color: %s' % color] * len(col)