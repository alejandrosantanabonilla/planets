import numpy as np

def convert_to_atom_list(data):
    """
    Converts a NumPy array of the format [['atom', x, y, z], ...] into a list of lists
    where each inner list is [atom_name, (x, y, z)].

    Args:
        data (numpy.ndarray): The input NumPy array.

    Returns:
        list: A list of lists in the desired atom format.
    """
    atom = []
    for row in data:
        atom_name = row[0]
        coords = tuple(map(float, row[1:]))  # Convert coordinates to floats and create a tuple
        atom.append([atom_name, coords])
    return atom

def generate_and_split_characters(combinations, shuffle=False):
    """
    Generates repeated strings based on given combinations, splits them into characters, and optionally shuffles the result.

    Args:
        combinations: A list of tuples, where each tuple contains (character, repetitions).
        shuffle (bool): Whether to shuffle the characters (default: False).

    Returns:
        A list of characters representing the split and potentially shuffled strings.
    """

    split_elements = []
    for char, reps in combinations:
        split_elements.extend([char] * reps)  # Create the list with repeated elements
        
    if shuffle:
        random.shuffle(split_elements)  # Shuffle in-place if needed

    return split_elements
