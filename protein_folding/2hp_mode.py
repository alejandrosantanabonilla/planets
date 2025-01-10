import random
import math
import matplotlib.pyplot as plt

def valid_move(struct: list) -> bool:
    """
    Checks if a protein structure is valid:
    - It is a self-avoiding walk (SAW).
    - Distances between consecutive elements are 1.

    Args:
        struct: The protein structure as a list of coordinates.

    Returns:
        True if the structure is valid, False otherwise.
    """

    # Check for self-avoiding walk (no repeated coordinates)
    if len(struct) != len(set(tuple(coord) for coord in struct)):
        return False

    # Check for unit distances between consecutive elements
    for i in range(1, len(struct)):
        if math.dist(struct[i - 1], struct[i]) != 1:
            return False

    return True

def tail_fold(struct: list) -> list:
    """
    Folds the tail of the protein by applying a random transformation
    (rotation, reflection, or diagonal move). Ensures that a change is always made.

    Args:
        struct: The protein structure as a list of coordinates.

    Returns:
        The updated structure with the tail folded.
    """

    if len(struct) < 3:
        return struct  # No need to fold if the structure is too short

    tail_start = random.randint(1, len(struct) - 2)  # Choose a random tail start point
    tail = struct[tail_start:]
    head = struct[:tail_start]
    original_tail = tail.copy()  # Store a copy of the original tail

    while tail == original_tail:  # Keep trying transformations until the tail changes
        # Choose a random transformation (rotation, reflection, or diagonal move)
        transformation = random.choice(["rotate90", "rotate180", "rotate270",
                                        "reflect_x", "reflect_y", "diagonal"])

        # Apply the chosen transformation to the tail
        if transformation.startswith("rotate"):
            # Calculate the pivot point for rotation (last monomer of the head)
            pivot = head[-1]

            # Apply the rotation to each monomer in the tail
            for i, monomer in enumerate(tail):
                # Translate the monomer so that the pivot is at the origin
                translated_monomer = [monomer[0] - pivot[0], monomer[1] - pivot[1]]

                # Apply the rotation
                if transformation == "rotate90":
                    rotated_monomer = [-translated_monomer[1], translated_monomer[0]]
                elif transformation == "rotate180":
                    rotated_monomer = [-translated_monomer[0], -translated_monomer[1]]
                elif transformation == "rotate270":
                    rotated_monomer = [translated_monomer[1], -translated_monomer[0]]

                # Translate the monomer back to its original position
                tail[i] = [rotated_monomer[0] + pivot[0], rotated_monomer[1] + pivot[1]]

        elif transformation == "reflect_x":
            # Reflect each monomer in the tail across the x-axis
            for i, monomer in enumerate(tail):
                tail[i] = [monomer[0], -monomer[1]]

        elif transformation == "reflect_y":
            # Reflect each monomer in the tail across the y-axis
            for i, monomer in enumerate(tail):
                tail[i] = [-monomer[0], monomer[1]]

        elif transformation == "diagonal":
            # Apply diagonal move to a random monomer in the tail
            try:
                diagonal_move(tail, random.randrange(1, len(tail) - 1))
            except ValueError:
                pass  # Ignore if no valid diagonal move is possible

    folded_struct = head + tail

    # If the move is not valid, try another tail_fold
    if not valid_move(folded_struct):
        return tail_fold(struct)  # Recursive call to try another move

    return folded_struct

def diagonal_move(struct: list, pos: int) -> list:
    """
    Moves the monomer at the given position along a diagonal,
    considering the previous and following monomers.

    Args:
        struct: The protein structure as a list of coordinates.
        pos: The position of the monomer to move.

    Returns:
        The updated structure with the monomer moved diagonally.
    """

    if not 0 < pos < len(struct) - 1:
        raise ValueError("Invalid position for diagonal move.")

    prev_monomer = struct[pos - 1]
    next_monomer = struct[pos + 1]
    curr_monomer = struct[pos]

    possible_moves = []

    # Check for possible diagonal moves based on neighbors
    if prev_monomer[0] == curr_monomer[0]:  # Previous monomer is vertical
        possible_moves.extend([
            [curr_monomer[0] + 1, curr_monomer[1] + 1],
            [curr_monomer[0] - 1, curr_monomer[1] + 1]
        ])
    if prev_monomer[1] == curr_monomer[1]:  # Previous monomer is horizontal
        possible_moves.extend([
            [curr_monomer[0] + 1, curr_monomer[1] + 1],
            [curr_monomer[0] + 1, curr_monomer[1] - 1]
        ])
    if next_monomer[0] == curr_monomer[0]:  # Next monomer is vertical
        possible_moves.extend([
            [curr_monomer[0] + 1, curr_monomer[1] - 1],
            [curr_monomer[0] - 1, curr_monomer[1] - 1]
        ])
    if next_monomer[1] == curr_monomer[1]:  # Next monomer is horizontal
        possible_moves.extend([
            [curr_monomer[0] - 1, curr_monomer[1] + 1],
            [curr_monomer[0] - 1, curr_monomer[1] - 1]
        ])

    # Remove duplicates and invalid moves (collisions)
    possible_moves = [move for move in possible_moves
                      if move not in struct and move != prev_monomer and move != next_monomer]

    if possible_moves:
        struct[pos] = random.choice(possible_moves)

    return struct

def is_valid_sequence(seq: str) -> bool:
    """
    Check if the protein sequence contains only H and/or P and its length is at least 3.

    Args:
        seq (str): The protein sequence (case-insensitive).

    Returns:
        bool: True if the sequence is valid, False otherwise.
    """

    if len(seq) < 3:
        print("The sequence is too short. It must be at least 3.")
        return False

    seq = seq.upper()  # Make the sequence uppercase for case-insensitivity
    unique_chars = set(seq)  # Collect unique characters

    # Check if the sequence contains only 'H' and/or 'P'
    if unique_chars.issubset({'H', 'P'}):
        return True
    else:
        return False

def plot_protein_structure(seq: str, struct: list) -> None:
    """
    Plots the protein structure on a 2D grid with lines connecting 
    nearest neighbors and H/P labels on the monomers.

    Args:
        seq (str): The protein sequence (string of 'H' and 'P').
        struct (list): The protein structure as a list of coordinates.
    """

    x_coords = [coord[0] for coord in struct]
    y_coords = [coord[1] for coord in struct]

    # Create the plot
    fig, ax = plt.subplots()

    # Plot lines between consecutive monomers (nearest neighbors)
    for i in range(1, len(struct)):
        ax.plot([x_coords[i-1], x_coords[i]], [y_coords[i-1], y_coords[i]], 'k-', linewidth=2)

    # Plot H/P labels on the monomers
    for i, (x, y) in enumerate(zip(x_coords, y_coords)):
        ax.text(x, y, seq[i], ha='center', va='center', fontsize=14)

    # Remove axis ticks and labels
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # Set title
    ax.set_title("Protein Structure")

    plt.show()

def linear_struct(seq: str, orientation: str = "horizontal") -> list:
    """
    Create a linear structure of the length of the input sequence.

    Args:
        seq (str): The protein sequence.
        orientation (str, optional): The orientation of the linear structure.
                                        Can be "horizontal" (default) or "vertical".

    Returns:
        list: The linear structure as a list of coordinates.
    """

    struct = []
    if orientation == "horizontal":
        for i in range(len(seq)):
            struct.append([i, 0])
    elif orientation == "vertical":
        for i in range(len(seq)):
            struct.append([0, i])
    else:
        raise ValueError("Invalid orientation. Choose 'horizontal' or 'vertical'.")

    return struct


def count_hh_non_bonded_neighbors(seq: str, struct: list) -> int:
    """
    Counts the number of non-bonded H-H neighbors (contacts) in a protein structure.

    Args:
        seq: The protein sequence (string of 'H' and 'P').
        struct: The protein structure as a list of coordinates.

    Returns:
        The number of non-bonded H-H neighbors.
    """
    count = 0
    for i in range(len(struct)):
        for j in range(i + 2, len(struct)):  # Avoid consecutive neighbors
            if seq[i] == 'H' and seq[j] == 'H' and math.dist(struct[i], struct[j]) == 1:
                count += 1
    return count

def calculate_energy(e: float, seq: str, struct: list) -> float:
    """
    Calculates the energy of the protein structure based on the number of 
    H-H non-bonded neighbors.

    Args:
        e: The energy contribution per H-H bond.
        seq: The protein sequence (string of 'H' and 'P').
        struct: The protein structure as a list of coordinates.

    Returns:
        The calculated energy of the structure.
    """
    all_bonds = count_hh_non_bonded_neighbors(seq, struct)
    energy = -e * all_bonds
    return energy

# Example usage
seq = "HPHPHPHPHPHPPHPHPPHPHPPHHP"  # Example protein sequence
e = 1.0  # Energy per H-H bond

if is_valid_sequence(seq):
    # Create a horizontal linear structure
    structure = linear_struct(seq, orientation="horizontal") 

    # Sequentially fold the structure 10 times
    for _ in range(10000):
        structure = tail_fold(structure.copy())  
        print("Structure:", structure)  # Print the structure after each fold

    # Count and print the H-H non-bonded neighbors
    hh_non_bonded_count = count_hh_non_bonded_neighbors(seq, structure)
    print("H-H Non-bonded neighbors:", hh_non_bonded_count)

    # Calculate and print the energy
    energy = calculate_energy(e, seq, structure)
    print("Energy:", energy)

    # Plot the structure after each fold
    plot_protein_structure(seq, structure)  
else:
    print("Invalid sequence provided.")
