import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers
from rdkit.Chem.rdmolops import CombineMols
from rdkit.Chem.rdchem import BondType

class Monomer:
    def __init__(self, id, generation, smiles, attachment_map_numbers, attachment_points_count=None):
        """
        Initializes a Monomer with its chemical structure and attachment points.

        Args:
            id (int): Unique identifier for the monomer.
            generation (int): The generation number.
            smiles (str): SMILES string for the monomer, including dummy atoms for attachment points.
                          Example: "c1([*:1])cc([*:2])cc([*:3])c1" for a core with 3 points.
            attachment_map_numbers (list): A list of atom map numbers corresponding to
                                           the dummy atoms that serve as attachment points.
                                           e.g., [1, 2] for "[*:1]" and "[*:2]".
            attachment_points_count (int, optional): The number of attachment points for *children*.
                                                      If None, derived from len(attachment_map_numbers)
                                                      after accounting for the parent connection.
        """
        self.id = id
        self.generation = generation
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            raise ValueError(f"Invalid SMILES string for monomer {id}: {smiles}")

        self.mol = Chem.AddHs(self.mol) # Add hydrogens before embedding and analysis
        # Generate initial 3D for monomer, without optimization yet for individual units
        AllChem.EmbedMolecule(self.mol, AllChem.ETKDGv3())

        self.attachment_map_numbers = attachment_map_numbers
        self.connections = {}  # Dictionary to store connected Monomer objects (child monomers)
                               # keyed by the *parent's* attachment_map_number.
        self.parent_connection_atom_idx = None # The atom index in *this* monomer that connects to its parent.
        self.global_atom_indices = {} # Maps original atom index in self.mol to global index in combined_mol

        self.attachment_atom_indices = self._get_attachment_atom_indices_revised()
        if len(self.attachment_atom_indices) != len(attachment_map_numbers):
             raise ValueError(f"Monomer {id} SMILES has {len(self.attachment_atom_indices)} valid "
                              f"attachment point atoms but {len(attachment_map_numbers)} map numbers were specified.")

        # This `attachment_points` refers to the number of *child* connections this monomer can make
        # For a branching monomer, this would be the total minus one for the parent connection
        # We will handle this more robustly in CustomMonomer
        self.attachment_points = attachment_points_count if attachment_points_count is not None \
                                 else len(self.attachment_map_numbers)

    def _get_attachment_atom_indices_revised(self):
        """Finds the real atom indices that dummy atoms are connected to."""
        indices = {}
        for atom in self.mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_num = atom.GetIntProp("molAtomMapNumber")
                if map_num in self.attachment_map_numbers:
                    # Find the atom directly bonded to this dummy atom
                    # This atom is the actual point of connection
                    for neighbor in atom.GetNeighbors():
                        indices[map_num] = neighbor.GetIdx()
                        break # Assume only one neighbor for dummy atom
        # Return indices in the same order as attachment_map_numbers
        return [indices[map_num] for map_num in self.attachment_map_numbers if map_num in indices]

class CustomMonomer(Monomer):
    def __init__(self, id, generation, smiles, parent_attachment_map_num, child_attachment_map_nums):
        """
        Initializes a Monomer specifically for dendrimer branching.

        Args:
            id (int): Unique identifier.
            generation (int): Generation number.
            smiles (str): SMILES string with dummy atoms.
            parent_attachment_map_num (int or None): Atom map number for connection to parent. None for core.
            child_attachment_map_nums (list): List of atom map numbers for connections to children.
        """
        all_map_numbers = child_attachment_map_nums
        if parent_attachment_map_num is not None:
            all_map_numbers = [parent_attachment_map_num] + child_attachment_map_nums

        super().__init__(id, generation, smiles, all_map_numbers,
                         attachment_points_count=len(child_attachment_map_nums))

        self.parent_attachment_map_num = parent_attachment_map_num
        self.child_attachment_map_nums = child_attachment_map_nums

        if self.generation == 0: # Core has no parent connection atom
             self.parent_connection_atom_idx = None
        else:
            # Find the atom index for the parent_attachment_map_num in *this* monomer
            parent_conn_found = False
            for map_num, atom_idx in zip(self.attachment_map_numbers, self.attachment_atom_indices):
                if map_num == self.parent_attachment_map_num:
                    self.parent_connection_atom_idx = atom_idx
                    parent_conn_found = True
                    break
            if not parent_conn_found:
                raise ValueError(f"Monomer {id} (Gen {generation}): Parent attachment map number "
                                 f"{parent_attachment_map_num} not found in SMILES or its connections.")

        # `self.attachment_points` (from base Monomer) now correctly represents number of child connections
        self.attachment_points = len(self.child_attachment_map_nums)

        # Store a more direct mapping for child attachment points to their atom indices
        self.child_attachment_atom_indices_map = {}
        for map_num, atom_idx in zip(self.attachment_map_numbers, self.attachment_atom_indices):
            if map_num in self.child_attachment_map_nums:
                self.child_attachment_atom_indices_map[map_num] = atom_idx


def generate_dendrimer_v2(core_data, branch_data, generations):
    """
    Generates a dendrimer structure with CustomMonomer objects.

    Args:
        core_data (dict): {'smiles': str, 'child_attachment_map_nums': list}
        branch_data (dict): {'smiles': str, 'parent_attachment_map_num': int, 'child_attachment_map_nums': list}
        generations (int): The number of generations to generate (core + N layers).

    Returns:
        A list of CustomMonomer objects representing the dendrimer.
    """
    monomer_count = 0
    core = CustomMonomer(monomer_count, 0,
                         core_data['smiles'],
                         parent_attachment_map_num=None, # Core has no parent attachment
                         child_attachment_map_nums=core_data['child_attachment_map_nums'])
    dendrimer = [core]
    monomer_count += 1

    for gen in range(1, generations + 1):
        new_monomers = []
        # Iterate over a *copy* of dendrimer to avoid modifying list while iterating
        for parent_monomer in list(dendrimer):
            if parent_monomer.generation == gen - 1:
                # Iterate through the parent's *child* attachment points
                current_child_connections_made = 0
                for parent_child_map_num in parent_monomer.child_attachment_map_nums:
                    if current_child_connections_made < parent_monomer.attachment_points: # Check against allowed children
                        child_monomer = CustomMonomer(monomer_count, gen,
                                                      branch_data['smiles'],
                                                      parent_attachment_map_num=branch_data['parent_attachment_map_num'],
                                                      child_attachment_map_nums=branch_data['child_attachment_map_nums'])

                        # Store connection: parent's child_attachment_map_num -> child monomer object
                        parent_monomer.connections[parent_child_map_num] = child_monomer
                        new_monomers.append(child_monomer)
                        monomer_count += 1
                        current_child_connections_made += 1
        dendrimer.extend(new_monomers)
    return dendrimer

def plot_dendrimer(dendrimer):
    """
    Plots a dendrimer structure using NetworkX and Matplotlib.

    Args:
        dendrimer: A list of Monomer objects representing the dendrimer
                   generated by the generate_dendrimer function.
    """
    graph = nx.Graph()
    for monomer in dendrimer:
        graph.add_node(monomer.id, generation=monomer.generation)
        for connected_monomer in monomer.connections.values():
            graph.add_edge(monomer.id, connected_monomer.id)

    pos = nx.kamada_kawai_layout(graph, scale=1.0) # Use kamada_kawai_layout

    node_colors = [monomer.generation for monomer in dendrimer]

    plt.figure(figsize=(10, 8)) # Adjust figure size for better visibility
    nx.draw(graph, pos, with_labels=False, node_size=100,
             node_color=node_colors, cmap=plt.cm.viridis)
    plt.title("Dendrimer Topological Structure")
    plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.viridis), ax=plt.gca(), label="Generation")
    plt.show()

def get_coordinates(dendrimer, scale=25.0):
    """
    Calculates the 2D coordinates of each monomer in the dendrimer.

    Args:
        dendrimer: A list of Monomer objects representing the dendrimer.
        scale (float): Scaling factor for the layout.

    Returns:
        A dictionary mapping monomer IDs to (x, y) coordinates.
    """
    graph = nx.Graph()
    for monomer in dendrimer:
        graph.add_node(monomer.id, generation=monomer.generation)
        for connected_monomer in monomer.connections.values():
            graph.add_edge(monomer.id, connected_monomer.id)

    pos = nx.kamada_kawai_layout(graph, scale=scale)
    return pos

def build_bonded_dendrimer_molecule(dendrimer, coordinates_2d=None):
    """
    Builds a single RDKit molecule with chemical bonds based on dendrimer connectivity,
    and attempts to arrange it based on 2D coordinates before 3D optimization.

    Args:
        dendrimer (list): List of Monomer objects.
        coordinates_2d (dict, optional): Dictionary mapping monomer IDs to (x, y) coordinates
                                        from get_coordinates, for initial 3D placement.

    Returns:
        Chem.Mol: The combined RDKit molecule with bonds, or None if an error occurs.
    """
    if not dendrimer:
        print("Dendrimer list is empty.")
        return None

    rw_mol = Chem.RWMol()
    # IMPORTANT: Add an empty conformer to rw_mol upfront to set atom positions
    main_conformer = Chem.Conformer()
    rw_mol.AddConformer(main_conformer, assignId=True)

    # Map (monomer_id, original_atom_idx) to new_rw_mol_atom_idx
    global_atom_map = {}

    # 1. Add all atoms from all monomers to the RWMol (excluding dummy atoms)
    # and store their new global indices.
    for monomer in dendrimer:
        # Create a temporary RWMol for the current monomer to safely remove dummy atoms
        # Using Chem.Mol(monomer.mol) creates a copy, which is good.
        temp_rw_monomer_mol = Chem.Mol(monomer.mol) # Start with Mol copy, then convert to RWMol if needed for atom removal

        # This part ensures that the conformer is retained if it exists
        if temp_rw_monomer_mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(temp_rw_monomer_mol, AllChem.ETKDGv3()) # Ensure 3D coords for positions

        # Now get the conformer from this temp molecule
        monomer_conf = temp_rw_monomer_mol.GetConformer()


        dummy_atom_indices_in_temp = sorted([atom.GetIdx() for atom in temp_rw_monomer_mol.GetAtoms() if atom.HasProp("molAtomMapNumber")], reverse=True)

        # To remove atoms from a mol and maintain conformer, it's safer to build a new mol
        # or use RWMol features for removal more carefully.
        # A simpler way to get "cleaned_monomer_mol" without explicit removal:
        # Iterate over atoms of the *original* monomer.mol and only add non-dummy ones.

        # Add atoms and bonds from the monomer molecule to the main rw_mol
        # and store global indices
        monomer_original_to_global_map = {}
        
        # Keep track of which atoms in the original monomer's molecule are NOT dummy atoms,
        # and thus will be added to the rw_mol.
        atoms_to_add_from_monomer = []
        for atom in temp_rw_monomer_mol.GetAtoms():
            if not atom.HasProp("molAtomMapNumber"):
                atoms_to_add_from_monomer.append(atom)

        if not atoms_to_add_from_monomer:
            print(f"Warning: Monomer {monomer.id} has no non-dummy atoms to add. Skipping.")
            continue # Skip this monomer if it's empty after considering dummy atoms

        # Add atoms and their positions to the main rw_mol's conformer
        for atom in atoms_to_add_from_monomer:
            new_idx = rw_mol.AddAtom(Chem.Atom(atom.GetAtomicNum()))
            # Copy properties (e.g., formal charge, isotope, chirality if present)
            rw_mol.GetAtomWithIdx(new_idx).SetFormalCharge(atom.GetFormalCharge())
            rw_mol.GetAtomWithIdx(new_idx).SetIsotope(atom.GetIsotope())
            
            # Set the atom's position in the main_conformer
            main_conformer.SetAtomPosition(new_idx, monomer_conf.GetAtomPosition(atom.GetIdx()))

            global_atom_map[(monomer.id, atom.GetIdx())] = new_idx
            monomer_original_to_global_map[atom.GetIdx()] = new_idx
            monomer.global_atom_indices[atom.GetIdx()] = new_idx

        # Add internal bonds for this monomer (only between non-dummy atoms)
        for bond in temp_rw_monomer_mol.GetBonds():
            bgn_idx_orig = bond.GetBeginAtomIdx()
            end_idx_orig = bond.GetEndAtomIdx()
            
            # Check if both atoms connected by the bond are non-dummy and were added to rw_mol
            if not (temp_rw_monomer_mol.GetAtomWithIdx(bgn_idx_orig).HasProp("molAtomMapNumber") or
                    temp_rw_monomer_mol.GetAtomWithIdx(end_idx_orig).HasProp("molAtomMapNumber")):
                
                # Ensure the atoms actually exist in our global map (i.e., they were added)
                if bgn_idx_orig in monomer_original_to_global_map and end_idx_orig in monomer_original_to_global_map:
                    rw_mol.AddBond(monomer_original_to_global_map[bgn_idx_orig],
                                   monomer_original_to_global_map[end_idx_orig],
                                   bond.GetBondType())
                else:
                    print(f"Warning: Skipping internal bond in monomer {monomer.id} because one or both atoms were not added (might be dummy or related issues).")

    # 2. Add inter-monomer bonds based on connections
    for parent_monomer in dendrimer:
        if not parent_monomer.connections:
            continue

        for parent_child_map_num, child_monomer in parent_monomer.connections.items():
            # Get the global atom index for the parent's attachment point
            # We need the original atom index in the parent's *initial* mol
            try:
                parent_conn_atom_orig_idx = parent_monomer.child_attachment_atom_indices_map[parent_child_map_num]
                parent_conn_atom_global_idx = global_atom_map[(parent_monomer.id, parent_conn_atom_orig_idx)]
            except KeyError:
                print(f"Error: Could not find global atom for parent {parent_monomer.id}'s attachment point {parent_child_map_num}. Skipping bond.")
                continue

            # Get the global atom index for the child's attachment point to its parent
            # child_monomer.parent_connection_atom_idx is already the original atom index in child's initial mol
            try:
                child_conn_atom_orig_idx = child_monomer.parent_connection_atom_idx
                child_conn_atom_global_idx = global_atom_map[(child_monomer.id, child_conn_atom_orig_idx)]
            except KeyError:
                print(f"Error: Could not find global atom for child {child_monomer.id}'s parent connection. Skipping bond.")
                continue

            # Form a single bond (you might need logic for double/aromatic etc. based on specific chemistry)
            if not rw_mol.GetBondBetweenAtoms(parent_conn_atom_global_idx, child_conn_atom_global_idx):
                rw_mol.AddBond(parent_conn_atom_global_idx, child_conn_atom_global_idx, BondType.SINGLE)
            else:
                print(f"Warning: Bond already exists between {parent_monomer.id}-{child_monomer.id}. Skipping.")

    # 3. Initial 3D placement based on 2D coordinates (if provided)
    if coordinates_2d and rw_mol.GetNumAtoms() > 0:
        # Use the already existing main_conformer
        conf = rw_mol.GetConformer(0)

        for monomer in dendrimer:
            if monomer.id in coordinates_2d and monomer.global_atom_indices:
                x_target, y_target = coordinates_2d[monomer.id]

                # Calculate the current geometrical center of this monomer's atoms in the rw_mol
                current_monomer_coords = []
                for orig_idx, global_idx in monomer.global_atom_indices.items():
                    if global_idx < conf.GetNumAtoms(): # Ensure index is valid in current conformer
                        current_monomer_coords.append(conf.GetAtomPosition(global_idx))

                if not current_monomer_coords:
                    continue

                current_centroid = Chem.rdGeometry.Point3D(
                    sum(p.x for p in current_monomer_coords) / len(current_monomer_coords),
                    sum(p.y for p in current_monomer_coords) / len(current_monomer_coords),
                    sum(p.z for p in current_monomer_coords) / len(current_monomer_coords)
                )

                # Calculate translation vector to move this group to the 2D layout target
                translation_vector = Chem.rdGeometry.Point3D(x_target - current_centroid.x,
                                                             y_target - current_centroid.y,
                                                             -current_centroid.z) # Flatten to Z=0

                # Apply translation to all atoms belonging to this monomer
                for orig_idx, global_idx in monomer.global_atom_indices.items():
                    if global_idx < conf.GetNumAtoms():
                        current_pos = conf.GetAtomPosition(global_idx)
                        new_pos = current_pos + translation_vector
                        conf.SetAtomPosition(global_idx, new_pos)
            else:
                if monomer.id not in coordinates_2d:
                    print(f"Warning: No 2D coordinates found for monomer {monomer.id}. Skipping initial placement.")
                if not monomer.global_atom_indices:
                    print(f"Warning: Monomer {monomer.id} has no global atom indices. Skipping initial placement.")


    # 4. Finalize and Sanitize
    mol = rw_mol.GetMol()
    if mol.GetNumAtoms() == 0:
        print("Final molecule is empty after processing monomers.")
        return None

    try:
        Chem.SanitizeMol(mol)
        print("Molecule sanitized successfully.")
    except Exception as e:
        print(f"Sanitization failed: {e}. Attempting to recover...")
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_PROPERTIES | Chem.SanitizeFlags.SANITIZE_SYMMRINGS)
            print("Sanitized with recovery flags.")
        except Exception as recover_e:
            print(f"Recovery sanitization also failed: {recover_e}. Molecule might be invalid.")
            return None

    # 5. Add Hs if not already added and Perform Energy Minimization (Crucial for a realistic 3D structure)
    try:
        mol = Chem.AddHs(mol)
        if mol.GetNumConformers() == 0:
             print("Embedding molecule for optimization (if no conformer exists after sanitization)...")
             AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        if AllChem.MMFFHasMolecule(mol):
            AllChem.MMFFOptimizeMolecule(mol)
            print("Molecule optimized using MMFF.")
        elif AllChem.UFFHasMolecule(mol): # Use UFF if MMFF is not applicable
            AllChem.UFFOptimizeMolecule(mol)
            print("Molecule optimized using UFF.")
        else:
            print("No suitable force field found for optimization.")

    except Exception as e:
        print(f"Energy optimization failed: {e}. The 3D structure might be unoptimized.")

    return mol

# --- Example Usage ---

# Monomer and CustomMonomer classes as defined previously

# Core with 3 child attachment points (e.g., a 1,3,5-substituted benzene)
core_data = {
    'smiles': "c1([*:1])cc([*:2])cc([*:3])c1",
    'child_attachment_map_nums': [1, 2, 3]
}

# Branching Monomer with 1 parent attachment point and 2 child attachment points
branch_data = {
    'smiles': "c1([*:1])cc([*:2])cc([*:3])c1",
    'parent_attachment_map_num': 1,
    'child_attachment_map_nums': [2, 3]
}

# Generate a dendrimer with 2 generations (core + 2 layers of branches)
dendrimer_monomers = generate_dendrimer_v2(core_data, branch_data, 5)
print(f"Generated {len(dendrimer_monomers)} monomers for the dendrimer.")

# Get 2D coordinates for spatial arrangement using NetworkX
dendrimer_coordinates = get_coordinates(dendrimer_monomers, scale=50.0) # Increased scale for larger structures

# Build the chemically bonded RDKit molecule
print("\nAttempting to build the bonded dendrimer molecule...")
final_dendrimer_mol = build_bonded_dendrimer_molecule(dendrimer_monomers, dendrimer_coordinates)

if final_dendrimer_mol:
    print(f"\nFinal bonded dendrimer molecule has {final_dendrimer_mol.GetNumAtoms()} atoms and {final_dendrimer_mol.GetNumBonds()} bonds.")
    print("\n--- Final Bonded Dendrimer Molecule (XYZ Block) ---")
    print(Chem.MolToXYZBlock(final_dendrimer_mol))

    try:
        writer = Chem.SDWriter('dendrimer.sdf')
        writer.write(final_dendrimer_mol)
        writer.close()
        print("\nSaved dendrimer to dendrimer.sdf (open with PyMOL, Avogadro, etc. to visualize 3D structure).")
    except Exception as e:
        print(f"Failed to save SDF file: {e}")

    plot_dendrimer(dendrimer_monomers)
else:
    print("Failed to build the bonded dendrimer molecule.")
