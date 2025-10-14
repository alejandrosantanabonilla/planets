import networkx as nx
import matplotlib.pyplot as plt
from typing import List, Dict, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # Suppress RDKit warnings
from openbabel import pybel as pb

from pysoftk.tools.utils_ob import ff_ob_relaxation
from pysoftk.format_printers.format_mol import Fmt

class Monomer:
    """Base class for a monomer unit."""
    def __init__(self, id: int, generation: int, smiles: str, attachment_map_numbers: List[int]):
        self.id = id
        self.generation = generation
        self.mol = Chem.MolFromSmiles(smiles)
        if not self.mol:
            raise ValueError(f"Invalid SMILES string for monomer {id}: {smiles}")

        self.mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(self.mol, AllChem.ETKDGv3())

        self.attachment_map_numbers = attachment_map_numbers
        self.connections: Dict[int, 'Monomer'] = {}
        self.parent_connection_atom_idx: Optional[int] = None
        self.global_atom_indices: Dict[int, int] = {}
        self.attachment_atom_indices = self._get_attachment_atom_indices()

    def _get_attachment_atom_indices(self) -> List[int]:
        """Finds the real atom indices that the dummy atoms are connected to."""
        indices = {}
        for atom in self.mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                map_num = atom.GetIntProp("molAtomMapNumber")
                if map_num in self.attachment_map_numbers:
                    indices[map_num] = atom.GetNeighbors()[0].GetIdx()
        return [indices[mn] for mn in self.attachment_map_numbers if mn in indices]

class CustomMonomer(Monomer):
    """A specialized Monomer for dendrimer branching units."""
    def __init__(self, id: int, generation: int, smiles: str, parent_attachment_map_num: Optional[int], child_attachment_map_nums: List[int]):
        all_maps = child_attachment_map_nums + ([parent_attachment_map_num] if parent_attachment_map_num else [])
        super().__init__(id, generation, smiles, all_maps)
        self.child_attachment_map_nums = child_attachment_map_nums
        self.child_attachment_atom_indices_map: Dict[int, int] = {}
        
        map_to_idx = dict(zip(self.attachment_map_numbers, self.attachment_atom_indices))
        if parent_attachment_map_num:
            self.parent_connection_atom_idx = map_to_idx.get(parent_attachment_map_num)
        for map_num in self.child_attachment_map_nums:
            self.child_attachment_atom_indices_map[map_num] = map_to_idx.get(map_num)

class SurfaceMonomer(Monomer):
    """A specialized Monomer for terminal surface groups."""
    def __init__(self, id: int, generation: int, smiles: str, parent_attachment_map_num: int):
        super().__init__(id, generation, smiles, [parent_attachment_map_num])
        self.parent_connection_atom_idx = self.attachment_atom_indices[0]

class DendrimerBuilder:
    """A class to construct, optimize, and manage dendrimer molecules."""

    def __init__(self, core_data: dict, branch_data: dict, generations: int, surface_data: Optional[dict] = None):
        self.core_data = core_data
        self.branch_data = branch_data
        self.generations = generations
        self.surface_data = surface_data
        
        self.monomer_graph: List[Monomer] = []
        self.optimized_mol: Optional[pb.Molecule] = None

    def build(self, force_field: str = "UFF", steps: int = 1000, ff_thr: float = 0.1):
        """Orchestrates the dendrimer construction with user-defined optimization parameters."""
        print("Step 1: Generating abstract dendrimer graph...")
        self._generate_graph()
        
        print("Step 2: Building and optimizing 3D molecule...")
        self.optimized_mol = self._build_and_optimize_mol(force_field, steps, ff_thr)
        
        if self.optimized_mol:
            print("Dendrimer construction complete.")
        else:
            print("Error: Dendrimer construction failed.")

    def get_optimized_mol(self) -> Optional[pb.Molecule]:
        """Returns the final, optimized Pybel molecule."""
        return self.optimized_mol

    def plot_topology(self, output_filename: str = "dendrimer_topology.png"):
        """Plots the topological graph of the dendrimer."""
        if not self.monomer_graph:
            print("Cannot plot topology. Generate the graph first by running .build()")
            return
            
        graph = nx.Graph()
        labels = {m.id: f"G{m.generation}" for m in self.monomer_graph}
        for m in self.monomer_graph:
            graph.add_node(m.id)
            for child in m.connections.values():
                graph.add_edge(m.id, child.id)

        pos = nx.kamada_kawai_layout(graph)
        node_colors = [m.generation for m in self.monomer_graph]

        fig, ax = plt.subplots(figsize=(10, 8))
        nx.draw(graph, pos, ax=ax, labels=labels, node_size=700,
                node_color=node_colors, cmap=plt.cm.viridis, font_size=8)
        
        ax.set_title("Dendrimer Topological Structure")
        plt.savefig(output_filename)
        print(f"Topology plot saved to '{output_filename}'")
        plt.close(fig)

    def _generate_graph(self):
        """Generates the abstract graph of connected monomer objects."""
        monomer_count = 0
        core = CustomMonomer(0, 0, self.core_data['smiles'], None, self.core_data['child_attachment_map_nums'])
        self.monomer_graph = [core]
        
        for gen in range(1, self.generations + 1):
            parents = [m for m in self.monomer_graph if m.generation == gen - 1]
            for parent in parents:
                if isinstance(parent, CustomMonomer):
                    for map_num in parent.child_attachment_map_nums:
                        monomer_count += 1
                        child = CustomMonomer(monomer_count, gen, **self.branch_data)
                        parent.connections[map_num] = child
                        self.monomer_graph.append(child)
        
        if self.surface_data:
            gen_to_cap = self.generations
            surface_parents = [m for m in self.monomer_graph if m.generation == gen_to_cap]
            for parent in surface_parents:
                if isinstance(parent, CustomMonomer):
                    for map_num in parent.child_attachment_map_nums:
                        monomer_count += 1
                        cap = SurfaceMonomer(monomer_count, parent.generation + 1, **self.surface_data)
                        parent.connections[map_num] = cap
                        self.monomer_graph.append(cap)
    
    def _build_and_optimize_mol(self, force_field: str, steps: int, ff_thr: float) -> Optional[pb.Molecule]:
        """Constructs, embeds, and optimizes the full molecule using pysoftk."""
        rw_mol = Chem.RWMol()
        
        # Build the molecular graph 
        for monomer in self.monomer_graph:
            mol_copy = Chem.Mol(monomer.mol)
            for atom in mol_copy.GetAtoms():
                if not atom.HasProp("molAtomMapNumber"):
                    new_idx = rw_mol.AddAtom(atom)
                    monomer.global_atom_indices[atom.GetIdx()] = new_idx
            for bond in mol_copy.GetBonds():
                b, e = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if b in monomer.global_atom_indices and e in monomer.global_atom_indices:
                    rw_mol.AddBond(monomer.global_atom_indices[b], monomer.global_atom_indices[e], bond.GetBondType())

        for parent in self.monomer_graph:
            if isinstance(parent, CustomMonomer):
                for map_num, child in parent.connections.items():
                    p_idx = parent.global_atom_indices[parent.child_attachment_atom_indices_map[map_num]]
                    c_idx = child.global_atom_indices[child.parent_connection_atom_idx]
                    rw_mol.AddBond(p_idx, c_idx, BondType.SINGLE)
        
        # Finalize graph and add hydrogens
        try:
            final_mol_rdkit = rw_mol.GetMol()
            Chem.SanitizeMol(final_mol_rdkit)
            final_mol_rdkit = Chem.AddHs(final_mol_rdkit, addCoords=True)
        except Exception as e:
            print(f"  Error during molecule finalization: {e}")
            return None

        # Generate a 3D conformer with RDKit
        try:
            print("    Embedding molecule to get initial 3D coordinates...")
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xf00d
            if AllChem.EmbedMolecule(final_mol_rdkit, params) == -1:
                 print("    ERROR: Conformer generation failed. Trying random coordinates.")
                 AllChem.EmbedMolecule(final_mol_rdkit, useRandomCoords=True)
        except Exception as e:
            print(f"    An exception occurred during RDKit embedding: {e}")
            return None

        # Convert to Pybel object for pysoftk utils functions
        mol_block = Chem.MolToMolBlock(final_mol_rdkit)
        mol_pybel = pb.readstring('mol', mol_block)

        # Optimize using pysoftk's ff_ob_relaxation function
        print(f"    Optimizing with pysoftk using {force_field} force field ({steps} steps, threshold {ff_thr})...")
        try:
            # The function ff_ob_relaxation takes steps as the 3rd argument (relax_iterations)
            optimized_mol_pybel = ff_ob_relaxation(mol_pybel, force_field, steps, ff_thr)
            print("    Optimization complete.")
        except Exception as e:
            print(f"    An exception occurred during pysoftk optimization: {e}")
            return None
            
        return optimized_mol_pybel


if __name__ == "__main__":
    # --- 1. Define Monomer Building Blocks ---
    core = {'smiles': "N([*:1])([*:2])[*:3]", 'child_attachment_map_nums': [1, 2, 3]}
    branch = {'smiles': "C(=O)N([*:1])CC(=O)N([*:2])CC(=O)N([*:3])", 
              'parent_attachment_map_num': 1, 'child_attachment_map_nums': [2, 3]}
    surface = {'smiles': "C(=O)O[*:1]", 'parent_attachment_map_num': 1}

    # --- 2. Initialize the Builder ---
    dendrimer_builder = DendrimerBuilder(
        core_data=core,
        branch_data=branch,
        generations=7,
        surface_data=surface
    )
    
    # --- 3. Run the Build with User-Defined Parameters ---
    # You can change these values to experiment with different optimizations
    dendrimer_builder.build(
        force_field="MMFF94", 
        steps=100, 
        ff_thr=0.001
    )
    
    # --- 4. Get the Final Molecule and Use It ---
    final_mol = dendrimer_builder.get_optimized_mol()
    
    if final_mol:
        # Save the final optimized structure to an XYZ file
        output_filename = f'dendrimer_G{dendrimer_builder.generations}_pysoftk_optimized.xyz'
        Fmt(final_mol).xyz_print(output_filename)
        print(f"\nSaved final 3D structure to '{output_filename}'")
        
        # Optional: Generate a topology plot
        dendrimer_builder.plot_topology()
    else:
        print("\nCould not generate the final molecule.")
