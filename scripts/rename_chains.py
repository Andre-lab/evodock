import sys
from Bio.PDB import PDBList, PDBIO, PDBParser

# pdbl = PDBList()

io = PDBIO()
parser = PDBParser()
# pdbl.retrieve_pdb_file(sys.argv[-1], pdir=".", file_format="pdb")

# pdb6gch.ent is the filename when retrieved by PDBList
structure = parser.get_structure(sys.argv[-1].replace(".pdb", ""), sys.argv[-1])

chains = [c.get_id() for c in structure.get_chains()]
print(chains)

renames = {chains[0]: "A", chains[1]: "B"}

for model in structure:
    for chain in model:
        old_name = chain.get_id()
        new_name = renames.get(old_name)
        if new_name:
            print(f"renaming chain {old_name} to {new_name}")
            chain.id = new_name
        else:
            print(f"keeping chain name {old_name}")

io.set_structure(structure)
io.save(sys.argv[-1].replace(".pdb", "_rename.pdb"))
