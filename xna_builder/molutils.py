from pathlib import Path
from itertools import product

from openbabel import openbabel as ob

def get_chirality_label(seq_linker: str) -> str:

    new_linker = "p" + seq_linker[1:].lower()  # First doesn't count
    chirality_label = "".join(char for char in new_linker if char in ("s", "r"))
    if not chirality_label:
        return ""
    return chirality_label


def get_diastereomers(config: dict) -> list[tuple[str, str, str]]:

    diastereomers = []

    seq_base = config["seq_base"]
    seq_sugar = config["seq_sugar"]
    seq_linker = "p" + config["seq_linker"][1:]

    uppercase_indexes = [i for i, char in enumerate(seq_linker) if char.isupper()]

    sr_positions = [
        i for i, char in enumerate(seq_linker) if char.lower() in {"s", "r"} and i != 0
    ]

    if not sr_positions or not config["diastereomers"]:
        diastereomers.append((seq_base, seq_sugar, seq_linker))
        return diastereomers

    if config["diastereomers"] == "all":
        for combination in product("sr", repeat=len(sr_positions)):
            new_linker = list(seq_linker)
            for pos, chirality in zip(sr_positions, combination):
                new_linker[pos] = chirality
            for i in uppercase_indexes:
                new_linker[i] = new_linker[i].upper()
            diastereomers.append((seq_base, seq_sugar, "".join(new_linker)))

    else:
        for chirality_label in config["diastereomers"]:
            new_linker = list(seq_linker)
            for pos, chirality in zip(sr_positions, chirality_label):
                new_linker[pos] = chirality
            for i in uppercase_indexes:
                new_linker[i] = new_linker[i].upper()
            diastereomers.append((seq_base, seq_sugar, "".join(new_linker)))

    return diastereomers


def renumber_atoms(pdb_in: Path, pdb_out: Path) -> None:

    mol = ob.OBMol()
    conv = ob.OBConversion()

    conv.ReadFile(mol, str(pdb_in))
    residues = sorted(ob.OBResidueIter(mol), key=lambda r: r.GetNum())
    idxs = [atom.GetIdx() for res in residues for atom in ob.OBResidueAtomIter(res)]

    mol.RenumberAtoms(idxs)

    conv.WriteFile(mol, str(pdb_out))


def remove_p5(olig: ob.OBMol) -> ob.OBMol:

    ATOMS_TO_STRIP = {"OP1", "OP2", "O1P", "O2P", "O5T", "H5T", "OP", "SP"}

    residues = list(ob.OBResidueIter(olig))

    first_res = min(residues, key=lambda r: r.GetNum())
    atoms = list(ob.OBResidueAtomIter(first_res))

    olig.BeginModify()
    for atom in atoms:
        atom_name = first_res.GetAtomID(atom).strip()
        if atom_name in ATOMS_TO_STRIP:
            olig.DeleteAtom(atom)
        elif atom_name == "P":
            atom.SetAtomicNum(1)
            atom.SetType("H")
            first_res.SetAtomID(atom, "HO5'")

    return olig
