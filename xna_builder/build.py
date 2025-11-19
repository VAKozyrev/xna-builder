import argparse
import shutil as sh
from pathlib import Path
from collections import defaultdict
from itertools import chain, product

import pnab
from openbabel import openbabel as ob

from xna_builder.constants import BLOCKS_DICT, RESNAMES_DICT
import xna_builder.utils as u


PROJECT_ROOT = Path(__file__).parent.parent.resolve()
BACKBONES_DIR = PROJECT_ROOT / "data" / "backbones"
BASES_DIR = PROJECT_ROOT / "data" / "bases"


def config_pnab_base(pnab_input: Path) -> pnab.pNAB:
    """Creates a pnab runner object with basic setup.

    Args:
        pnab_input (Path): Path to the pnab configuration .yaml file.

    Returns:
        pnab.pNAB: Configured pnab runner object.
    """
    PNAB = pnab.pNAB(pnab_input)

    PNAB.options["Base methyl_cytosine"] = {
        "code": "M",
        "file_path": f"{BASES_DIR}/5methylcytosine.pdb",
        "linker": [2, 1],
        "name": "M",
        "pair_name": "G",
        "align": True,
    }

    PNAB.options["Backbone"] = {}

    return PNAB


def build_1d_sequence(bb: str, seq_base: str, pnab_input: Path, seed: int) -> Path:
    """Generates a single-stranded oligonucleotide sequence with specified base sequence and backbone type.

    Args:
        bb (str): Backbone type, must be a key in `BLOCKS_DICT` dictionary.
        seq_base (str): Sequence of nucleobases.
        pnab_input (Path): Path to pnab configuration .yaml file.

    Returns:
        Path: Path to the generated .pdb file.
    """
    PNAB = config_pnab_base(pnab_input)

    bb_pdb = BLOCKS_DICT[bb]["bb_pdb"]
    bb_link = BLOCKS_DICT[bb]["bb_link"]
    bb_base_link = BLOCKS_DICT[bb]["bb_base_link"]

    PNAB.options["Backbone"]["file_path"] = f"{BACKBONES_DIR}/{bb_pdb}"
    PNAB.options["Backbone"]["interconnects"] = bb_link
    PNAB.options["Backbone"]["linker"] = bb_base_link

    PNAB.options["RuntimeParameters"]["strand"] = list(seq_base)
    PNAB.options["RuntimeParameters"]["build_strand"] = [
        True,  # Only 1 sequence
        True,
        False,
        False,
        False,
        False,
    ]

    PNAB.options["RuntimeParameters"]["seed"] = seed
    PNAB.run(number_of_cpus=1)

    out_pdb_path = Path(f"%i_%i.pdb" % (PNAB.results[0, 0], PNAB.results[0, 1]))

    return out_pdb_path


def combine_1d_sequences(
    residues_dict: dict[str, list[int]], sequences_dict: dict[str, Path]
) -> ob.OBMol:

    olig = ob.OBMol()
    conv = ob.OBConversion()

    residues = {}
    for bb in residues_dict.keys():
        seq = ob.OBMol()
        conv.ReadFile(seq, str(sequences_dict[bb]))
        for res_idx in residues_dict[bb]:
            res = ob.OBMol(seq)
            res.BeginModify()
            atoms = list(ob.OBMolAtomIter(res))
            for atom in atoms:
                if atom.GetResidue().GetNum() != res_idx:
                    res.DeleteAtom(atom)
            residues[res_idx] = res

    for res_idx in residues.keys():
        olig += residues[res_idx]

    old_atoms = list(
        chain.from_iterable(
            ob.OBMolAtomIter(residues[res_idx]) for res_idx in residues.keys()
        )
    )
    new_atoms = list(ob.OBMolAtomIter(olig))

    for source_atom, new_atom in zip(old_atoms, new_atoms):
        source_res = source_atom.GetResidue()
        new_res = new_atom.GetResidue()

        new_res.SetAtomID(new_atom, source_res.GetAtomID(source_atom))
        new_res.SetName(source_res.GetName())
        new_res.SetNum(source_res.GetNum())

    return olig


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


def get_chirality_label(seq_linker: str) -> str:

    new_linker = "p" + seq_linker[1:].lower()  # First doesn't count
    chirality_label = "".join(char for char in new_linker if char in ("s", "r"))
    if not chirality_label:
        return ""
    return chirality_label


def add_o3p_bonds(olig: ob.OBMol) -> ob.OBMol:

    residues = sorted(list(ob.OBResidueIter(olig)), key=lambda r: r.GetNum())

    def find_atom(residue: ob.OBResidue, name: str) -> ob.OBAtom | None:
        for atom in ob.OBResidueAtomIter(residue):
            if residue.GetAtomID(atom).strip() == name:
                return atom
        return None

    for res1, res2 in zip(residues[:-1], residues[1:]):
        o3 = find_atom(res1, "O3'")
        p = find_atom(res2, "P")
        if o3 and p:
            olig.AddBond(o3.GetIdx(), p.GetIdx(), 1)

    return olig


def remove_p5(olig: ob.OBMol) -> ob.OBMol:

    _ATOMS_TO_STRIP = {"OP1", "OP2", "O1P", "O2P", "O5T", "H5T", "OP", "SP", "H01"}

    residues = list(ob.OBResidueIter(olig))

    first_res = min(residues, key=lambda r: r.GetNum())
    atoms = list(ob.OBResidueAtomIter(first_res))

    olig.BeginModify()
    for atom in atoms:
        atom_name = first_res.GetAtomID(atom).strip()
        if atom_name in _ATOMS_TO_STRIP:
            olig.DeleteAtom(atom)
        elif atom_name == "P":
            atom.SetAtomicNum(1)
            atom.SetType("H")
            first_res.SetAtomID(atom, "HO5'")

    return olig


def rename_residues(olig: ob.OBMol, seq_base, seq_sugar, seq_linker) -> ob.OBMol:
    """
    Rename each residues in `olig` according by the provided base, sugar,
    and linker sequences.
    """
    seq_linker = seq_linker.replace("r", "s")
    if len(seq_linker) > 1 and seq_linker[1] == "s":
        seq_linker = "s" + seq_linker[1:]  # 5 terminus determided by the next linker
    residues = sorted(ob.OBResidueIter(olig), key=lambda r: r.GetNum())

    for sug, link, base, res in zip(seq_sugar, seq_linker, seq_base, residues):
        res_name = RESNAMES_DICT[f"{link.lower()}{sug}"] + base
        if link.isupper():
            res_name = res_name + "H"
        res.SetName(res_name)

    res5_name = residues[0].GetName().strip()
    if res5_name.endswith("H"):
        residues[0].SetName(res5_name[:-1] + "5")
    else:
        residues[0].SetName(res5_name + "5")

    # Handle the 3' terminus
    res3_name = residues[-1].GetName().strip()
    if res3_name.endswith("H"):
        residues[-1].SetName(res3_name[:-1] + "3")
    else:
        residues[-1].SetName(res3_name + "3")

    return olig



def build_oligonucleotide(job_config: dict) -> None:
    
    sequences_dict = {}  # Stores generated sequences

    PNAB_INPUT = "/home/lunet/cgvk2/xna-builder/data/DNA.yaml"

    conv = ob.OBConversion()

    residues_dict = defaultdict(list)
    name = job_config.get("name")
    seq_base: str = job_config.get("seq_base")
    seq_backbone: list[str] = job_config.get("seq_backbone")
    save_dir = Path(job_config.get("save_dir")).resolve()

    save_dir.mkdir(parents=True, exist_ok=True)

    for i, bb in enumerate(seq_backbone):
        residues_dict[bb].append(i + 1)

        for bb in residues_dict.keys():  
            if bb in sequences_dict:
                continue

            with u.in_temp_dir():

                source_file = build_1d_sequence(bb, seq_base, PNAB_INPUT, seed=42)
                destination_file =  save_dir / f"{bb}_{source_file}"
                sh.copy(source_file, destination_file)

            sequences_dict[bb] = destination_file
            
        # I.2. Generate oligomer by combining residues from different sequences
        combined_seq = combine_1d_sequences(residues_dict, sequences_dict)

        oligomer = add_o3p_bonds(combined_seq)
        #oligomer = isimo.remove_p5(oligomer)
        #oligomer = isimo.rename_residues(oligomer, seq_base, seq_sugar, seq_linker)

        oligomer_name = f"{name}.pdb"

        save_path = save_dir / oligomer_name
        conv.WriteFile(oligomer, str(save_path))


def main() -> None:
    """Main entry point for the script."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        "-c",
        type=Path,
        required=True,
        help="Path to the JSON configuration file.",
    )
    args = parser.parse_args()

    configs = u.read_config(args.config)

    for config in configs: 
        build_oligonucleotide(config)