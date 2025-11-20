import argparse
from pathlib import Path
from collections import defaultdict
from itertools import product

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
        False,
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

    for bb in residues_dict.keys():
        seq = ob.OBMol()
        conv.ReadFile(seq, str(sequences_dict[bb]))
        bv = ob.OBBitVec(seq.NumAtoms() + 1)

        for atom in ob.OBMolAtomIter(seq):
            if atom.GetResidue().GetNum() in residues_dict[bb]:
                bv.SetBitOn(atom.GetIdx())

        seq.CopySubstructure(olig, bv, None)

    return olig


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
        if o3 and p and not olig.GetBond(o3, p):
            olig.AddBond(o3.GetIdx(), p.GetIdx(), 1)

    return olig


def rename_residues(olig: ob.OBMol, seq_backbone: list[str], seq_base: str) -> ob.OBMol:
    """
    Rename each residues in olig according by the provided base and backbone sequences
    """
    residues = sorted(ob.OBResidueIter(olig), key=lambda r: r.GetNum())

    for res, bb, base in zip(residues, seq_backbone, seq_base):
        res_name = RESNAMES_DICT[bb] + base
        res.SetName(res_name)

    res5_name = residues[0].GetName().strip()
    residues[0].SetName(res5_name + "5")

    res3_name = residues[-1].GetName().strip()
    residues[-1].SetName(res3_name + "3")

    return olig


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


def build_oligonucleotide(
    seq_base: str, seq_backbone: list[str], pnab_input: Path
) -> ob.OBMol:

    sequences_dict = {}
    residues_dict = defaultdict(list)

    for i, bb in enumerate(seq_backbone):
        residues_dict[bb].append(i + 1)

    for bb in residues_dict.keys():
        if bb in sequences_dict:
            continue

        with u.in_temp_dir():

            seq_pdb = build_1d_sequence(bb, seq_base, pnab_input, seed=42)
            sequences_dict[bb] = seq_pdb

    oligomer = combine_1d_sequences(residues_dict, sequences_dict)
    oligomer = add_o3p_bonds(oligomer)
    oligomer = rename_residues(oligomer, seq_backbone, seq_base)

    return oligomer


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

    conv = ob.OBConversion()

    for config in configs:
        name = config.get("name")
        seq_base = config.get("seq_base")
        seq_backbone = config.get("seq_backbone")
        save_dir = Path(config.get("save_dir"))

        save_dir.mkdir(parents=True, exist_ok=True)

        pnab_input = config.get("pnab_input", PROJECT_ROOT / "data" / "DNA.yaml")

        oligomer = build_oligonucleotide(seq_base, seq_backbone, pnab_input)

        conv.WriteFile(oligomer, str(save_dir / f"{name}.pdb"))
