import argparse
from pathlib import Path
from collections import defaultdict

from openbabel import openbabel as ob

from xna_builder import OligoConfig, config_pnab, build_sequence
from xna_builder.molutils import renumber_atoms
from xna_builder.constants import RESNAMES_DICT


def combine_sequences(
    residues_dict: dict[str, list[int]], sequences_dict: dict[str, ob.OBMol]
) -> ob.OBMol:

    olig = ob.OBMol()

    for bb in residues_dict.keys():
        seq = sequences_dict[bb]
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

    configs = OligoConfig.from_json(args.config)

    conv = ob.OBConversion()

    for config in configs:

        config.save_dir.mkdir(parents=True, exist_ok=True)
        olig_pdb = config.save_dir / f"{config.name}.pdb"

        residues_dict = defaultdict(list)  # Maps backbone types to res idxs
        sequences_dict = {}  # Maps sequences pdb to res idxs

        for i, bb in enumerate(config.seq_backbone, start=1):
            residues_dict[bb].append(i)

        for bb in residues_dict.keys():
            if bb in sequences_dict:
                continue

            PNAB = config_pnab(
                bb, config.seq_base, config.pnab_config, **config.pnab_kwargs
            )
            sequences_dict[bb] = build_sequence(PNAB, config.save_dir / f"{bb}_seq.pdb")

        oligomer = combine_sequences(residues_dict, sequences_dict)
        oligomer = add_o3p_bonds(oligomer)
        oligomer = rename_residues(oligomer, config.seq_backbone, config.seq_base)
    
        conv.WriteFile(oligomer, str(olig_pdb))

        renumber_atoms(olig_pdb, olig_pdb)


if __name__ == "__main__":
    main()
