import shutil as sh
from pathlib import Path
from typing import Any

import pnab
from openbabel import openbabel as ob

from xna_builder.constants import BLOCKS_DICT
import xna_builder.utils as utils


PROJECT_ROOT = Path(__file__).parent.parent.resolve()
BACKBONES_DIR = PROJECT_ROOT / "data" / "backbones"
BASES_DIR = PROJECT_ROOT / "data" / "bases"


def config_pnab(
    bb: str, seq_base: str, pnab_config: Path, **pnab_kwargs: dict[str, Any]
) -> pnab.pNAB:
    """
    Creates and configures a pnab runner object for a single oligonucleotide strand.

    Args:
        bb (str): The backbone type key (e.g., 'pD').
        seq_base (str): The nucleobase sequence (e.g., 'ATGC').
        pnab_config (Path): Path to the pnab configuration (.yaml) file.
        **pnab_kwargs (Dict[str, Any]): Keyword arguments passed
            directly to PNAB.options["RuntimeParameters"].

    Returns:
        pnab.pNAB: Configured pnab runner object.
    """
    PNAB = pnab.pNAB(pnab_config)

    PNAB.options["Base methyl_cytosine"] = {
        "code": "M",
        "file_path": str(BASES_DIR / "5methylcytosine.pdb"),
        "linker": [2, 1],
        "name": "M",
        "pair_name": "G",
        "align": True,
    }

    bb_pdb = BLOCKS_DICT[bb]["bb_pdb"]
    bb_link = BLOCKS_DICT[bb]["bb_link"]
    bb_base_link = BLOCKS_DICT[bb]["bb_base_link"]

    PNAB.options["Backbone"] = {}
    PNAB.options["Backbone"]["file_path"] = str(BACKBONES_DIR / bb_pdb)
    PNAB.options["Backbone"]["interconnects"] = bb_link
    PNAB.options["Backbone"]["linker"] = bb_base_link
    PNAB.options["RuntimeParameters"]["strand"] = list(seq_base)

    for key, value in pnab_kwargs.items():
        PNAB.options["RuntimeParameters"][key] = value

    return PNAB


def build_sequence(PNAB: pnab.pNAB, out_pdb: Path) -> ob.OBMol:
    """
    Executes the pnab run process in a temporary directory 
    and reads the resulting PDB file into an OBMol object.

    Args:
        PNAB (pnab.pNAB): pnab runner object.
        out_pdb (Path): The final destination path for the generated PDB file.

    Returns:
        ob.OBMol: The Open Babel molecule object loaded from the generated PDB file.
    """
    seq = ob.OBMol()
    conv = ob.OBConversion()

    with utils.in_temp_dir():
        PNAB.run(number_of_cpus=1)
        result_pdb = Path(f"%i_%i.pdb" % (PNAB.results[0, 0], PNAB.results[0, 1]))

        sh.copy(result_pdb, out_pdb)

    conv.ReadFile(seq, str(out_pdb))

    return seq

