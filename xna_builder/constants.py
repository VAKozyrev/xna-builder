"""
Module contains contants used across the isimo project.
"""

SUGARS = ("D", "F", "X", "L", "M", "E", "C")
BASES = ("A", "C", "T", "G", "U", "M")
LINKERS = ("p", "s", "r")

BLOCKS_DICT = {
    "pD": {
        "bb_pdb": "DNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [12, 13],
    },
    "rD": {
        "bb_pdb": "RpDNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [12, 13],
    },
    "sD": {
        "bb_pdb": "SpDNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [12, 13],
    },
    "pX": {
        "bb_pdb": "RNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "sX": {
        "bb_pdb": "SpRNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "rX": {
        "bb_pdb": "RpRNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "pL": {
        "bb_pdb": "LNA_backbone.pdb",
        "bb_link": [11, 1],
        "bb_base_link": [6, 7],
    },
    "rL": {
        "bb_pdb": "RpLNA_backbone.pdb",
        "bb_link": [11, 1],
        "bb_base_link": [6, 7],
    },
    "sL": {
        "bb_pdb": "SpLNA_backbone.pdb",
        "bb_link": [11, 1],
        "bb_base_link": [6, 7],
    },
    "pM": {
        "bb_pdb": "2MO_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "rM": {
        "bb_pdb": "Rp2MO_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "sM": {
        "bb_pdb": "Sp2MO_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "pE": {
        "bb_pdb": "MOE_backbone.pdb",
        "bb_link": [1, 31],
        "bb_base_link": [6, 30],
    },
    "rE": {
        "bb_pdb": "RpMOE_backbone.pdb",
        "bb_link": [1, 31],
        "bb_base_link": [6, 30],
    },
    "sE": {
        "bb_pdb": "SpMOE_backbone.pdb",
        "bb_link": [1, 31],
        "bb_base_link": [6, 30],
    },
    "pC": {
        "bb_pdb": "CET_backbone.pdb",
        "bb_link": [11, 1],
        "bb_base_link": [6, 7],
    },
    "rC": {
        "bb_pdb": "RpCET_backbone.pdb",
        "bb_link": [11, 1],
        "bb_base_link": [6, 7],
    },
    "sC": {
        "bb_pdb": "SpCET_backbone.pdb",
        "bb_link": [11, 1],
        "bb_base_link": [6, 7],
    },
    "pF": {
        "bb_pdb": "FNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "sF": {
        "bb_pdb": "SpFNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
    "rF": {
        "bb_pdb": "RpFNA_backbone.pdb",
        "bb_link": [10, 1],
        "bb_base_link": [13, 14],
    },
}

RESNAMES_DICT = {
    "pD": "D",
    "rD": "S",
    "sD": "S",
    "pX": "X",
    "rX": "R",
    "sX": "R",
    "pL": "L",
    "rL": "N",
    "sL": "N",
    "pM": "M",
    "rM": "T",
    "sM": "T",
    "pE": "E",
    "rE": "O",
    "sE": "O",
    "pC": "C",
    "rC": "I",
    "sC": "I",
    "pF": "F",
    "rF": "V",
    "sF": "V",
}

FRAGMENTS = {
    "a": ("O3'", "C3'"),
    "a-B": ("O3'", "C3'"),
    "b": ("P", "O3'"),
    "c": ("P", "O5'"),
    "d": ("O5'", "C5'"),
    "w": ("O3'", "C3'"),
    "x": ("P", "O3'"),
    "y": ("P", "O5'"),
    "z": ("O5'", "C5'"),
}
# DIHEDRALS = {
#     "alpha": O3' (i-1) - P (i) - 05' (i) - C5' (i)
#     \item$\beta$: P (i) - 05' (i) - C5' (i) - C4' (i)
#     "gamma": O5' (i) - C5' (i) - C4' (i) - C3' (i)
#     \item$\delta$: C5' (i) - C4' (i) - C3'(i) - O3' (i)
#     \item$\epsilon$: C4' (i) - C3'(i) - O3' (i) - P (i+1)
#     \item$\zeta$: C3'(i) - O3' (i) - P (i+1) - O5' (i+1)
# }


#     nu0: C4'- 04' - C1'- C2'
#     nu1: O4' - C1' - C2' - C3'
#     nu2: C1' - C2' - C3' - C4'
#     nu3: C2' - C3' - C4' - 04'
#     nu4: C3' - C4' - O4' - C1'


# The orientation of the nucleobase relative to the sugar ring is defined by the glycosidic torsion angle, $\chi$ (chi). Its definition depends on the type of nucleobase:
# Purines (Adenine, Guanine): $\chi$ is defined by the atoms O4' - C1' - N9 - C4. Pyrimidines(Cytosine, Thymine, Uracil): $\chi$ is defined by the atoms O4' - C1' - N1 - C2.
