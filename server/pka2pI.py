from pdbmender.formats import read_pdb_line
import numpy as np
from copy import copy
import pandas as pd
from scipy.spatial.distance import pdist, squareform

anionic_aas = ["CTR", "GLU", "ASP", "CYS", "TYR"]
cationic_aas = ["NTR", "HIS", "ARG", "LYS"]


def pkas_2_titcurve(pdbfile, pKas, pH_range=[-20, 20], pH_step=0.5):
    n_args = 0
    chains = []
    with open(pdbfile) as f:
        for line in f:
            if line.startswith("ATOM "):
                cols = read_pdb_line(line)
                atomname = cols[0]
                resname = cols[2]
                chain = cols[3]

                if atomname == "C":
                    if resname == "ARG":
                        n_args += 1

                if chain not in chains:
                    chains.append(chain)

    pKas = copy(pKas)
    for chain in chains:
        pKas += [
            ["", "NTR", "", 7.99],
            ["", "CTR", "", 2.9],
        ]

    tit_x = []
    tit_y = []
    for pH in np.arange(pH_range[0], pH_range[1] + 1, pH_step):

        charge = n_args
        for tit_res in pKas:
            _, resname, _, res_pka = tit_res
            partial_occ = 1.0 / (10 ** (pH - res_pka) + 1.0)
            if resname in cationic_aas:
                partial_charge = partial_occ
            elif resname in anionic_aas:
                partial_charge = (1 - partial_occ) * -1

            charge += partial_charge

        tit_x.append(pH)
        tit_y.append(round(charge, 2))

    return tit_x, tit_y


def titcurve_2_pI(tit_x, tit_y):
    titration_curve = [(pH, tit_y[i]) for i, pH in enumerate(tit_x)]
    i = 0
    pH, prot = titration_curve[i]
    while prot > 0.0 and i < len(titration_curve) - 1:
        i += 1
        pH, prot = titration_curve[i]
    if i == 0 or i == len(titration_curve) - 1:
        return pH
    else:
        last_pH, last_prot = titration_curve[i - 1]
        isoelectric_point = pH - (prot * (last_pH - pH)) / (last_prot - prot)
        return isoelectric_point


def exclude_cys(pdbfile, pkas):
    to_exclude = identify_cys_cys_bridges(pdbfile)

    final_pkas = []

    for res in pkas:
        if res[1] == "CYS":
            add_trigger = True
            for exclude_cys in to_exclude:
                if exclude_cys[0] == res[0] and exclude_cys[1] == res[2]:
                    add_trigger = False
                    break
            if add_trigger:
                final_pkas.append(res)
        else:
            final_pkas.append(res)

    return final_pkas


def identify_cys_cys_bridges(pdbfile):
    cys_atoms = []
    with open(pdbfile) as f:
        for line in f:
            if line.startswith("ATOM "):
                cols = read_pdb_line(line)
                resname = cols[2]
                if resname == "CYS":
                    cys_atoms.append(cols)
    df_cys = pd.DataFrame(
        cys_atoms,
        columns=["anumb", "aname", "resname", "chain", "resnumb", "x", "y", "z"],
    )

    cys_res = (
        df_cys[["chain", "resnumb"]]
        .value_counts()
        .reset_index(name="count")[["chain", "resnumb"]]
        .values.tolist()
    )

    dists = squareform(pdist(df_cys[["x", "y", "z"]]))

    to_exclude = []
    for i, i_cys in enumerate(cys_res):
        i_indices = df_cys.query(
            f"chain == '{i_cys[0]}' and resnumb == {i_cys[1]}"
        ).index.values.tolist()

        mindist = 9999
        mindist_residues = -1
        for atom_i in i_indices:
            dists_i = dists[atom_i]
            for ii in range(len(cys_res)):
                if ii == i:
                    continue
                ii_cys = cys_res[ii]
                ii_indices = df_cys.query(
                    f"chain == '{ii_cys[0]}' and resnumb == {ii_cys[1]}"
                ).index.values.tolist()

                dists_i_atom_vs_ii_cys = min(dists_i[ii_indices])

                if dists_i_atom_vs_ii_cys < mindist:
                    mindist = dists_i_atom_vs_ii_cys
                    mindist_residues = ii_cys

        if mindist < 3.5:
            to_exclude.append(i_cys)

    return to_exclude


if __name__ == "__main__":
    pkas = [
        ["A", "LYS", 1, 10.26],
        ["A", "CYS", 6, 10.83],
        ["A", "GLU", 7, 3.56],
        ["A", "LYS", 13, 11.16],
        ["A", "HIS", 15, 6.39],
        ["A", "ASP", 18, 2.92],
        ["A", "TYR", 20, 10.11],
        ["A", "TYR", 23, 9.85],
        ["A", "CYS", 30, 11.87],
        ["A", "LYS", 33, 10.47],
        ["A", "GLU", 35, 4.1],
        ["A", "ASP", 48, 2.47],
        ["A", "ASP", 52, 2.31],
        ["A", "TYR", 53, 11.86],
        ["A", "CYS", 64, 13.02],
        ["A", "ASP", 66, 1.94],
        ["A", "CYS", 76, 12.91],
        ["A", "CYS", 80, 12.25],
        ["A", "ASP", 87, 2.77],
        ["A", "CYS", 94, 14.06],
        ["A", "LYS", 96, 11.03],
        ["A", "LYS", 97, 11.14],
        ["A", "ASP", 101, 3.46],
        ["A", "CYS", 115, 14.18],
        ["A", "LYS", 116, 10.13],
        ["A", "ASP", 119, 2.9],
        ["A", "CYS", 127, 11.55],
    ]
    tit_x, tit_y = pkas_2_titcurve(
        "pdbs/202211240952597681597710917568829924855.pdb", pkas
    )

    pI = titcurve_2_pI(tit_x, tit_y)

    exclude_cys("pdbs/202211240952597681597710917568829924855.pdb", pkas)
