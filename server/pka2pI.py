from pypka.clean.cleaning import get_pdb_Hs
from pdbmender.formats import read_pdb_line, read_pqr_line, new_pdb_line
from pdbmender.utils import mend_pdb, add_tautomers
from pdbmender.postprocess import fix_structure_states
import numpy as np
from copy import copy
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import os

anionic_aas = ["CTR", "GLU", "ASP", "CYS", "TYR"]
cationic_aas = ["NTR", "HIS", "ARG", "LYS"]

TITRABLETAUTOMERS = {
    "LYS": 3,
    "HIS": 2,
    "ASP": 4,
    "GLU": 4,
    "SER": 3,
    "THR": 3,
    "CYS": 3,
    "CTR": 4,
    "NTR": 3,
    "TYR": 2,
}


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

    pH_list = np.arange(pH_range[0], pH_range[1] + 1, pH_step).tolist()
    tit_x, tit_y = get_tit_curve(pKas, n_args, pH_list)
    return tit_x, tit_y


def get_tit_curve(pKas, n_args, pH_list):
    tit_x = []
    tit_y = []
    for pH in pH_list:

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
    if len(cys_atoms) == 0:
        return []
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
        # mindist_residues = -1
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
                    # mindist_residues = ii_cys

        if mindist < 3.5:
            to_exclude.append(i_cys)

    return to_exclude


def clean_pdb(f_pdb, new_f_pdb):
    atom_lines = []
    with open(f_pdb) as f:
        for line in f:
            if line.startswith("ATOM "):
                atom_lines.append(line)
            elif line.startswith("END"):
                break

    with open(new_f_pdb, "w") as f_new:
        f_new.write("".join(atom_lines))


def pkas_2_pdb(
    subID, inputpdbfilename, outputfilename, outputfilepH, outputfilenaming, results
):
    resnames = {}
    with open(inputpdbfilename) as f:
        for line in f:
            if line.startswith("ATOM "):
                cols = read_pdb_line(line)
                resname = cols[2]
                resnumb = cols[4]
                chain = cols[3]
                if chain not in resnames:
                    resnames[chain] = {}
                resnames[chain][resnumb] = resname

    terminal_offset = 5000

    sites = {}
    chains_res = {}
    for chain, resname, resnumb, pKa in results["pKas"]:
        termini_resname = ""
        if resname in ("NTR", "CTR"):
            resnumb += terminal_offset
            termini_resname = resnames[chain][resnumb]

        tit_x, tit_y = get_tit_curve(
            [[chain, resname, resnumb, pKa]], 0, [outputfilepH]
        )
        tit_curve = {x: y for x, y in zip(tit_x, tit_y)}

        charge = tit_curve[outputfilepH]
        if resname in cationic_aas:
            avg_prot = charge
            if avg_prot >= 0.5:
                state_prob = avg_prot
                most_prob_taut = TITRABLETAUTOMERS[resname] + 1
            else:
                state_prob = 1 - avg_prot
                most_prob_taut = 0 + 1
        elif resname in anionic_aas:
            avg_prot = 1 - charge * -1
            if avg_prot <= 0.5:
                state_prob = 1 - avg_prot
                most_prob_taut = TITRABLETAUTOMERS[resname] + 1
            else:
                state_prob = avg_prot
                most_prob_taut = 0 + 1

        taut_prob = -1.0

        if chain not in sites:
            sites[chain] = []
        sites[chain].append(
            (
                resname,
                resnumb,
                termini_resname,
                most_prob_taut,
                state_prob,
                taut_prob,
                {outputfilepH: avg_prot},
            )
        )
        if chain not in chains_res:
            chains_res[chain] = {}
        chains_res[chain][resnumb] = resname

    cys_bridges = {}
    cys_bridges_list = identify_cys_cys_bridges(inputpdbfilename)
    for chain, resnumb in cys_bridges_list:
        if chain not in cys_bridges:
            cys_bridges[chain] = []
        cys_bridges[chain].append(resnumb)

    ff_family = outputfilenaming
    if ff_family == "gromos_cph":
        ff_family = "GROMOS"
    logfile = "LOG_pdb2pqr"
    pdb2pqr_out = f"pdb2pqr_out_{subID}.pqr"
    mend_pdb(
        inputpdbfilename,
        pdb2pqr_out,
        ff_family,
        ff_family,
        logfile=logfile,
        hopt=False,
    )

    mainHs = get_pdb_Hs(pdb2pqr_out, chains_res)

    outputpqr = f"addhtaut_out_{subID}.pqr"
    add_tautomers(
        pdb2pqr_out,
        chains_res,
        ff_family,
        outputpqr,
        terminal_offset=terminal_offset,
    )

    last_anumb = 0
    with open(outputpqr) as f:
        delphi_input_content = []
        for line in f:
            if line.startswith("ATOM "):
                (
                    aname,
                    anumb,
                    resname,
                    chain,
                    resnumb,
                    x,
                    y,
                    z,
                    _,
                    _,
                ) = read_pqr_line(line)

            last_anumb += 1

            new_line = new_pdb_line(
                last_anumb, aname, resname, resnumb, x, y, z, chain=chain
            )
            delphi_input_content.append(new_line)
        delphi_input_content = delphi_input_content

    tit_atoms = []
    other_atoms = []
    for line in delphi_input_content:
        if line.startswith("ATOM "):
            line_parts = read_pdb_line(line)
            chain = line_parts[3]

            if chain in chains_res.keys():
                anumb = line_parts[1]
                resnumb = line_parts[4]
                if resnumb in chains_res[chain].keys():
                    tit_atoms.append(anumb)
                else:
                    other_atoms.append(anumb)

    # import pickle
    #
    # with open(f"/home/pedror/PypKa-Server-Back/server/{subID}.pickle", "wb") as handle:
    #     pickle.dump(
    #         [sites, tit_atoms],
    #         handle,
    #         protocol=pickle.HIGHEST_PROTOCOL,
    #     )
    import json

    with open("coco", "w") as f:
        f.write(
            str(subID)
            + "\n"
            + str(inputpdbfilename)
            + "\n"
            + str(outputfilename)
            + "\n"
            + str(outputfilepH)
            + "\n"
            + str(outputfilenaming)
            + "\n"
            + json.dumps(results)
        )

    fix_structure_states(
        outputfilename,
        outputfilepH,
        outputfilenaming,
        sites,
        tit_atoms,
        other_atoms,
        cys_bridges,
        delphi_input_content,
        [],
        terminal_offset,
        inputpdbfilename,
        mainchain_Hs=mainHs,
    )

