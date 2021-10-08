import numpy as np
import pandas as pd
from ortho_seq_code.constants_orthoseqs import *


def get_seq_info(seqf, alphbt_input, molecule):
    with open(seqf) as f:
        seq = f.readlines()
    seq_series_rm = pd.Series(seq).str.replace("\n", "")
    seq_series_nospace = seq_series_rm.str.replace(" ", "")
    seq_series = seq_series_nospace[seq_series_nospace != ""]
    if molecule == 'RNA'
        seq_series = pd.Series(seq).str.replace("T", "U")
    sites = max(seq_series.str.len())
    pop_size = len(seq_series)
    seq_list = list(np.unique(list("".join(list(seq_series)))))
    for i in seq_series:
        if i == "\n":
            pop_size -= 1
            seq.remove(i)
    # Autopadding lowercase n's
    if len(min(seq, key=len)) != len(max(seq, key=len)):
        incomplete_seq_series = seq_series[seq_series.str.len() < sites]
        while len(incomplete_seq_series) > 0:
            incomplete_seq_series = seq_series[seq_series.str.len() < sites]
            seq_series[incomplete_seq_series.index] += "n"
        seq_series += "\n"
        seq = list(seq_series)

    # Custom alphabet
    seq_oneline = "".join(seq_series)
    seq_list = list(seq_oneline)
    if alphbt_input is not None:
        alphbt = alphbt_input.upper()
        if alphbt == "PROTEIN_PNP":
            alphbt_input = "RNDEQHKSTY"
        elif alphbt == "ESSENTIAL":
            alphbt_input = "ILVFWHKTM"
        elif alphbt == "FREQUENCY_11AA":
            alphbt_input = "Y,G,D,V,S,A,F,R,L,PTWNEM"
        elif alphbt == "FREQUENCY_9AA":
            alphbt_input = "YF,G,DE,VLI,ST,A,RK,PWMN"
        elif alphbt == "ALBERTS":
            alphbt_input = "KRH,DE,AVLIPFMWGC"
        elif alphbt == "SIGMA":
            alphbt_input = "AILMV,FYW,NQCST,KRH,DE,G"
        elif alphbt == "HBOND":
            alphbt_input = "NQSTDERKYHW"
        elif alphbt == "HYDROPHOBICITY":
            alphbt_input = "LIFWVM,CYA,TEGSQD"
        if "," in alphbt_input:
            alphbt = alphbt_input.upper()
            # Adding on remaining letters as the last group
            alphbt_excluded = np.array([i for i in alphbt if i != ","])
            if "protein" in molecule:
                alphbt_last_group = "".join(
                    np.setdiff1d(
                        np.array(PROTEIN_ALPHABETS).ravel(), np.array(alphbt_excluded)
                    )
                )
            elif "RNA" in molecule:
                alphbt_last_group = "".join(
                    np.setdiff1d(
                        np.array(RNA_ALPHABETS).ravel(), np.array(alphbt_excluded)
                    )
                )
            else:
                alphbt_last_group = "".join(
                    np.setdiff1d(
                        np.array(DNA_ALPHABETS).ravel(), np.array(alphbt_excluded)
                    )
                )
            alphbt += "," + str(alphbt_last_group)
            custom_aa = alphbt.split(",")
            if "" in custom_aa:
                custom_aa.remove("")
            # Assign group names to the group
            alphbt_count = 0
            aa_dict = dict()
            for i in range(len(custom_aa)):
                aa_dict[str(alphbt_count)] = list(np.unique(list(custom_aa[i])))
                if aa_dict[str(alphbt_count)] == []:
                    del aa_dict[str(alphbt_count)]
                    alphbt_count -= 1
                alphbt_count += 1
            if "n" in seq_list:
                aa_dict[str(alphbt_count)] = ["n"]
                custom_aa.append("n")
                alphbt_count += 1
            # Replaces amino acids with groups
            for i in range(len(seq_list)):
                for j in range(alphbt_count):
                    if seq_list[i] in aa_dict[str(j)]:
                        seq_list[i] = str(list(aa_dict.keys())[j])
            seq_list_sub = seq_list
            alphabets = list(aa_dict.keys())

        else:
            alphabets = sorted(list(alphbt_input))
            alphabets_other = np.setdiff1d(np.array(seq_list), np.array(alphabets))
            if len(alphabets_other) > 0 and list(alphabets_other) != ["n"]:
                if "protein" in molecule:
                    alphbt_last_group = list(
                        np.setdiff1d(
                            np.array(PROTEIN_ALPHABETS).ravel(), np.array(alphabets)
                        )
                    )
                elif "RNA" in molecule:
                    alphbt_last_group = list(
                        np.setdiff1d(
                            np.array(RNA_ALPHABETS).ravel(), np.array(alphabets)
                        )
                    )
                else:
                    alphbt_last_group = list(
                        np.setdiff1d(
                            np.array(DNA_ALPHABETS).ravel(), np.array(alphabets)
                        )
                    )
                seq_list_sub = []
                for i in range(len(seq_list)):
                    if seq_list[i] in alphbt_last_group:
                        seq_list_sub.append("z")
                    else:
                        seq_list_sub.append(seq_list[i])
            if "z" in seq_list_sub and "z" not in alphabets:
                alphabets.append("z")
            if "n" in seq_list_sub and "n" not in alphabets:
                alphabets.append("n")
            custom_aa = alphabets
        seq_adj = "".join(seq_list_sub)
        seq = [seq_adj[i : i + sites] for i in range(0, len(seq_adj), sites)]
    else:
        alphabets = list(np.unique(seq_list))
        custom_aa = None
    while "" in alphabets:
        alphabets.remove("")
    while " " in alphabets:
        alphabets.remove(" ")
    while "\n" in alphabets:
        alphabets.remove("\n")
    dm = len(alphabets)
    return [dm, sites, pop_size, seq, seq_series, alphabets, custom_aa]
