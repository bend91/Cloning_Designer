from collections import namedtuple as nt
import re
import matplotlib.pyplot as plt


fasta_nt = nt("FastaData", ["name", "sequence"])
primer_nt = nt("PrimerData", ["sequence", "tm"])


def fasta_parser(filepath)->list:
    """
    Parses data in a file in fasta format, setting the header defined with a ">" as the name and the content immediately under as the sequence in a named tuple format.
    Currently sets data until another header is found, will result in multiple named_tuples with the same name but different sequence if there are "\n" values in the sequence
    """
    data_list = []
    with open(filepath, "r") as f:
        a = f.readline()
        while a:
            if a[0] == ">":
                header = a[1:].strip()
            elif a[0] != "\n":
                data_list.append(fasta_nt(header, a.strip()))
            a = f.readline()
    return data_list


def calculate_annealing_temperature(seq)->float:
    """
    Calculates the estimated annealing temperature of the given DNA sequence to a complementary strand
    :param seq: string of a DNA sequence comprising A, C, G, T
    """
    seq = seq.upper()
    if len(seq) > 13:
        tm = 64.9 + 41 * (seq.count("G") + seq.count("C") - 16.4) / len(seq)
    else:
        tm = (seq.count("A") + seq.count("T")) * 2 + (seq.count("G") + seq.count("C")) * 4
    return round(tm, 2)


def reverse_complement(seq)->str:
    seq = seq.lower()
    return seq[::-1].replace('a', "T").replace('c', "G").replace('g', "C").replace('t', "A")


def identify_primers(seq, length_cutoff: int=40, tm_min: float=55.0, tm_max: float=70.0)->list:
    """
    Identifies potential primer sequences of the given sequence defined by a maximum length cut off and a minimum and maximum annealing temperature
    :param seq: string of DNA sequence containing only A, C, G, T
    :param length_cutoff: integer defining the maximum length of the primer sequence
    :param tm_min: float defining the minimum annealing temperature for a primer sequence to be valid
    :param tm_max: float defining the maximum annealing temperature for a primer sequence to be valid
    """
    seq = seq.upper()
    primer_list = []
    for i in range(length_cutoff):
        bind1 = seq[:i]
        binding = re.match("[A|T][G|C]", bind1[-2:])
        tm = calculate_annealing_temperature(bind1)
        if binding and (tm_max > tm > tm_min):
            primer_list.append(primer_nt(bind1, tm))
        elif tm > tm_max:
            if len(primer_list) == 0:
                print("No suitable primer sequences found with these parameters, consider increasing length_cutoff or tm_max")
            break
    return primer_list


def primer_designer(data_list)->list:
    """
    Identifies a list of primers that bind to the sense and antisense of the given list of sequences and returns a list of a dictionary with Forward and Reverse primers
    :param data_list: list of fasta_nt namedtuples
    """
    seq1 = data_list[0].sequence
    re_seq = reverse_complement(seq1)
    primer_list1 = identify_primers(re_seq)
    primer_list = [{"Reverse": primer_list1}]
    for i, seq in enumerate(data_list[1:], 1):
        primer_list.append({
            "Forward": identify_primers(seq.sequence),
            "Reverse": identify_primers(reverse_complement(seq.sequence))
            })
    return primer_list
    ...


def identify_optimal_annealing_temp(primer_list, opt_tm: float=None):
    """
    Takes a list of dicts of lists of primers with each lower level list being equivalent primers with varying sequences and tms and identifies a set of primers, one from each list, with the most similar tms
    :param primer_list: list of dicts of primers
    :opt_tm: define a custom optimal tm to aim for, if this is set then there will be more variability in primer annealing temp.
    """
    if opt_tm is None:
        min_tms = []
        for primers in primer_list:
            for x in list(primers.keys()):
                min_tms.append(min(primer.tm for primer in primers[x]))
        max_min_tm = max(min_tms)
    else:
        max_min_tm = opt_tm
    opt_primers = {"Forward": [], "Reverse": []}
    for i, primers in enumerate(primer_list):
        for key in list(primers.keys()):
            opt_primer = min(primers[key], key=lambda p: abs(p.tm - max_min_tm))
            opt_primers[key].append(opt_primer)
    return opt_primers
    ...


def design_hifi_primers(fasta_data: list, primer_list: list, forward_key: str="Forward", reverse_key: str="Reverse"):
    hifi_primers = {forward_key: [], reverse_key: []}
    i = 0
    while i + 1 < len(fasta_data):
        fprimer = primer_list[forward_key][i]
        rprimer = primer_list[reverse_key][i]
        hifi_primers[forward_key].append(i)
        hifi_primers[reverse_key].append(i)

        fhifi_sequence = fasta_data[i].sequence[-8:] + fprimer.sequence
        rhifi_sequence = reverse_complement(fasta_data[i+1].sequence)[-8:] + rprimer.sequence

        hifi_primers[forward_key][i] = primer_nt(fhifi_sequence, calculate_annealing_temperature(fhifi_sequence))
        hifi_primers[reverse_key][i] = primer_nt(rhifi_sequence, calculate_annealing_temperature(rhifi_sequence))

        i += 1
    return hifi_primers

