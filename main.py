from functions import *


def main():
    fasta_filepath = input("Enter path to fasta file containing parts to clone together, in order: ")
    fasta_data = fasta_parser(fasta_filepath)
    primer_list = primer_designer(fasta_data)
    opt_primers = identify_optimal_annealing_temp(primer_list)
    hifi_primers = design_hifi_primers(fasta_data, opt_primers)
    print(hifi_primers)
    for i, data in enumerate(fasta_data):
        print(data.name)
        if i == len(fasta_data):
            break
        if i > 0:
            print("Forward: ", hifi_primers["Forward"][i-1])
        if i < len(fasta_data) - 1:
            print("Reverse: ", hifi_primers["Reverse"][i])
    # with open("cloning_plan.txt", "w") as f:
        # for i in fasta_data:
        # f.writelines(lines: Iterable[str])
    ...


if __name__ == "__main__":
    main()
