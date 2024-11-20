from functions import *


strategy_list = [
    "HiFi"
]

strategy = {
    "HiFi": design_hifi_primers,
}


def main():
    fasta_filepath = input("Enter path to fasta file containing parts to clone together, in order, otherwise type test to use the test data: ")
    if fasta_filepath == "test":
        fasta_data = fasta_parser('test_data/test.fasta')
    else:
        fasta_data = fasta_parser(fasta_filepath)


    primer_list = primer_designer(fasta_data)
    if input("Define a custom annealing temperature? If no then annealing temperature will be based on the least variation between primer sequences\n (Y/N): ") in "Yesyes":
        opt_tm = float("Enter annealing temperature to aim for: ")
    else:
        opt_tm = None
    opt_primers = identify_optimal_annealing_temp(primer_list, opt_tm)
    print("Enter the strategy you would like to use for cloning: \n")
    [print(f"[{i}] {x}\n") for i, x in enumerate(strategy_list)]
    chosen_strategy = strategy_list[int(input("> "))]
    designed_primers = strategy[chosen_strategy](fasta_data, opt_primers)

    if input("Enter filepath to write plan to (otherwise will just be saved in this directory as cloning_plan.txt) (Y/N)? ") in "Yesyes":
        save_file_dir = input("Enter path to save the file to: ")
        save_file_name = input("Enter the name you would like for the file: ")
        if save_file_dir[0] != "/":
            save_file_dir = save_file_dir + "/"
        save_file = save_file_dir + save_file_name
    else:
        save_file = "cloning_plan.txt"

    with open(save_file, "w") as f:
        for i, data in enumerate(fasta_data):
            f.write(f"{data.name}\n")
            if i == len(fasta_data):
                break
            if i > 0:
                f.write(f" - Forward: {designed_primers['Forward'][i-1].sequence}\n")
            if i < len(fasta_data) - 1:
                f.write(f" - Reverse: {designed_primers['Reverse'][i].sequence}\n")
            f.write("\n")
        f.write(f"Optimal PCR annealing temperature is {round(sum([x.tm for x in opt_primers['Forward']] + [x.tm for x in opt_primers['Reverse']]) / len(opt_primers['Forward'] + opt_primers['Reverse']), 2) + 5}\n")
    ...


if __name__ == "__main__":
    main()
