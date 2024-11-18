#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <cmath>
#include <fstream>
#include "SequenceData.h"
#include "PrimerData.h"
#include "SequenceSection.h"
#include "FastaFile.h"


float averageSequencePrimerTm(std::vector <SequenceSection> sequences){
    float total_tm = 0;
    for (std::size_t i=0; i<sequences.size(); i++){
        total_tm += sequences[i].getAveragePrimerTm();
    }
    return total_tm / sequences.size();
}


int main(){
    int i, num_t;
    float optimal_tm = 0.0;
    std::string filename, f_primer_sequence, r_primer_sequence;
    std::string save_file = "cloning_plan.txt";
    std::vector<PrimerData> primer_list;
    std::vector<SequenceSection> data;
    std::ofstream CloningPlan(save_file);
    FastaFile *fastafile;

    std::cout << "Enter path to fasta file containing parts to clone together, in order (type \"test\" to run a test: " << std::endl;
    // std::cin reads up to the first whitespace character, whereas getline will read the whole line
    std::getline(std::cin, filename);
    std::cout << "filename: " << filename << std::endl;
    if (filename == "test"){
        filename = "../test_data/test.fasta";
    }
    fastafile = new FastaFile(filename);
    data = fastafile->getData();
    std::cout << "Sequence parts of the fasta file: " << std::endl;
    for (i=0; i<data.size(); i++){
        std::cout << data[i].getName() << ": " << data[i].getSequence() << std::endl;
    }
    std::cout << std::endl;
    i=0;
    optimal_tm =  averageSequencePrimerTm(data);


    for (i=0; i<data.size(); i++){
        CloningPlan << data[i].getName() << std::endl;
        data[i].calculateOptimalPrimers(optimal_tm);
        if (i < data.size() - 1){
            r_primer_sequence = data[i+1].getRCSequenceEnd(8) + data[i].getOptimalPrimer("Reverse").getSequence();
        } else{
            r_primer_sequence = "";
        }
        if (i > 0){
            f_primer_sequence = data[i-1].getSequenceEnd(8) + data[i].getOptimalPrimer("Forward").getSequence();
        } else {
            f_primer_sequence = "";
        }

        data[i].setHiFiPrimers(f_primer_sequence, r_primer_sequence);
        CloningPlan << "- Forward: " << data[i].getHiFiPrimer("Forward").getSequence() << std::endl;
        CloningPlan << "- Reverse: " << data[i].getHiFiPrimer("Reverse").getSequence() << std::endl;
        CloningPlan << std::endl;
    }
    CloningPlan << "Optimal annealing temperature is: " << optimal_tm + 5.0 << "oC" << std::endl;
    CloningPlan.close();
    std::cout << "Cloning plan has been saved at " << save_file << std::endl;
    return 0;
}
