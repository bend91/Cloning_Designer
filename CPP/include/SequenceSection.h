#ifndef SEQUENCESECTION_H
#define SEQUENCESECTION_H

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include "SequenceData.h"
#include "PrimerData.h"

class SequenceSection: public SequenceData{
protected:
    std::vector <PrimerData> primer_list;
    std::vector <PrimerData> reverse_primer_list;
    // std::map <std::string, PrimerData> optimal_primers;
    std::vector <PrimerData> optimal_primers;
    std::vector <PrimerData> hifi_primers;
    float average_primer_tm;

    void identify_primers(){
        identify_primers_from_sequence(sequence, primer_list);
        identify_primers_from_sequence(reverse_complement, reverse_primer_list);
        averagePrimerTm();
    }

    void averagePrimerTm(){
        float average_f_primer_tm = calculateAverageTm(primer_list);
        float average_r_primer_tm = calculateAverageTm(reverse_primer_list);
        average_primer_tm = (average_f_primer_tm + average_r_primer_tm) / 2;
    }

    float calculateAverageTm(const std::vector <PrimerData>& primers){
        float total_tm = 0;
        for (std::size_t i = 0; i<primers.size(); i++){
            total_tm += primers[i].getTM();
        }
        return total_tm / primers.size();
    }

    void identify_primers_from_sequence(const std::string &seq, std::vector <PrimerData> &primers, int length_cutoff=40, float tm_min=55.0, float tm_max=70.0){
        for (int i=5; i<length_cutoff; i++){
            std::string primer_seq = seq.substr(0, i);
            int sequence_length = i - 1;
            int sequence_length_m1 = i - 2;
            if ((primer_seq[sequence_length] == 'G' || primer_seq[sequence_length] == 'C') && (primer_seq[sequence_length_m1] == 'T' || primer_seq[sequence_length_m1] == 'A')){
                PrimerData bind1("primer", primer_seq);
                float tm = bind1.getTM();
                if ((tm_min < tm) && (tm < tm_max)){
                    primers.push_back(bind1);
                } else if(tm > tm_max){
                    break;
                }
            }
        }
    }

public:
    SequenceSection(): SequenceData(){};
    SequenceSection(std::string n): SequenceData(n){}
    SequenceSection(std::string n, std::string seq): SequenceData(n, seq){
        identify_primers();
    }

    std::vector <PrimerData> getForwardPrimers() const{
        return primer_list;
    }
    std::vector <PrimerData> getReversePrimers() const{
        return reverse_primer_list;
    }


    float getAveragePrimerTm() const{
        return average_primer_tm;
    }

    void setSequence(std::string seq) override {
        SequenceData::setSequence(seq);
        identify_primers();
    }

    void calculateOptimalPrimers(float optimal_tm){
        optimal_primers.insert(optimal_primers.begin(), findOptimalPrimer(primer_list, optimal_tm));
        optimal_primers.insert(optimal_primers.begin() + 1, findOptimalPrimer(reverse_primer_list, optimal_tm));
    }

    // TODO - find the absolute difference
    PrimerData findOptimalPrimer(const std::vector <PrimerData>& primers, float optimal_tm){
        PrimerData optimal_primer;
        float tm_difference = 100;
        for (std::size_t i=0; i<primers.size(); i++){
            if (primers[i].getTM() - optimal_tm < tm_difference){
                tm_difference = std::abs(primers[i].getTM() - optimal_tm);
                optimal_primer = primers[i];
            }
        }
        return optimal_primer;
    }

    PrimerData getOptimalPrimer(std::string key){
        if (optimal_primers.size() < 2){
            throw std::runtime_error("Primers haven't been calculated yet");
        }
        else{
            switch(key[0]){
            case 'F': return optimal_primers[0];
            case 'R': return optimal_primers[1];
            default: throw std::invalid_argument("Invalid option for primer key, use Forward or Reverse\n");
            }

        }
    }

    void setHiFiPrimers(std::string forward_sequence, std::string reverse_sequence){
        std::string f_name = "F-" + getName();
        std::string r_name = "R-" + getName();
        PrimerData f_primer(f_name, forward_sequence);
        PrimerData r_primer(r_name, reverse_sequence);
        hifi_primers.insert(hifi_primers.begin(), f_primer);
        hifi_primers.insert(hifi_primers.begin() + 1, r_primer);
    }

    PrimerData getHiFiPrimer(std::string key){
        switch(key[0]){
        case 'F': return hifi_primers[0];
        case 'R': return hifi_primers[1];
        default: throw std::invalid_argument("Invalid option for primer key, use Forward or Reverse");
        }
    }
};

#endif
