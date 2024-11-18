#ifndef SEQUENCEDATA_H
#define SEQUENCEDATA_H

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

class SequenceData{
public:
    SequenceData();
    SequenceData(std::string n){
        name = n;
    };
    SequenceData(std::string n, std::string seq){
        name = n;
        setSequence(seq);
    }
    ~SequenceData();
    void process_counts(){
        a_count = 0;
        c_count = 0;
        g_count = 0;
        t_count = 0;
        for (std::size_t i=0; i<sequence.size(); i++){
            switch(sequence[i]){
            case 'A': a_count++; break;
            case 'C': c_count++; break;
            case 'G': g_count++; break;
            case 'T': t_count++; break;
            default: break;
            }
        }
    }

    void processReverseComplement(){
        int i, j;
        reverse_complement.reserve(sequence.length());
        for (i=sequence.length() - 1; i >= 0; i--){
            char complement_base;
            switch(sequence[i]){
                case 'A': complement_base = 'T'; break;
                case 'C': complement_base = 'G'; break;
                case 'G': complement_base = 'C'; break;
                case 'T': complement_base = 'A'; break;
                default: throw std::invalid_argument("Invalid character in DNA sequence\n");
            }
            reverse_complement.push_back(complement_base);
        }
    }

    virtual void setSequence(std::string seq){
        int i;
        // Capitalise the sequence - lowercase ascii values are over 92 and are 32 places away from the uppercase ascii values
        for (i=0; i<seq.length(); i++){
            if (seq[i] > 92){
                seq[i] = seq[i] - 32;
            }
        }
        // set sequence to seq
        sequence = seq;
        // perform useful functions on the sequence
        process_counts();
        processReverseComplement();
    }
    int aCount() const{
        return a_count;
    }
    int tCount() const{
        return t_count;
    }
    int cCount() const{
        return c_count;
    }
    int gCount() const{
        return g_count;
    }
    std::string getSequence() const{
        return sequence;
    }

    std::string getSequenceEnd(int length){
        return sequence.substr(sequence.length() - length, length);
    }

    std::string getName() const{
        return name;
    }
    std::string getRCSequence() const{
        return reverse_complement;
    }

    std::string getRCSequenceEnd(int length){
        return reverse_complement.substr(reverse_complement.length() - length, length);
    }
protected:
    std::string sequence;
    std::string reverse_complement;
    std::string name;
    int a_count;
    int c_count;
    int g_count;
    int t_count;
};

#endif
