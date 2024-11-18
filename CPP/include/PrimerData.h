#ifndef PRIMERDATA_H
#define PRIMERDATA_H

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include "SequenceData.h"

class PrimerData: public SequenceData{
protected:
    float tm;
public:
    PrimerData(): SequenceData(){};
    PrimerData(std::string n): SequenceData(n){}
    PrimerData(std::string n, std::string seq): SequenceData(n, seq) {
        calculate_annealing_temp();
    }
    ~PrimerData();
    void calculate_annealing_temp(){
        if (sequence.length() > 13){
            tm = 64.9 + 41 * ((g_count + c_count) - 16.4) / sequence.length();
        } else {
            tm = ((a_count + t_count) * 2) + ((g_count + c_count) * 4);
        }
    };
    float getTM() const{
        return tm;
    }
};

#endif
