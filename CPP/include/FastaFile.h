#ifndef FASTAFILE_H
#define FASTAFILE_H

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include "SequenceSection.h"

class FastaFile{
public:
    FastaFile();
    FastaFile(std::string fname){
        filename = fname;
        data = readFastaFile(filename);
        int i;
        for (i=0; i<data.size(); i++){
            headers.push_back(data[i].getName());
            content.push_back(data[i].getSequence());
        }
    };
    std::string getFileName(){
        return filename;
    };
    std::vector<SequenceSection> getData(){
        return data;
    }
    std::vector<std::string> getHeaders(){
        return headers;
    }
    std::vector<std::string> getContent(){
        return content;
    }
    std::vector<SequenceSection> readFastaFile(const std::string fname){
        std::ifstream fasta_file(fname);
        if (!fasta_file.is_open()){
            std::cerr << "Error, can't open file" << std::endl;
        }
        std::vector <SequenceSection> data;
        SequenceSection *sequence_data;
        std::string line;
        while (std::getline(fasta_file, line)){
            if (line.empty()){continue;}
            if (line[0] == '>'){
                sequence_data = new SequenceSection(line);
            } else if (line[0] != '\n') {
                sequence_data->setSequence(line);
                data.push_back(*sequence_data);
            } else { };
        };
        return data;
    };
protected:
    std::string filename;
    std::vector <SequenceSection> data;
    std::vector<std::string> content;
    std::vector<std::string> headers;
};
#endif
