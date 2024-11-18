



// std::vector<std::vector <std::string> > readFastaFileVector(const std::string fname){
//     std::ifstream fasta_file(fname);
//     if (!fasta_file.is_open()){
//         std::cerr << "Error, can't open file" << std::endl;
//     }
//     std::vector <std::vector<std::string> > data;
//     std::vector <std::string> headers;
//     std::vector <std::string> content;
//     std::string line;
//     while (std::getline(fasta_file, line)){
//         if (line.empty()){continue;}
//         if (line[0] == '>'){
//             headers.push_back(line);
//         } else if (line[0] != '\n') {
//             content.push_back(line);
//         } else { };
//     };
//     data.push_back(headers);
//     data.push_back(content);

//     return data;
// };





