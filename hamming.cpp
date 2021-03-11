#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <time.h>
#include <fstream>



void program_options(int argc, char* argv[], int& pattern_length, std::string& text_file_name, int& text_length){
    namespace po = boost::program_options;
    po::options_description alldesc("Allowed options");
    alldesc.add_options()
                    ("help", "produce help message")
                    ("pattern_length", po::value<int>(&pattern_length)->default_value(9), "the length of pattern")
                    ("text_file_name", po::value<std::string>(&text_file_name)->default_value("data/Escherichia_coli_strain_FORC_028.fasta"), "file storing the input sequence")
                    ("text_length", po::value<int>(&text_length)->default_value(25), "the length of text");
    po::positional_options_description pos;
    pos.add("pattern_length", 1);
    pos.add("text_file_name", 1);
    pos.add("text_length", 1);

    po::options_description all;
    all.add(alldesc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc,argv,all),vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << alldesc << "\n";
    }

    if (vm.count("pattern_length")) {
        std::cout << "pattern_length was set to " 
        << vm["pattern_length"].as<int>() << ".\n";
    } else {
        std::cout << "pattern_length was not set.\n";
    }

    if (vm.count("text_file_name")) {
        std::cout << "text_file_name was set to " 
        << vm["text_file_name"].as<std::string >() << ".\n";
    } else {
        std::cout << "text_file_name was not set.\n";
    }

    if (vm.count("text_length")) {
        std::cout << "text_length was set to " 
        << vm["text_length"].as<int>() << ".\n";
    } else {
        std::cout << "pattern_length was not set.\n";
    }
}

int hamming_distance(std::string str1, std::string str2){
    int hamming_distance = 0;
    for(int i = 0; i < str1.length(); i++){
        if (str1.substr(i, 1).compare(str2.substr(i,1)) )
            hamming_distance += 1;
    }
    return hamming_distance;
}

int count_hamming_one(std::string text, std::string pattern){
    int counter = 0;
    int rest = 0;
    int pattern_length = pattern.length();
    for(int i = 0; i < text.length() - pattern_length + 1; i++){
        if(hamming_distance(text.substr(i, pattern_length), pattern) == 0){
            counter += 1;
        }
    }
    return counter;
}

std::string read_file(std::string file_name, unsigned int max_length){
    std::ifstream file(file_name);
    std::string str, result;
    unsigned int current_length = 0; 
    // avoid header line in the file
    std::getline(file, str);
    while(std::getline(file, str)){
        if(current_length + str.length() > max_length){
            result += str.substr(0, max_length - current_length);
            break;
        }
        current_length += str.length();
        result += str;
    }
    return result;
}

int main(int argc, char** argv){
    int pattern_length = 0;
    std::string text_file_name = "";
    int text_length = 0;
    program_options(argc, argv, pattern_length, text_file_name, text_length);
    

    std::string text_letter = read_file(text_file_name, text_length);
    std::string pattern_letter = read_file(text_file_name, pattern_length);

    auto duration_total = 0;
    for(int i = 0; i < 10; i++){
        auto start = std::chrono::high_resolution_clock::now();
        count_hamming_one(text_letter, pattern_letter);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( end - start ).count();
        duration_total += duration;
    }
    std::cout << "hamming_1:   " << count_hamming_one(text_letter, pattern_letter) << std::endl;

    std::cout << std::setprecision(9) << "milliseconds:  " << duration_total/10 << std::endl;


    return 0;
}