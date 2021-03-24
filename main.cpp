#include <iostream>
#include <algorithm>  
#include <string>
#include <map>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <boost/program_options.hpp>
#include <time.h>

const unsigned int base = 5;

unsigned int * block_length_array;

unsigned int * reversed_block_length_array;

// how many pieces of the block divided into.
unsigned int number_blocks;

bool debug_mode = false;


std::map<std::string, std::string> create_letterMap(){
    std::map<std::string, std::string> letterMap;
    letterMap["A"] = "1";
    letterMap["C"] = "2";
    letterMap["G"] = "3";
    letterMap["T"] = "4";
    return letterMap;
}
const std::map<std::string, std::string> letterMap = create_letterMap();

// convert the letter string to number string
std::string letter_to_number(std::string content){
    std::string number_string = "";
    for(int i = 0; i < content.size(); ++i){
        number_string.append(letterMap.at(content.substr(i,1)));
    }
    return number_string;
}

// compute the fingerprint value for each block text
unsigned int fingerprint_computation(std::string block){
    unsigned int result = 0;
    for(int i = 0; i < block.size(); ++i){
        result += std::stoi(block.substr(i, 1)) * pow(base, block.size() -i - 1);
    }
    return result;
}

// compute the fingerprint array for text
void fingerprint_array_computation(std::string factor, unsigned int * fingerprint_array_pointer, std::string * text_array_pointer){
    unsigned int start_index = 0;
    for(int i = 0 ; i < number_blocks; ++i){
        text_array_pointer[i] = factor.substr(start_index, block_length_array[i]);
        fingerprint_array_pointer[i] = fingerprint_computation(text_array_pointer[i]);
        start_index += block_length_array[i];

    }
}

unsigned int find_first_index_linear(std::string text_factor, std::string pattern_factor, unsigned int block_length){
    for(int i = 0; i < block_length; ++i){
        if(text_factor.substr(i,1) != pattern_factor.substr(i,1)){
            return i;
        }
    }
    return 0;
}

unsigned int find_last_index_linear(std::string text_factor, std::string pattern_factor, unsigned int block_length){
    for(int i = block_length - 1; i>= 0 ; --i){
        if(text_factor.substr(i,1) != pattern_factor.substr(i,1)){
            return i;
        }
    }
    return 0;
}

// find the first different index in a block (binary search)
unsigned int find_first_index(unsigned int text_fingerprint, unsigned int pattern_fingerprint, unsigned int block_length){
    unsigned int start = 0, end = block_length - 1;
    unsigned int pivot = (start + end)/2;
    unsigned int text_temp = text_fingerprint, pattern_temp = pattern_fingerprint;
    unsigned int text_first, pattern_first;
    while(start != end ){
        text_first = text_temp / (unsigned int)(pow(base, end - pivot));
        pattern_first = pattern_temp/ (unsigned int)(pow(base, end - pivot));

        if(text_first == pattern_first){
            text_temp = text_temp - text_first * (unsigned int)pow(base, end - pivot);
            pattern_temp = pattern_temp - pattern_first * (unsigned int)pow(base, end - pivot);
            start = pivot + 1;

        }
        else{
            text_temp = text_first;
            pattern_temp = pattern_first;
            end = pivot;

        }
        pivot = (start + end)/2;
    }
    return start;
}

// find the last different index in a block (binary search)
unsigned int find_last_index(unsigned int text_fingerprint, unsigned int pattern_fingerprint, unsigned int block_length){
    unsigned int start = 0, end = block_length - 1;;
    unsigned int pivot = (start + end)/2;
    unsigned int text_temp = text_fingerprint, pattern_temp = pattern_fingerprint;
    unsigned int text_last, pattern_last;
    while(start != end){
        text_last = text_temp % (unsigned int)(pow(base, end - pivot));
        pattern_last = pattern_temp % (unsigned int)(pow(base, end - pivot));

        if(text_last == pattern_last){
            text_temp = (text_temp - text_last) / (unsigned int)(pow(base, end - pivot));
            pattern_temp = (pattern_temp - pattern_last) / (unsigned int)(pow(base, end - pivot));
            end = pivot;
        }
        else{
            text_temp = text_last;
            pattern_temp = pattern_last;
            start = pivot + 1;

        }
        pivot = (start + end)/2;
    }
    return start;
}

// identify first and last different blocks' index (linear search)
void identify_different_blocks(unsigned int * first, unsigned int * last, bool * exist_difference,
    unsigned int * text_fingerprint_array, unsigned int * pattern_fingerprint_array){
    * exist_difference = false;
    * first = 0;
    * last = 0;
    for(int i = 0; i < number_blocks; ++i){
        if(text_fingerprint_array[i] != pattern_fingerprint_array[i]){
            * first = i;
            * exist_difference = true;
            break;
        }
    }
    for(int i = number_blocks - 1; i >= 0; --i){
        if(text_fingerprint_array[i] != pattern_fingerprint_array[i]){
            * last = i;
            break;
        }
    }
}

void block_length_array_generation(std::string pattern, unsigned int block_length){
    block_length_array = new unsigned int[number_blocks];
    reversed_block_length_array = new unsigned int[number_blocks];

    for(int i = 0; i < number_blocks; ++i){
        block_length_array[i] = std::min(block_length, (unsigned int)(pattern.size() - i*block_length));
        reversed_block_length_array[number_blocks - 1 - i] = block_length_array[i];
    }
}

// determine whether is position is legal to switch
bool legal_position(unsigned int first_block, unsigned int first_off, unsigned int last_block, unsigned int last_off){
    if (first_block < last_block)
        return true;
    if(first_block == last_block && first_off < last_off)
        return true;
    return false;
}

//reverse the fingerprint block
unsigned int reverse_block(unsigned int block){
    unsigned int result = 0;
    while(block > 0){
        result = result * 5 + block % 5;
        block = block / 5;
    }
    return result;
}

//reverse the fingerprint array for each entry.
void reverse_array_entry(unsigned int * original, unsigned int * reversed){
    for(unsigned int i = 0; i < number_blocks; ++i){
        reversed[number_blocks - i - 1] = reverse_block(original[i]);
    }
}

void retrieve_array(unsigned int * text, unsigned int * result, unsigned int first_block_index, unsigned int first_offsite,
                    unsigned int last_block_index, unsigned int last_offsite, unsigned int * length_array){
    unsigned int index = 0;
    if(last_block_index == first_block_index){
        unsigned int remainder = text[last_block_index] % (unsigned int)pow(base, length_array[last_block_index] - first_offsite);
        result[0] = remainder / (unsigned int)pow(base, length_array[last_block_index] - last_offsite - 1);
    }
    else{
        while(index <= last_block_index - first_block_index){
            if(index == 0){
                result[index] =text[first_block_index + index]%(unsigned int)(pow(base, length_array[first_block_index + index] - first_offsite));
                index ++;
                continue;
            }
            if(index == last_block_index - first_block_index){

                result[index] =text[first_block_index + index]/(unsigned int)(pow(base, length_array[first_block_index + index] - last_offsite - 1));
                index ++;
                continue;
            }
            result[index] = text[first_block_index + index];
            index ++ ;
        }
    }
}

bool two_arrays_equal(unsigned int * arr1, unsigned int * arr2, unsigned int length){
    for(unsigned int i = 0; i < length; i++){
        if (arr1[i] != arr2[i])
            return false;
    }
    return true;
}
// based on the block_index and offsite to compare the middle part of text_fingerprint_array and
// the middle part of the reversal_pattern_fingerprint_array

bool reversal_exists(unsigned int * text, unsigned int * reversed_pattern, 
                     unsigned int first_block_index, unsigned int first_offsite,
                     unsigned int last_block_index, unsigned int last_offsite){
    // based on index and offsite to retrieve the part of text_fingerprint_array.
    unsigned int temp_text_fingerprint_array[last_block_index - first_block_index + 1];
    retrieve_array(text, temp_text_fingerprint_array, first_block_index, first_offsite, last_block_index, last_offsite, block_length_array);

    unsigned int temp_reversal_pattern_array[last_block_index - first_block_index + 1];
    unsigned int reversed_first_offsite = block_length_array[last_block_index] - last_offsite - 1;
    unsigned int reversed_last_offsite = block_length_array[first_block_index] - first_offsite - 1;
    retrieve_array(reversed_pattern, temp_reversal_pattern_array, number_blocks - last_block_index - 1, 
                reversed_first_offsite, number_blocks - first_block_index - 1,
                reversed_last_offsite, reversed_block_length_array);

    if(first_block_index == last_block_index){
        return two_arrays_equal(temp_text_fingerprint_array, temp_reversal_pattern_array,last_block_index - first_block_index + 1);    
    }
    // off 
    unsigned int different_offsite = 0;
    unsigned int move_number = 0;
    int index = 0;
    unsigned int border_block_leave_base = 0;
    unsigned int move_base = 0;
    unsigned int leave_base = 0;
    unsigned int temp_store = 0;
    // move text array forward
    if(first_offsite > reversed_first_offsite){
        temp_store = 0;
        different_offsite = first_offsite - reversed_first_offsite;
        index = last_block_index - first_block_index;
        border_block_leave_base = pow(base, last_offsite + 1 - different_offsite);
        move_base = pow(base, different_offsite);
        leave_base = pow(base, block_length_array[first_block_index + 1] - different_offsite);

        while(index >= 0 ){
            if(index == last_block_index - first_block_index){
                temp_store = temp_text_fingerprint_array[index];
                move_number = temp_store / border_block_leave_base;
                temp_text_fingerprint_array[index] = temp_store % border_block_leave_base;
                index --;
                continue;
            }
            if(index == 0){
                temp_text_fingerprint_array[index] = temp_text_fingerprint_array[index] * move_base + move_number;
                index --;
                continue;
            }
            temp_store = temp_text_fingerprint_array[index];
            temp_text_fingerprint_array[index] = temp_store % leave_base;
            temp_text_fingerprint_array[index] = temp_text_fingerprint_array[index] * move_base + move_number;
            move_number = temp_store / leave_base;
            index --;  
        }
    }
    // move reversed pattern array forward
    else{
        temp_store = 0;
        different_offsite = reversed_first_offsite - first_offsite;
        index = last_block_index - first_block_index;
        border_block_leave_base = pow(base,reversed_last_offsite + 1 - different_offsite);
        move_base = pow(base, different_offsite);
        leave_base = pow(base, block_length_array[first_block_index + 1] - different_offsite);
        while(index >= 0){
            if(index == last_block_index - first_block_index){
                temp_store = temp_reversal_pattern_array[index];
                // the front part which we want to move it to previous block. 
                // the length should be equal to different_offsite.
                move_number = temp_store / border_block_leave_base;
                temp_reversal_pattern_array[index] = temp_store % border_block_leave_base;
                index --;
                continue;
            }
            if(index  == 0){
                temp_reversal_pattern_array[index] = temp_reversal_pattern_array[index] * move_base + move_number;
                index --;
                continue;
            }
            temp_store = temp_reversal_pattern_array[index];
            temp_reversal_pattern_array[index] = temp_store % leave_base;
            temp_reversal_pattern_array[index] = temp_reversal_pattern_array[index] * move_base + move_number;
            move_number = temp_store / leave_base;
            index --;       
        }
    }
    return two_arrays_equal(temp_text_fingerprint_array, temp_reversal_pattern_array,last_block_index - first_block_index + 1);    
}

// slide the window
void slide_window(unsigned int * text_fingerprint_array, std::string * text_number_array, std::string last_letter){
    unsigned int index = 0;
    unsigned int front = text_fingerprint_array[index]/(unsigned int)(pow(base, block_length_array[index] - 1));
    unsigned int end = 0;
    while(index < number_blocks - 1){
        // remove front;
        text_fingerprint_array[index] -= front * (unsigned int)(pow(base, block_length_array[index] - 1));
        // shift the rest;
        text_fingerprint_array[index] = text_fingerprint_array[index] * base;
        end = text_fingerprint_array[index + 1]/(unsigned int)(pow(base, block_length_array[index + 1] - 1));
        
        text_fingerprint_array[index] += end;
        front = end;
        text_number_array[index] = text_number_array[index].substr(1, block_length_array[index] - 1) + 
                                   text_number_array[index + 1].substr(0,1);
        index++;
    }
    // last block
    text_fingerprint_array[index] -= front * (unsigned int)(pow(base, block_length_array[index] - 1));
    text_fingerprint_array[index] = text_fingerprint_array[index] * base;
    text_fingerprint_array[index] += std::stoi(last_letter);
    text_number_array[index] = text_number_array[index].substr(1, block_length_array[index] - 1) + last_letter;
}

std::string convert_array_to_string(std::string* array, int length){
    std::string result = "";
    for(int i = 0; i < length; ++i){
        result += array[i] + "  ";
    }
    return result;
}


unsigned int count_pattern(std::string text, std::string pattern){
    unsigned int block_length = sqrt(pattern.size());
    unsigned int count = 0;

    number_blocks = pattern.size()/block_length;
    pattern.size() % block_length > 0 ? number_blocks += 1: NULL;
    // initial and compute the each block's length;
    block_length_array_generation(pattern, block_length);
    
    // compute fingerprint array and pattern_number_array for pattern

    unsigned int pattern_fingerprint_array [number_blocks];
    std::string pattern_number_array[number_blocks]; 
    fingerprint_array_computation(pattern, pattern_fingerprint_array, pattern_number_array);

    //compute the reversal fingerprint array for pattern
    unsigned int reversal_pattern_fingerprint_array[number_blocks];
    reverse_array_entry(pattern_fingerprint_array, reversal_pattern_fingerprint_array);

    unsigned int text_fingerprint_array [number_blocks];
    std::string text_number_array [number_blocks];

    // the first and last blocks with the different fingerprint number
    unsigned int first_block_index = 0, last_block_index = 0;
    // first offsite index in first different block.
    // last offsite index in last different block.
    unsigned int first_offsite = 0, last_offsite = 0;

    bool exist_difference = false;
    unsigned int * first = &first_block_index, * last = &last_block_index;
    bool * exist_pointer = & exist_difference;
    unsigned int index = 0;
    // initial part.
    fingerprint_array_computation(text.substr(index, pattern.size()), text_fingerprint_array, text_number_array);
    // loop part.
    while(true){
        identify_different_blocks(first, last, exist_pointer, text_fingerprint_array, pattern_fingerprint_array);
        if(exist_difference){
            // find two corresponding offsite in the blocks.
            first_offsite = find_first_index_linear(text_number_array[first_block_index], 
                            pattern_number_array[first_block_index], block_length_array[first_block_index]);
            last_offsite = find_last_index_linear(text_number_array[last_block_index], 
                            pattern_number_array[last_block_index], block_length_array[last_block_index]);
            
            if (reversal_exists(text_fingerprint_array, reversal_pattern_fingerprint_array, 
                            first_block_index, first_offsite, last_block_index, last_offsite))
            {
                count += 1;
            }
        }
        else{
            count += 1;
        }
        index ++;
        if(index >= text.size() - pattern.size() + 1){
            break;
        }
        //slide window
        if(debug_mode){
            std::cout << "index:      " << index << std::endl;
            std::cout << "count:      " << count << std::endl;
            std::cout << first_block_index << "   " << first_offsite << "   " << last_block_index << "   " << last_offsite << std::endl;
            std::cout << "pattern:    " << convert_array_to_string(pattern_number_array, number_blocks) << std::endl;
            std::cout << "text:       " << convert_array_to_string(text_number_array, number_blocks) << std::endl;
            std::cout << "------------------------------" << std::endl;
        }
        
        slide_window(text_fingerprint_array, text_number_array, text.substr(pattern.size() + index - 1, 1));
    }
    first = NULL, last = NULL, exist_pointer = NULL;

    return count;
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

void program_options(int argc, char* argv[], int& pattern_length, std::string& text_file_name, int& text_length, bool& debug_mode){
    namespace po = boost::program_options;
    po::options_description alldesc("Allowed options");
    alldesc.add_options()
                    ("help", "produce help message")
                    ("pattern_length", po::value<int>(&pattern_length)->default_value(20), "the length of pattern")
                    ("text_file_name", po::value<std::string>(&text_file_name)->default_value("data/Escherichia_coli_strain_FORC_028.fasta"), "file storing the input sequence")
                    ("text_length", po::value<int>(&text_length)->default_value(25), "the length of text")
                    ("debug_mode", po::value<bool>(&debug_mode)->default_value(false), "print all factors");

    po::positional_options_description pos;
    pos.add("pattern_length", 1);
    pos.add("text_file_name", 1);
    pos.add("text_length", 1);
    pos.add("debug_mode", 1);

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

    if (vm.count("debug_mode")) {
        std::cout << "debug mode was set to " 
        << vm["debug_mode"].as<bool>() << ".\n";
    } else {
        std::cout << "debug_mode was not set.\n";
    }
}

int main(int argc, char** argv){
    /**
    int pattern_length;
    std::string text_file_name;
    int text_length;
    program_options(argc, argv, pattern_length, text_file_name, text_length, debug_mode);

    std::string text_letter = read_file(text_file_name, text_length);
    std::string pattern_letter = read_file(text_file_name, pattern_length);

    std::string text_number = letter_to_number(text_letter);
    std::string pattern_number = letter_to_number(pattern_letter);

    clock_t start, end;
    start = clock();
    unsigned int count = count_pattern(text_number, pattern_number);
    end = clock();

    std::cout << "final count:   " <<count << std::endl;
    float time_taken = end - start;
    std::cout << std::setprecision(9) << "milliseconds:  " << time_taken << std::endl;

    **/

    
}