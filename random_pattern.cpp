#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <vector>

unsigned int prime_base = 0;

const unsigned int base = 5;

unsigned int * block_length_array;

unsigned int * reversed_block_length_array;

// how many pieces of the block divided into.
unsigned int number_blocks;

unsigned int * pattern_fingerprint_array;
std::string * pattern_number_array;

unsigned int * reversal_pattern_fingerprint_array;
std::string * reversal_pattern_number_array;

unsigned int * text_fingerprint_array;
std::string * text_number_array;

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

//generate all prime numbers
bool is_prime(int n) 
{ 
    // Corner case 
    if (n <= 1) 
        return false; 
  
    // Check from 2 to n-1 
    for (int i = 2; i < n; i++) 
        if (n % i == 0) 
            return false; 
  
    return true; 
}

void prime_base_generation(unsigned int n){
    std::vector<unsigned int> prime_vector;
    for (int i = 2; i <= n; i++){
        if(is_prime(i))
            prime_vector.push_back(i);
    }
    srand ( time(NULL) );
    prime_base = prime_vector.at(rand() % prime_vector.size());
}

void block_length_array_generation(std::string pattern, unsigned int block_length){
    block_length_array = new unsigned int[number_blocks];
    reversed_block_length_array = new unsigned int[number_blocks];

    for(int i = 0; i < number_blocks; ++i){
        block_length_array[i] = std::min(block_length, (unsigned int)(pattern.size() - i * block_length));
        reversed_block_length_array[number_blocks - 1 - i] = block_length_array[i];
    }

}

unsigned int pow_mod(unsigned int base, unsigned int times){
    unsigned int result = 1;
    for(int i = 0; i < times; ++i){
        result = result * base;
        result = result % prime_base;
    }
    return result;
}

// compute the fingerprint value for each block text
unsigned int fingerprint_computation(std::string block){
    unsigned int result = 0;
    for(int i = 0 ; i < block.size(); ++i){
        result += std::stoi(block.substr(i,1)) * pow_mod(base, block.size() - i - 1);
        result = result % prime_base;
    }
    return result;
}

void fingerprint_array_computation(std::string factor, unsigned int * fingerprint_array_pointer, std::string * text_array_pointer){
    unsigned int start_index = 0;
    for(int i = 0; i < number_blocks; ++i){
        text_array_pointer[i] = factor.substr(start_index, block_length_array[i]);
        fingerprint_array_pointer[i] = fingerprint_computation(text_array_pointer[i]);
        start_index += block_length_array[i];
    }
}

void reverse_array_entry(std::string * original_array, std::string * reverse_array, unsigned int * number_array){
    // make a copy of the orginal array; and reverse the array
    for(unsigned int i = 0; i < number_blocks; ++i){
        reverse_array[i] = original_array[number_blocks - i - 1];
    }

    // reverse each entry of the array;
    for(unsigned int i = 0; i < number_blocks; ++i){
        std::reverse(reverse_array[i].begin(), reverse_array[i].end());
    }

    // compute the fingerprint for the ververse_array
    for(unsigned int  i = 0; i < number_blocks; ++i){
        number_array[i] = fingerprint_computation(reverse_array[i]);
    }
    
}

int find_the_first_difference(std::string s1, std::string s2){
    for(int i = 0; i < s1.length(); ++i){
        if(s1.substr(i,1) != s2.substr(i,1)){
            return i;
        }
    }
    return -1;
}

void identify_different_blocks(unsigned int * first, unsigned int * last, bool * exist_difference){
    * exist_difference = false;
    * first = 0;
    * last = 0;
    // find the first block index
    for(int i = 0; i < number_blocks; ++i){
        // if the fingerprint of two blocks are equal to each. we will check the real content.
        if(text_fingerprint_array[i] == pattern_fingerprint_array[i]){
            if(text_number_array[i] != pattern_number_array[i]){
                * first = i;
                * exist_difference = true;
                break;
            }  
        }
        else{
            * first = i;
            * exist_difference = true;
            break;
        }
    }
    // fnd the last block index
    for(int i = number_blocks - 1; i >= 0; --i){
        if(text_fingerprint_array[i] == pattern_fingerprint_array[i]){
            if(text_number_array[i] != pattern_number_array[i]){
                * last = i;
                break;
            }  
        }
        else{
            * last = i;
            break;
        }
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
    for(int i = block_length - 1; i >= 0; --i){
        if(text_factor.substr(i,1) != pattern_factor.substr(i,1)){
            return i;
        }
    }
    return 0;
}

void retrieve_array(unsigned int * fingerprint_array, std::string * number_array, 
                    unsigned int * fingerprint_result, std::string * number_result,
                    unsigned int first_block_index, unsigned int first_offsite,
                    unsigned int last_block_index, unsigned int last_offsite, unsigned int * length_array){
    unsigned int index = 0;

    while(index <= last_block_index - first_block_index){
        if(index == 0){
            number_result[index] = number_array[first_block_index + index].substr(first_offsite, length_array[first_block_index + index] - first_offsite);
            fingerprint_result[index] = fingerprint_computation(number_result[index]);
            index ++;
            continue;
        }
        if(index == last_block_index - first_block_index){
            number_result[index] = number_array[first_block_index + index].substr(0, last_offsite + 1);
            fingerprint_result[index] = fingerprint_computation(number_result[index]);
            index ++;
            continue;
        }
        number_result[index] = number_array[first_block_index + index];
        fingerprint_result[index] = fingerprint_array[first_block_index + index];
        index ++;
    }
}

bool reversal_exists(unsigned int first_block_index, unsigned int first_offsite, 
                    unsigned int last_block_index, unsigned int last_offsite){
    
    //based on the index and offsite to retrieve the part of the text_finger_array.
    unsigned int temp_text_fingerprint_array[last_block_index - first_block_index + 1];
    std::string temp_text_number_array[last_block_index - first_block_index + 1];
    retrieve_array(text_fingerprint_array, text_number_array, temp_text_fingerprint_array, temp_text_number_array, 
                    first_block_index, first_offsite, last_block_index, last_offsite, block_length_array);
    
    

    unsigned int temp_reversal_pattern_fingerprint_array[last_block_index - first_block_index + 1];
    std::string temp_reversal_pattern_number_array[last_block_index - first_block_index + 1];
    unsigned int reversed_first_offsite = block_length_array[last_block_index] - last_offsite - 1;
    unsigned int reversed_last_offsite = block_length_array[first_block_index] - first_offsite - 1;
    retrieve_array(reversal_pattern_fingerprint_array, reversal_pattern_number_array, temp_reversal_pattern_fingerprint_array,
                    temp_reversal_pattern_number_array, number_blocks - last_block_index - 1, reversed_first_offsite,
                    number_blocks - first_block_index - 1, reversed_last_offsite, reversed_block_length_array);

    // off site
    unsigned int different_offsite = 0;
    int index = 0;
    unsigned int temp_fingerprint_store = 0;
    std::string temp_number_store;
    unsigned int move_fingerprint = 0;
    std::string move_number;


    //move temp text array forward
    if(first_offsite > reversed_first_offsite){
        temp_fingerprint_store = 0;
        temp_number_store = "";
        move_fingerprint = 0;
        temp_number_store = "";
        different_offsite = first_offsite - reversed_first_offsite;
        index = last_block_index - first_block_index;

        while(index >= 0){
            if(index == last_block_index - first_block_index){
                move_number = temp_text_number_array[index].substr(0, different_offsite);
                move_fingerprint = fingerprint_computation(move_number);
                
                temp_text_number_array[index] = temp_text_number_array[index].substr(different_offsite, temp_text_number_array[index].length() - different_offsite);
                
                temp_text_fingerprint_array[index] = temp_text_fingerprint_array[index] + prime_base - move_fingerprint * pow_mod(base, last_offsite + 1 - different_offsite);
                temp_text_fingerprint_array[index] = temp_text_fingerprint_array[index] % prime_base;

                index -- ;
                continue;
            }
            if(index == 0){
                temp_text_number_array[index] = temp_text_number_array[index] + move_number;

                temp_text_fingerprint_array[index] = (temp_text_fingerprint_array[index] * pow_mod(base, different_offsite)) % prime_base + move_fingerprint;
                temp_text_fingerprint_array[index] = temp_text_fingerprint_array[index] % prime_base;

                index --;
                continue;
            }
            temp_number_store = temp_text_number_array[index];
            temp_fingerprint_store = temp_text_fingerprint_array[index];

            temp_text_number_array[index] = temp_number_store.substr(different_offsite, temp_number_store.length() - different_offsite) + move_number;
            
            temp_text_fingerprint_array[index] = (temp_fingerprint_store + prime_base) - (fingerprint_computation(temp_number_store.substr(0, different_offsite)) * pow_mod(base, block_length_array[index + first_block_index] - different_offsite)) % prime_base;
            temp_text_fingerprint_array[index] = temp_text_fingerprint_array[index] * pow_mod(base, different_offsite) + move_fingerprint;
            temp_text_fingerprint_array[index] = temp_text_fingerprint_array[index] % prime_base;

            move_number = temp_number_store.substr(0, different_offsite);
            move_fingerprint = fingerprint_computation(move_number);
            index --;
        }
    }

    //move temp reversed pattern array forward
    else{
        temp_fingerprint_store = 0;
        temp_number_store = "";
        move_fingerprint = 0;
        temp_number_store = "";
        different_offsite = reversed_first_offsite - first_offsite;
        index = last_block_index - first_block_index;

        while(index >= 0){
            if(index == last_block_index - first_block_index){
                move_number = temp_reversal_pattern_number_array[index].substr(0, different_offsite);
                move_fingerprint = fingerprint_computation(move_number);
                
                temp_reversal_pattern_number_array[index] = temp_reversal_pattern_number_array[index].substr(different_offsite, temp_reversal_pattern_number_array[index].length() - different_offsite);
                
                temp_reversal_pattern_fingerprint_array[index] = temp_reversal_pattern_fingerprint_array[index] + prime_base - move_fingerprint * pow_mod(base, reversed_last_offsite + 1 - different_offsite);
                temp_reversal_pattern_fingerprint_array[index] = temp_reversal_pattern_fingerprint_array[index] % prime_base;

                index -- ;
                continue;
            }
            if(index == 0){
                temp_reversal_pattern_number_array[index] = temp_reversal_pattern_number_array[index] + move_number;

                temp_reversal_pattern_fingerprint_array[index] = (temp_reversal_pattern_fingerprint_array[index] * pow_mod(base, different_offsite)) % prime_base + move_fingerprint;
                temp_reversal_pattern_fingerprint_array[index] = temp_reversal_pattern_fingerprint_array[index] % prime_base;

                index --;
                continue;

            }
            temp_number_store = temp_reversal_pattern_number_array[index];
            temp_fingerprint_store = temp_reversal_pattern_fingerprint_array[index];

            temp_reversal_pattern_number_array[index] = temp_number_store.substr(different_offsite, temp_number_store.length() - different_offsite) + move_number;
            
            temp_reversal_pattern_fingerprint_array[index] = (temp_fingerprint_store + prime_base) - (fingerprint_computation(temp_number_store.substr(0, different_offsite)) * pow_mod(base, reversed_block_length_array[index + first_block_index] - different_offsite)) % prime_base;
            temp_reversal_pattern_fingerprint_array[index] = temp_reversal_pattern_fingerprint_array[index] * pow_mod(base, different_offsite) + move_fingerprint;
            temp_reversal_pattern_fingerprint_array[index] = temp_reversal_pattern_fingerprint_array[index] % prime_base;

            move_number = temp_number_store.substr(0, different_offsite);
            move_fingerprint = fingerprint_computation(move_number);
            index --;
        }
    }

    for(unsigned int i = 0; i <= last_block_index - first_block_index; ++i){
        if(temp_text_fingerprint_array[i] == temp_reversal_pattern_fingerprint_array[i]){
            if(temp_text_number_array[i] != temp_reversal_pattern_number_array[i]){
                return false;
            }
        }
        else{
            return false;
        }
    }
    return true;
}

//need to update text_fingerprint_array text_number_array;
void slide_window(std::string last_letter){
    unsigned int index = 0;
    std::string letter = "";
    std::string front = "";
    while(index < number_blocks - 1){
        letter = text_number_array[index + 1].substr(0,1);
        front = text_number_array[index].substr(0,1);
        text_number_array[index] = text_number_array[index].substr(1, block_length_array[index] - 1) + letter;

        text_fingerprint_array[index] = text_fingerprint_array[index] + prime_base;
        text_fingerprint_array[index] -= (std::stoi(front) * pow_mod(base, block_length_array[index] - 1)) % prime_base ;
        text_fingerprint_array[index] = (text_fingerprint_array[index] * base + std::stoi(letter)) % prime_base;

        index ++;
    }
    // last_block
    front = text_number_array[index].substr(0,1);
    text_number_array[index] = text_number_array[index].substr(1, block_length_array[index] - 1) + last_letter;

    text_fingerprint_array[index] = text_fingerprint_array[index] + prime_base;
    text_fingerprint_array[index] -= (std::stoi(front) * pow_mod(base, block_length_array[index] - 1)) % prime_base;
    text_fingerprint_array[index] = (text_fingerprint_array[index] * base + std::stoi(last_letter)) % prime_base;
}

unsigned int count_pattern(std::string text, std::string pattern){
    unsigned int block_length = sqrt(pattern.size());
    unsigned int count = 0;
    
    number_blocks = pattern.size()/block_length;
    pattern.size() % block_length > 0 ? number_blocks += 1: NULL;

    // initial and compute the each block's length;
    block_length_array_generation(pattern, block_length);
    
    // compute fingerprint array for pattern
    pattern_fingerprint_array = new unsigned int[number_blocks]; 
    pattern_number_array = new std::string[number_blocks];
    fingerprint_array_computation(pattern, pattern_fingerprint_array, pattern_number_array);

    //compute the reversal fingerprint array for pattern
    reversal_pattern_fingerprint_array = new unsigned int[number_blocks];
    reversal_pattern_number_array = new std::string[number_blocks];
    reverse_array_entry(pattern_number_array, reversal_pattern_number_array, reversal_pattern_fingerprint_array);

    text_fingerprint_array = new unsigned int[number_blocks];
    text_number_array = new std::string[number_blocks];

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
        identify_different_blocks(first, last, exist_pointer);
        if(exist_difference){
            // find two corresponding offsite in the blocks.
            first_offsite = find_first_index_linear(text_number_array[first_block_index], 
                            pattern_number_array[first_block_index], block_length_array[first_block_index]);
            last_offsite = find_last_index_linear(text_number_array[last_block_index], 
                            pattern_number_array[last_block_index], block_length_array[last_block_index]);
            if (reversal_exists(first_block_index, first_offsite, last_block_index, last_offsite))
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
        slide_window(text.substr(pattern.size() + index - 1, 1));
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


int main(int argc, char** argv){
    std::string text_letter = "ACGTAGCTTGCACAGT";
    std::string pattern_letter = "ACGT";

    std::string text_number = letter_to_number(text_letter);
    std::string pattern_number = letter_to_number(pattern_letter);
    
    std::string file_name = "data/Escherichia_coli_strain_FORC_028.fasta";
    std::string input_string = read_file(file_name, 400);
    //std::cout << input_string.length() << std::endl;

    // randomly generate the prime number based on the input.
    prime_base_generation(1000);



    unsigned int count = count_pattern(text_number, pattern_number);
    std::cout << count << std::endl;


    return 0;
}