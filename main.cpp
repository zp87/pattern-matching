#include <iostream>
#include <algorithm>  
#include <string>
#include <map>
#include <math.h>

const unsigned int base = 5;

unsigned int * block_length_array;


// how many pieces of the block divided into.
unsigned int number_blocks;


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
void fingerprint_array_computation(std::string factor, unsigned int * fingerprint_array_pointer){
    unsigned int start_index = 0;
    for(int i = 0; i < number_blocks; ++i){
        fingerprint_array_pointer[i] = fingerprint_computation(
            factor.substr(start_index, block_length_array[i]));
            start_index += block_length_array[i];
    }
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

// identify first and last different positions
void identify_different_blocks(unsigned int * first, unsigned int * last, bool * exist_difference,
    unsigned int * text_fingerprint_array, unsigned int * pattern_fingerprint_array){
    * exist_difference = false;
    * first = 0;
    * last = 0;
    for(unsigned int i = 0; i < number_blocks; ++i){
        if(text_fingerprint_array[i] != pattern_fingerprint_array[i]){
            if(! *exist_difference){
                * first = i;
            }
            * last = i;
            * exist_difference = true;
        }
    }
}

void block_length_array_generation(std::string pattern, unsigned int block_length){
    block_length_array = new unsigned int[number_blocks];
    for(int i = 0; i < number_blocks; ++i){
        block_length_array[i] = std::min(block_length, (unsigned int)(pattern.size() - i*block_length));
    }
}

// check wheter two arrays contain same values in order.
bool array_same(unsigned int * arr1, unsigned int * arr2){
    for(unsigned int i = 0; i < number_blocks; ++i){
        if(arr1[i] != arr2[i])
            return false;
    }
    return true;
}

// determine whether is position is legal to switch
bool legal_position(unsigned int first_block, unsigned int first_off, unsigned int last_block, unsigned int last_off){
    if (first_block < last_block)
        return true;
    if(first_block == last_block && first_off < last_off)
        return true;
    return false;
}

// compute the number for the letter based on offsite and block length
unsigned int single_number(unsigned int fingerprint, unsigned int offsite, unsigned int block_length){
    fingerprint = fingerprint % (unsigned int)pow(base, block_length - offsite);
    return fingerprint / ((unsigned int)pow(base, block_length - offsite - 1));
}

// reverse the text_fingerprint_array based on the block index and offsite.
void reverse_fingerprint_array(unsigned int * reversal, 
            unsigned int first_block_index, unsigned int first_offsite,
            unsigned int last_block_index, unsigned int last_offsite){

    bool indicator = legal_position(first_block_index, first_offsite, last_block_index, last_offsite);
    unsigned int front = 0, end = 0;

    while(indicator){
        // compute the correspondings number for the letter.
        front = single_number(reversal[first_block_index], first_offsite, block_length_array[first_block_index]);
        end = single_number(reversal[last_block_index], last_offsite, block_length_array[last_block_index]);
    
        // swith the position
        reversal[first_block_index] -= front * (unsigned int)pow(base, block_length_array[first_block_index] - first_offsite - 1);
        reversal[first_block_index] += end * (unsigned int)pow(base, block_length_array[first_block_index] - first_offsite - 1);
        reversal[last_block_index] -= end * (unsigned int)pow(base, block_length_array[last_block_index] - last_offsite - 1);
        reversal[last_block_index] += front * (unsigned int)pow(base, block_length_array[last_block_index] - last_offsite - 1);

        // compute the new position
        first_offsite += 1;
        if(first_offsite >= block_length_array[first_block_index]){
            first_block_index += 1;
            first_offsite -= block_length_array[first_block_index];
        }
        if(last_offsite > 0){
            last_offsite -= 1;
        }
        else{
            last_offsite = block_length_array[first_block_index] - 1;
            last_block_index -= 1;
        }
        indicator = legal_position(first_block_index, first_offsite, last_block_index, last_offsite);
    }    

}

//check two arrays are equal.
bool two_arrays_equal(unsigned int * arr1, unsigned int * arr2){
    for(unsigned int i = 0; i < number_blocks; ++i){
        if(arr1[i] != arr2[i])
            return false;
    }
    return true;
}

// array copy arr1 -> arr2
void array_copy(unsigned int * arr1, unsigned int * arr2){
    for(unsigned int i = 0 ; i < number_blocks; ++i){
        arr2[i] = arr1[i];
    }
}

// slide the window
void slide_window(unsigned int * text_fingerprint_array, std::string last_letter){
    //std::cout<< last_letter << std::endl;
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
        index ++;
    }
    // last block
    text_fingerprint_array[index] -= front * (unsigned int)(pow(base, block_length_array[index] - 1));
    text_fingerprint_array[index] = text_fingerprint_array[index] * base;
    text_fingerprint_array[index] += std::stoi(last_letter);
}


unsigned int count_pattern(std::string text, std::string pattern){
    unsigned int block_length = sqrt(pattern.size());
    unsigned int count = 0;
    
    number_blocks = pattern.size()/block_length;
    pattern.size() % block_length > 0 ? number_blocks += 1: NULL;

    // initial and compute the each block's length;
    block_length_array_generation(pattern, block_length);
    
    // compute fingerprint array for pattern
    unsigned int pattern_fingerprint_array [number_blocks]; 
    fingerprint_array_computation(pattern, pattern_fingerprint_array);

    unsigned int text_fingerprint_array [number_blocks];
    unsigned int reversal_text_fingerprint_array[number_blocks];
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
    fingerprint_array_computation(text.substr(index, pattern.size()), text_fingerprint_array);
    array_copy(text_fingerprint_array, reversal_text_fingerprint_array);
    // loop part.
    while(true){
        identify_different_blocks(first, last, exist_pointer, text_fingerprint_array, pattern_fingerprint_array);
        if(exist_difference){
            // find two corresponding offsite in the blocks.
            first_offsite = find_first_index(text_fingerprint_array[first_block_index], 
                            pattern_fingerprint_array[first_block_index], block_length_array[first_block_index]);
            last_offsite = find_last_index(text_fingerprint_array[last_block_index], 
                            pattern_fingerprint_array[last_block_index], block_length_array[last_block_index]);
        
            // reversal the letters and update the text_fingerprint_array!
            reverse_fingerprint_array(reversal_text_fingerprint_array, first_block_index, 
                                first_offsite, last_block_index, last_offsite);

            // check the reversal_text_fingerprint_array is equal to pattern_fingerprint_array 
            if(two_arrays_equal(reversal_text_fingerprint_array, pattern_fingerprint_array))
                count += 1;
        }
        index ++;
        if(index >= text.size() - pattern.size() + 1){
            break;
        }
        //slide window
        slide_window(text_fingerprint_array, text.substr(pattern.size() + index - 1, 1));
    }

    first = NULL, last = NULL, exist_pointer = NULL;

    return count;
}



int main(){
    std::string text_letter = "ACTTGCTGACAACTGCAC";
    std::string pattern_letter = "ACTTGTCAACAGTCGCA";

    std::string text_number = letter_to_number(text_letter);
    std::string pattern_number = letter_to_number(pattern_letter);
    
    unsigned int count = count_pattern(text_number, pattern_number);
    std::cout << count << std::endl;
    return 0;
}