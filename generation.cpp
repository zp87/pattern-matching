#include <iostream>
#include <string>
#include <time.h>
#include <vector>


std::string pattern_generation(int pattern_length){
    std::string pattern = "";
    srand ( time(0) );
    for(int i = 0; i < pattern_length; i++){
        pattern = pattern + std::to_string(rand() % 4);
    }
    return pattern;
}
bool range_check(int index, std::vector<int> index_vector, int width){
    for(int i = 0; i < index_vector.size(); i++){
        if(index > index_vector.at(i) - width && index < index_vector.at(i) + width){
            return true;
        }
    }
    return false;
}

std::string reverse_pattern(std::string pattern){
    // random two numbers;
    srand ( time(0) );

    int start = rand() % pattern.size();
    int end = 0;
    while(true){
        end = rand() % pattern.size();
        if (end < start - 1 || end > start + 1)
            break;
    }
    if (start > end){
        int temp = start;
        start = end;
        end = temp;
    }
    std::string substring = pattern.substr(start, end - start);
    std::reverse(substring.begin(), substring.end());
    pattern.replace(start, end-start, substring);
    return pattern;
}

std::string generation(std::string* text, int required_total, int required_pattern_length){
    std::string pattern = pattern_generation(required_pattern_length);
    std::string reversed_pattern = "";
    int start_index = 0;
    std::vector<int> index_vector;
    srand ( time(0) );
    int number = 0;
    while(number < required_total){
        start_index = rand() % (text -> length() - required_pattern_length );
        if(!range_check(start_index, index_vector, required_pattern_length)){
            reversed_pattern = reverse_pattern(pattern);
            if(pattern.compare(reversed_pattern) == 0){
                continue;
            }
            index_vector.push_back(start_index);
            text -> replace(start_index, required_pattern_length, reversed_pattern);
            number += 1;
        }
    }
    return pattern;


}
