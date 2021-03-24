cmake_minimum_required(VERSION 2.6)

project(pattern_matching)
find_package(Boost REQUIRED COMPONENTS program_options)

include_directories(SYSTEM ${BOOST_INCLUDE_DIR})
    
#add_executable(pattern_matching main.cpp )
#target_link_libraries(pattern_matching ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(pattern_matching_random random_pattern.cpp )
target_link_libraries(pattern_matching_random ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(pattern_matching_random_align random_pattern_align.cpp )
target_link_libraries(pattern_matching_random_align ${Boost_PROGRAM_OPTIONS_LIBRARY})


#add_executable(hamming hamming.cpp)
#target_link_libraries(hamming ${Boost_PROGRAM_OPTIONS_LIBRARY})

# add_executable(generation generation.cpp)
