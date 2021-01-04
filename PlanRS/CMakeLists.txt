cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)

project(PlanRS VERSION 1.0)

add_compile_options(-std=c++11)

find_package(OMPL)
if (NOT OMPL_FOUND)
  message(SEND_ERROR "Failed to find OMPL! Install OMPL and try again!")
endif ()

include_directories(
    include
    ${OMPL_INCLUDE_DIRS})
    
include_directories(/include/ompl-1.5.1/)

# Add support in Boost::Python for std::shared_ptr
# This is a hack that replaces boost::shared_ptr related code with std::shared_ptr.
# Proper support for std::shared_ptr was added in Boost 1.63.
if (NOT Boost_VERSION_STRING)
    if (Boost_VERSION_MAJOR)
        set(Boost_VERSION_STRING
            "${Boost_VERSION_MAJOR}.${Boost_VERSION_MINOR}.${Boost_VERSION_PATCH}")
    else()
        set(Boost_VERSION_STRING
            "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_PATCH_VERSION}")
    endif()
endif()
if(Boost_VERSION_STRING VERSION_LESS "1.63.0")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/external")
endif()

add_library(ReedsSheppPlanner src/ReedsSheppPlanner.cpp)
target_link_libraries(ReedsSheppPlanner ${OMPL_LIBRARIES})

add_executable(main src/main.cpp)
target_link_libraries(main ReedsSheppPlanner)