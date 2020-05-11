cmake_minimum_required(VERSION 3.15)
project(ISDBSCAN_R_master)

include_directories(/usr/local/lib/R/site-library/RcppArmadillo/include)
include_directories(/usr/local/lib/R/site-library/beachmat/include)
include_directories(/usr/local/lib/R/include /usr/local/lib/R/site-library/Rcpp/include)
include_directories(flann/src/cpp)

link_directories(flann/build/lib)

set(CMAKE_CXX_STANDARD 14)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")

# ------------------------------------------------------------------------------
# Set the compiler options (extra).
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic-errors")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")



add_executable(ISDBSCAN_Testing
        src/isdbscan.cpp
        src/influence_space.cpp
        src/clustering.cpp
        src/exact_knn.cpp
        src/approx_knn.cpp
        src/exact_knn.cpp
        src/testing.cpp src/inflo.cpp src/inflo.h src/stratifier.cpp src/stratifier.h)

target_link_libraries(ISDBSCAN_Testing flann_cpp_s)
target_link_libraries(ISDBSCAN_Testing lz4)
