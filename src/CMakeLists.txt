cmake_minimum_required(VERSION 3.20)

add_subdirectory(${CMAKE_SOURCE_DIR}/extern/pybind11 ${CMAKE_BINARY_DIR}/pybind11-build)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP)

pybind11_add_module(sqd_hpc_py
    cpp/module.cpp
    cpp/sqd_sbd_binding.cpp
)

target_include_directories(sqd_hpc_py PRIVATE
    ${CMAKE_SOURCE_DIR}/extern/qiskit-addon-sqd-hpc/include
    ${CMAKE_SOURCE_DIR}/extern/sbd/include
    ${CMAKE_SOURCE_DIR}/extern/boost    
)

target_link_libraries(sqd_hpc_py PRIVATE
    MPI::MPI_CXX
    BLAS::BLAS
    LAPACK::LAPACK
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(sqd_hpc_py PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS sqd_hpc_py
        LIBRARY DESTINATION sqd_hpc_py
)
