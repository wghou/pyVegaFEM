cmake_minimum_required(VERSION 3.4...3.18)
project(pyVegaFEM)

INCLUDE_DIRECTORIES(
    ${pyVegaFEM_SOURCE_DIR}
    ${pyVegaFEM_SOURCE_DIR}/VegaFEM
    ${pyVegaFEM_SOURCE_DIR}/python
)

# hrpd dir
set(VegaFEM_SOURCE_DIR 
    ${pyVegaFEM_SOURCE_DIR}/VegaFEM 
    ${pyVegaFEM_SOURCE_DIR}/VegaFEM/minivector 
    ${pyVegaFEM_SOURCE_DIR}/VegaFEM/ElasticForceFEM
    )

add_subdirectory(pybind11)

subdirs(VegaFEM)
subdirs(python)