#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_sqd_sbd(py::module_ &m);

PYBIND11_MODULE(sqd_hpc_py, m) {
    m.doc() = "SQD-HPC + SBD python wrapper";
    init_sqd_sbd(m);
}
