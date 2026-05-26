// bindings.cpp — pybind11 surface exposing the C++ core to Python.
//
// The exposed module is `bhrt_trisect_core` and is imported by the
// Python package in python/bhrt_trisect/__init__.py when available;
// when not built, the pure-Python implementation is used transparently.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bhrt_trisect.hpp"

namespace py = pybind11;

PYBIND11_MODULE(bhrt_trisect_core, m) {
    m.doc() = "C++/CUDA core for the BHRT trisection pipeline.";

    py::class_<bhrt::Triangulation>(m, "Triangulation")
        .def(py::init<>())
        .def("size", &bhrt::Triangulation::size)
        .def("n_vertices", &bhrt::Triangulation::nVertices)
        .def("n_edges", &bhrt::Triangulation::nEdges)
        .def("n_triangles", &bhrt::Triangulation::nTriangles)
        .def("n_tetrahedra", &bhrt::Triangulation::nTetrahedra)
        .def("euler", &bhrt::Triangulation::eulerCharacteristic)
        .def("is_closed", &bhrt::Triangulation::isClosed)
        .def("isosig", &bhrt::Triangulation::isoSig)
        .def("edge_degree_isosig", &bhrt::Triangulation::edgeDegreeIsoSig);

    py::class_<bhrt::TSColouring>(m, "TSColouring")
        .def_readonly("colour", &bhrt::TSColouring::colour)
        .def_readonly("is_c", &bhrt::TSColouring::is_c)
        .def_readonly("is_ts", &bhrt::TSColouring::is_ts)
        .def_readonly("audit", &bhrt::TSColouring::audit);

    m.def("is_ts_tricolouring", &bhrt::isTsTricolouring,
          py::arg("T"), py::arg("colour"), py::arg("collapse_budget") = 10000);
    m.def("enumerate_ts_tricolourings", &bhrt::enumerateTsTricolourings,
          py::arg("T"), py::arg("collapse_budget") = 10000,
          py::arg("stop_at_first") = false);

    py::class_<bhrt::CellularSurface>(m, "CellularSurface")
        .def_readonly("num_vertices", &bhrt::CellularSurface::numVertices)
        .def_readonly("edges", &bhrt::CellularSurface::edges)
        .def_readonly("faces", &bhrt::CellularSurface::faces)
        .def("genus", &bhrt::CellularSurface::genus);

    py::class_<bhrt::CutCurve>(m, "CutCurve")
        .def_readonly("edge_sequence", &bhrt::CutCurve::edge_sequence)
        .def_readonly("role", &bhrt::CutCurve::role);

    py::class_<bhrt::TrisectionDiagram>(m, "TrisectionDiagram")
        .def_readonly("surface", &bhrt::TrisectionDiagram::surface)
        .def_readonly("alpha", &bhrt::TrisectionDiagram::alpha)
        .def_readonly("beta", &bhrt::TrisectionDiagram::beta)
        .def_readonly("gamma", &bhrt::TrisectionDiagram::gamma)
        .def_readonly("audit", &bhrt::TrisectionDiagram::audit);

    m.def("extract_diagram", &bhrt::extractDiagram,
          py::arg("T"), py::arg("colouring"));

    py::class_<bhrt::SearchConfig>(m, "SearchConfig")
        .def(py::init<>())
        .def_readwrite("beam_width", &bhrt::SearchConfig::beam_width)
        .def_readwrite("excess_height", &bhrt::SearchConfig::excess_height)
        .def_readwrite("time_limit_seconds", &bhrt::SearchConfig::time_limit_seconds)
        .def_readwrite("max_iterations", &bhrt::SearchConfig::max_iterations)
        .def_readwrite("seed", &bhrt::SearchConfig::seed)
        .def_readwrite("use_gpu_scorer", &bhrt::SearchConfig::use_gpu_scorer)
        .def_readwrite("checkpoint_dir", &bhrt::SearchConfig::checkpoint_dir);

    m.def("beam_search", &bhrt::beamSearch,
          py::arg("start"), py::arg("config"));
}
