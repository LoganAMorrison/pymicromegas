#include "pymicromegas.hpp"
#include "ewsb.hpp"
#include "results.hpp"
#include "settings.hpp"
#include "sugra.hpp"
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pymicromegas, m) { // NOLINT

  auto sugra = py::class_<SugraParameters>(m, "SugraParameters");
  auto ewsb = py::class_<EwsbParameters>(m, "EwsbParameters");
  auto results = py::class_<MicromegasResults>(m, "MicromegasResults");
  auto settings = py::class_<MicromegasSettings>(m, "MicromegasSettings");

  py::add_ostream_redirect(m, "ostream_redirect");

  define_sugra_parameters(&sugra);
  define_ewsb_parameters(&ewsb);
  define_results(&results);
  define_settings(&settings);
}
