#include "ewsb.hpp"
#include "execute.hpp"
#include "results.hpp"
#include "settings.hpp"
#include "sugra.hpp"
#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(softsusy, m) { // NOLINT

  py::module_::import("pymicromegas");

  static const char *BASE_DOC_STRING_SUGRA = R"pbdoc(
Run micromegas with the parameters define at the GUT scale and settings using softsusy as an RGE backend.

Parameters
----------
params: Union[SugraParameters, List[SugraParameters]]
  Parameter object containing the model parameters defined at the GUT scale.
settings: MicromegasSettings, optional
  Settings object containing parameters to run micromegas with. Default is
  MicromegasSettings();

Returns
-------
results: MicromegasResults
  Results object containing requested results.
)pbdoc";

  static const char *BASE_DOC_STRING_EWSB = R"pbdoc(
Run micromegas with the parameters define at the EW scale and settings using softsusy as an RGE backend.

Parameters
----------
params: Union[EwsbParameters, List[EwsbParameters]]
  Parameter object containing the model parameters defined at the EW scale.
settings: MicromegasSettings, optional
  Settings object containing parameters to run micromegas with. Default is
  MicromegasSettings();

Returns
-------
results: MicromegasResults
  Results object containing requested results.
)pbdoc";

  m.doc() = "Python interface to micromegas with suspect as an RGE backend.";

  m.def(
      "softsusy",
      [](const SugraParameters &params,
         const MicromegasSettings &settings = MicromegasSettings()) {
        return execute(settings, params);
      },
      BASE_DOC_STRING_SUGRA, py::arg("params"),
      py::arg("settings") = MicromegasSettings(),
      py::return_value_policy::move);

  m.def(
      "softsusy",
      [](const std::vector<SugraParameters> &params,
         const MicromegasSettings &settings = MicromegasSettings()) {
        return execute(settings, params);
      },
      BASE_DOC_STRING_SUGRA, py::arg("params"),
      py::arg("settings") = MicromegasSettings(),
      py::return_value_policy::move);

  m.def(
      "softsusy",
      [](const EwsbParameters &params,
         const MicromegasSettings &settings = MicromegasSettings()) {
        return execute(settings, params);
      },
      BASE_DOC_STRING_SUGRA, py::arg("params"),
      py::arg("settings") = MicromegasSettings(),
      py::return_value_policy::move);

  m.def(
      "softsusy",
      [](const std::vector<EwsbParameters> &params,
         const MicromegasSettings &settings = MicromegasSettings()) {
        return execute(settings, params);
      },
      BASE_DOC_STRING_SUGRA, py::arg("params"),
      py::arg("settings") = MicromegasSettings(),
      py::return_value_policy::move);
}
