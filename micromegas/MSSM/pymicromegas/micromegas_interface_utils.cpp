#include "micromegas.hpp"
#include "micromegas_interface.hpp"
#include <fmt/format.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void define_micromegas_utils(py::module_ *ll) {

  using micromegas::Micromegas;

  // ========================================================================
  // ---- CDM1 --------------------------------------------------------------
  // ========================================================================

  static const char *sort_odd_particles_doc = R"pbdoc(
  Sorts the odd particles with increasing masses.
 
  This routine fills the text parameters CDM1 and CDM2 with the names of the
  lightest odd particle starting with one and two tildes respectively and
  assigns the value of the mass of the lightest odd particle in each sector to
  the global parameters Mcdm1 and Mcdm2. For models with only one DM candidate,
  micrOMEGAs will set CDM2=NULL and Mcdm2=0. This routine returns a non zero
  error code for a wrong set of parameters, for example parameters for which
  some constraint cannot be calculated. The name of the corresponding
  constraint is written in txt. This routine has to be called after a
  reassignment of any input parameter.
)pbdoc";

  ll->def(
      "sort_odd_particles", []() { return Micromegas::sort_odd_particles(); },
      sort_odd_particles_doc);

  // ========================================================================
  // ---- CDM1 --------------------------------------------------------------
  // ========================================================================

  ll->def(
      "clean_decay_table", []() { return Micromegas::clean_decay_table(); },
      "Clears the decay table.");
}
