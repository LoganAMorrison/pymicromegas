#include "micromegas.hpp"
#include "micromegas_interface.hpp"
#include <fmt/format.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void define_micromegas_get_set(py::module_ *ll) {

  using micromegas::Micromegas;
  using micromegas::MssmParameter;

  // ========================================================================
  // ---- CDM1 --------------------------------------------------------------
  // ========================================================================

  ll->def(
      "cdm1", []() { return Micromegas::get_cdm1(); },
      "Name of the 1st DM particle");

  // ========================================================================
  // ---- CDM2 --------------------------------------------------------------
  // ========================================================================

  ll->def(
      "cdm2", []() { return Micromegas::get_cdm2(); },
      "Name of the 2nd DM particle");

  // ========================================================================
  // ---- MCDM1 -------------------------------------------------------------
  // ========================================================================

  ll->def(
      "mcdm1", []() { return Micromegas::get_mcdm1(); },
      "Mass of the 1st DM particle");

  // ========================================================================
  // ---- MCDM2 -------------------------------------------------------------
  // ========================================================================

  ll->def(
      "mcdm2", []() { return Micromegas::get_mcdm2(); },
      "Mass of the 2nd DM particle");

  // ========================================================================
  // ---- VVDecay -----------------------------------------------------------
  // ========================================================================

  static const char *set_v_decay_doc = R"pbdoc(
  Set the value of VZdecay and VWdecay.

  Parameters
  ----------
  flag: int
    If 0, off-shell Z/W decays are turned off. If 1, 3-body
    final-states are included. If 2, Z's and W's are included in
    coannihilations as well->
  )pbdoc";

  ll->def(
      "set_v_decay", [](int flag) { return Micromegas::set_v_decay(flag); },
      set_v_decay_doc, py::arg("flag"));

  // ========================================================================
  // ---- assignVal ---------------------------------------------------------
  // ========================================================================

  static const char *assign_val_doc = R"pbdoc(
  Assign value val to parameter name. Returns a non-zero value if it cannot
  recognize a parameter name.

  Parameters
  ----------
  name: str
    Name of the parameter to modify.
  val: float
    New value of parameter.

  Returns
  -------
  err: int
    If non-zero, then parameter was not recognized. 
  )pbdoc";

  ll->def(
      "assign_val",
      [](const std::string &name, double val) {
        return Micromegas::assign_val(name, val);
      },
      assign_val_doc, py::arg("name"), py::arg("val"));

  // ========================================================================
  // ---- assignValW --------------------------------------------------------
  // ========================================================================
  static const char *assign_val_w_doc = R"pbdoc(
  Assign value val to parameter name. Returns a non-zero value and writes
  error message if it cannot recognize a parameter name.

  Parameters
  ----------
  name: str
    Name of the parameter to modify.
  val: float
    New value of parameter.

  Returns
  -------
  err: int
    If non-zero, then parameter was not recognized. 
  )pbdoc";

  ll->def(
      "assign_val_w",
      [](const std::string &name, double val) {
        return Micromegas::assign_val_w(name, val);
      },
      assign_val_w_doc, py::arg("name"), py::arg("val"));

  // ========================================================================
  // ---- assignVal (MSSM) --------------------------------------------------
  // ========================================================================

  static const char *mssm_assign_val_doc = R"pbdoc(
  Assign value val to parameter name. Returns a non-zero value if it cannot
  recognize a parameter name.

  Parameters
  ----------
  param: MssmParameter
    Enum value of the parameter.
  val: float
    New value of parameter.

  Returns
  -------
  err: int
    If non-zero, then parameter was not recognized. 
  )pbdoc";

  ll->def(
      "mssm_assign_val",
      [](MssmParameter param, double val) {
        return Micromegas::mssm_assign_val(param, val);
      },
      mssm_assign_val_doc, py::arg("param"), py::arg("val"));

  // ========================================================================
  // ---- assignValW (MSSM) -------------------------------------------------
  // ========================================================================

  static const char *mssm_assign_val_w_doc = R"pbdoc(
  Assign value val to parameter name. Returns a non-zero value and writes
  error message if it cannot recognize a parameter name.

  Parameters
  ----------
  param: MssmParameter
    Enum value of the parameter.
  val: float
    New value of parameter.

  Returns
  -------
  err: int
    If non-zero, then parameter was not recognized. 
  )pbdoc";

  ll->def(
      "mssm_assign_val_w",
      [](MssmParameter param, double val) {
        return Micromegas::mssm_assign_val_w(param, val);
      },
      mssm_assign_val_w_doc, py::arg("param"), py::arg("val"));

  // ========================================================================
  // ---- findVal -----------------------------------------------------------
  // ========================================================================

  static const char *find_val_doc = R"pbdoc(
  Finds the value of variable name and assigns it to parameter val. Returns
  NaN if parameter was not recognized.

  Parameters
  ----------
  name: str
    Name of the parameter to query.

  Returns
  -------
  val: int
    Value of the parameter. NaN if parameter was not recognized.
  )pbdoc";

  ll->def(
      "find_val",
      [](const std::string &name) { return Micromegas::find_val(name); },
      find_val_doc, py::arg("name"));

  // ========================================================================
  // ---- findValW ----------------------------------------------------------
  // ========================================================================

  static const char *find_val_w_doc = R"pbdoc(
  Finds the value of variable name and assigns it to parameter val. Returns
  NaN and print error message if parameter was not recognized.

  Parameters
  ----------
  name: str
    Name of the parameter to query.

  Returns
  -------
  val: int
    Value of the parameter. NaN if parameter was not recognized.
  )pbdoc";

  ll->def(
      "find_val_w",
      [](const std::string &name) { return Micromegas::find_val_w(name); },
      find_val_w_doc, py::arg("name"));

  // ========================================================================
  // ---- findVal (MSSM) ----------------------------------------------------
  // ========================================================================

  static const char *mssm_find_val_doc = R"pbdoc(
  Finds the value of variable name and assigns it to parameter val. Returns
  NaN if parameter was not recognized.

  Parameters
  ----------
  param: MssmParameter
    Enum value of the parameter to query.

  Returns
  -------
  val: int
    Value of the parameter. NaN if parameter was not recognized.
  )pbdoc";

  ll->def(
      "mssm_find_val",
      [](MssmParameter param) { return Micromegas::mssm_find_val(param); },
      mssm_find_val_doc, py::arg("name"));

  // ========================================================================
  // ---- findValW (MSSM) ---------------------------------------------------
  // ========================================================================

  static const char *mssm_find_val_w_doc = R"pbdoc(
  Finds the value of variable name and assigns it to parameter val. Returns
  NaN and prints error message if parameter was not recognized.

  Parameters
  ----------
  param: MssmParameter
    Enum value of the parameter to query.

  Returns
  -------
  val: int
    Value of the parameter. NaN if parameter was not recognized.
  )pbdoc";

  ll->def(
      "mssm_find_val_w",
      [](MssmParameter param) { return Micromegas::mssm_find_val_w(param); },
      mssm_find_val_w_doc, py::arg("param"));

  // ========================================================================
  // ---- pNum --------------------------------------------------------------
  // ========================================================================

  ll->def(
      "name_to_pdg",
      [](const std::string &name) { return Micromegas::name_to_pdg(name); },
      "Returns the PDG code of the particle with given name.");

  // ========================================================================
  // ---- pMass -------------------------------------------------------------
  // ========================================================================

  ll->def(
      "name_to_mass",
      [](const std::string &name) { return Micromegas::name_to_mass(name); },
      "Returns the mass of the particle with given name.");

  // ========================================================================
  // ---- pdg2Name ----------------------------------------------------------
  // ========================================================================

  ll->def(
      "pdg_to_name", [](int pdg) { return Micromegas::pdg_to_name(pdg); },
      "Returns the mass of the particle with given PDG code.");
}
