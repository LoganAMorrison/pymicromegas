#include "pymicromegas.hpp"
#include <iomanip>
#include <sstream>

namespace py = pybind11;

static constexpr int DEFUALT_PRECISION = 5;

void define_sugra_parameters(pybind11::class_<SugraParameters> *sugra) {
  // =====================
  // ---- Constructor ----
  // =====================

  static const char *sugra_init_default_doc = R"pbdoc(
  Construct a Parameters object with parameters defined at the GUT scale.
  )pbdoc";

  sugra->def(py::init<>(), sugra_init_default_doc);

  static const char *sugra_init_full_doc = R"pbdoc(
  Construct a Parameters object with parameters defined at the GUT scale.

  Parameters
  ----------
  m0: float
    Universal scalar (sfermion and Higgs boson) masses: 
      m₀ = MQᵢ(Mgut) = MuRᵢ(Mgut) = MdRᵢ(Mgut) = MLᵢ(Mgut) = MlRᵢ(Mgut) 
        = MHu(Mgut) = MHd(Mgut).
  mhf: float
    Unification of gaugino (bino, wino and gluino) masses: 
      M₁(Mgut) = M₂(Mgut) = M₃(Mgut) = mhf.
  a0: float
    Universal trilinear couplings:
      Auᵤᵥ(Mgut) = Adᵤᵥ(Mgut) = Alᵤᵥ(Mgut) = A₀ δᵤᵥ.
  tb: float
    Tangent beta: tan(β) = v₂ / v₁ (defined at EW scale.)
  sgn: float
    Sign of the coefficient of the Higgs superpotential: μ Hu Hd 
    (defined at EW scale.)
  )pbdoc";

  sugra->def(py::init<double, double, double, double, double>(),
             sugra_init_full_doc, py::arg("m0"), py::arg("mhf"), py::arg("a0"),
             py::arg("tb"), py::arg("sgn"));

  // ==============
  // ---- Repr ----
  // ==============

  sugra->def("__repr__", [](const SugraParameters &p) {
    std::stringstream stream;
    stream << std::scientific << std::setprecision(DEFUALT_PRECISION);
    stream << "SugraParameters(";
    stream << "m0 = " << p.get_m0() << ", ";
    stream << "mhf = " << p.get_mhf() << ", ";
    stream << "a0 = " << p.get_a0() << ", ";
    stream << "sgn = " << p.get_sgn();
    stream << ")";
    return stream.str();
  });

  // ==========================
  // ---- Sugra Parameters ----
  // ==========================

  sugra->def_property("m0", &SugraParameters::get_m0, &SugraParameters::set_m0,
                      "Universal scalar (sfermion and Higgs boson) masses: "
                      "MQᵢ(Mgut) = MuRᵢ(Mgut) = MdRᵢ(Mgut) = MLᵢ(Mgut) = "
                      "MlRᵢ(Mgut) = MHu(Mgut) = MHd(Mgut) = m₀.");

  sugra->def_property("mhf", &SugraParameters::get_mhf,
                      &SugraParameters::set_mhf,
                      "Unification of gaugino (bino, wino and gluino) masses: "
                      "M₁(Mgut) = M₂(Mgut) = M₃(Mgut) = mhf.");

  sugra->def_property(
      "a0", &SugraParameters::get_mhf, &SugraParameters::set_mhf,
      "Universal trilinear couplings: Auᵤᵥ(Mgut) = Adᵤᵥ(Mgut) = "
      "Alᵤᵥ(Mgut) = A₀ δᵤᵥ.");

  sugra->def_property("tb", &SugraParameters::get_tb, &SugraParameters::set_tb,
                      "Tangent beta: tan(β) = v₂ / v₁");

  sugra->def_property(
      "sgn", &SugraParameters::get_sgn, &SugraParameters::set_sgn,
      "Sign of the coefficient of the Higgs superpotential μ Hu Hd");
}
