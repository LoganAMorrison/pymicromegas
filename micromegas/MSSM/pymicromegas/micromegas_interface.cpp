#include "micromegas_interface.hpp"
#include "micromegas.hpp"
#include <fmt/format.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(micromegas, ll) { // NOLINT
  using micromegas::DirectDetectionAmps;
  using micromegas::DirectDetectionResults;
  using DDExp = micromegas::DirectDetectionExperiment;
  using micromegas::Micromegas;
  using micromegas::MssmParameter;

  py::add_ostream_redirect(ll, "ostream_redirect");

  ll.doc() = "Module providing low-level access to the micromegas library.";

  // ==================================
  // ---- Define Top-Level Classes ----
  // ==================================
  // using pymicromegas::EwsbParameters;
  // using pymicromegas::SugraParameters;

  // auto ewsb = py::class_<EwsbParameters>(m, "EwsbParameters");
  // pymicromegas::define_ewsb_parameters(&ewsb);

  // auto sugra = py::class_<SugraParameters>(m, "SugraParameters");
  // pymicromegas::define_sugra_parameters(&sugra);

  // ====================================
  // ---- Define Low-Level Interface ----
  // ====================================

  define_micromegas_direct_detection(&ll);
  define_micromegas_get_set(&ll);
  define_micromegas_utils(&ll);
  define_micromegas_pheno(&ll);

  // clang-format off
  py::enum_<MssmParameter>(ll, "MssmParameter", py::arithmetic(),
                           "Enumeration of the MSSM parameters")
      .value("Mu", MssmParameter::Mu, "Higgs Superpotential term: ℒ ⊃ μ Hu Hd")
      .value("Mg1", MssmParameter::Mg1, "U(1)_Y gaugino (Bino) mass: ℒ ⊃ -1/2 M₁ B B")
      .value("Mg2", MssmParameter::Mg2, "SU(2) gaugino (Wino) mass: ℒ ⊃ -1/2 M₂ Wᵃ Wₐ")
      .value("Mg3", MssmParameter::Mg3, "SU(3) gaugino (Gluino) mass: ℒ ⊃ -1/2 M₃ Gᵃ Gₐ")
      .value("Ml1", MssmParameter::Ml1, "Mass of the 1st gen slepton doublet: ℒ ⊃ -ML₁ * ⟨L₁, L₁⟩")
      .value("Ml2", MssmParameter::Ml2, "Mass of the 2nd gen slepton doublet: ℒ ⊃ -ML₂ * ⟨L₂, L₂⟩")
      .value("Ml3", MssmParameter::Ml3, "Mass of the 3rd gen slepton doublet: ℒ ⊃ -ML₃ * ⟨L₃, L₃⟩")
      .value("Mr1", MssmParameter::Mr1, "Mass of the 1st gen RH slepton : ℒ ⊃ -Ml₁ * |lR₁|²")
      .value("Mr2", MssmParameter::Mr2, "Mass of the 2nd gen RH slepton : ℒ ⊃ -Ml₂ * |lR₂|²")
      .value("Mr3", MssmParameter::Mr3, "Mass of the 3rd gen RH slepton : ℒ ⊃ -Ml₃ * |lR₃|²")
      .value("Mq1", MssmParameter::Mq1, "Mass of the 1st gen squark doublet: ℒ ⊃ -MQ₁ * ⟨Q₁, Q₁⟩")
      .value("Mq2", MssmParameter::Mq2, "Mass of the 2nd gen squark doublet: ℒ ⊃ -MQ₂ * ⟨Q₂, Q₂⟩")
      .value("Mq3", MssmParameter::Mq3, "Mass of the 3rd gen squark doublet: ℒ ⊃ -MQ₃ * ⟨Q₃, Q₃⟩")
      .value("Mu1", MssmParameter::Mu1, "Mass of the 1st gen RH up-squark : ℒ ⊃ -Mu₁ * |uR₁|²")
      .value("Mu2", MssmParameter::Mu2, "Mass of the 2nd gen RH up-squark : ℒ ⊃ -Mu₂ * |uR₂|²")
      .value("Mu3", MssmParameter::Mu3, "Mass of the 3rd gen RH up-squark : ℒ ⊃ -Mu₃ * |uR₃|²")
      .value("Md1", MssmParameter::Md1, "Mass of the 1st gen RH down-squark : ℒ ⊃ -Md₁ * |dR₁|²")
      .value("Md2", MssmParameter::Md2, "Mass of the 2nd gen RH down-squark : ℒ ⊃ -Md₂ * |dR₂|²")
      .value("Md3", MssmParameter::Md3, "Mass of the 3rd gen RH down-squark : ℒ ⊃ -Md₃ * |dR₃|²")
      .value("Mh3", MssmParameter::Mh3, "Mass of the psuedo-scalar Higgs boson MA")
      .value("Tb", MssmParameter::Tb, "Tangent beta: tan(β) = v₂ / v₁")
      .value("At", MssmParameter::At, "Trilinear top coupling: ℒ ⊃ At Yt tR ⟨Hu, Q₃⟩")
      .value("Ab", MssmParameter::Ab, "Trilinear bottom coupling: ℒ ⊃ Ab Yb bR ⟨Hd, Q₃⟩")
      .value("Al", MssmParameter::Al, "Trilinear tau coupling: ℒ ⊃ Al Yl lR ⟨Hd, Q₃⟩")
      .value("Au", MssmParameter::Au, "Trilinear 1st/2nd gen up-quark couplings:ℒ ⊃ Au [ Yu uR ⟨Hu, Q₁⟩ + Yc cR ⟨Hu, Q₂⟩ ]")
      .value("Ad", MssmParameter::Ad, "Trilinear 1st/2nd gen down-quark couplings: ℒ ⊃ Ad [ Yd dR ⟨Hd, Q₁⟩ + Ys sR ⟨Hd, Q₂⟩ ]")
      .value("Alfsmz", MssmParameter::Alfsmz, "alpha_QCD(MZ)")
      .value("Mz", MssmParameter::Mz, "Z-boson mass.  PDG value   91.1876 +/- 0.0021")
      .value("Mw", MssmParameter::Mw, "W-boson mass.  PDG value   80.385  +/- 0.015")
      .value("Mtp", MssmParameter::Mtp, "top quark pole mass")
      .value("Mbmb", MssmParameter::Mbmb, "Mb(Mb) running mass")
      .value("Mcmc", MssmParameter::Mcmc, "Mc(Mc) running mass")
      .value("Gg", MssmParameter::Gg, "Drived by CalcHEP")
      .value("Q", MssmParameter::Q, "QCD scale for running quark masses")
      .value("Ee", MssmParameter::Ee, "PDG value 4*pi/EE^2=128 corresponds to EE=0.3134")
      .value("Ml", MssmParameter::Ml, "mass of tau-lepton")
      .value("Mq", MssmParameter::Mq, "mass for light quarks")
      .value("Am", MssmParameter::Am, "No idea what this represents")
      .export_values();
  // clang-format on

  // ========================================
  // ---- Setting parameters + Calculate ----
  // ========================================

  ll.def(
      "mssm_ewsb", []() { return Micromegas::mssm_ewsb(); },
      "Calculates the masses of Higgs and supersymmetric "
      "particles in the MSSM including one-loop corrections "
      "starting from weak scale input parameters.");

  ll.def(
      "mssm_sugra",
      [](double tb, double gMG1, double gMG2, double gMG3, double gAl,
         double gAt, double gAb, double sgn, double gMHu, double gMHd,
         double gMl1, double gMl3, double gMr1, double gMr3, double gMq1,
         double gMq3, double gMu1, double gMu3, double gMd1, double gMd3) {
        return Micromegas::mssm_sugra(tb, gMG1, gMG2, gMG3, gAl, gAt, gAb, sgn,
                                      gMHu, gMHd, gMl1, gMl3, gMr1, gMr3, gMq1,
                                      gMq3, gMu1, gMu3, gMd1, gMd3);
      },
      R"pbdoc(
  Derives the parameters at the electroweak symmetry breaking scale assuming
  that all input parameters except `tb` and `sgn` are defined at the GUT
  scale.

  Parameters
  ----------
  tb: float
    Tangent beta: tan(β) = v₂ / v₁
  MG1: float
    Soft supersymmetry breaking U(1)_Y gaugino (Bino) mass: ℒ ⊃ -1/2 M₁ B B
  MG2: float 
    Soft supersymmetry breaking SU(2) gaugino (Wino) mass: ℒ ⊃ -1/2 M₂ Wᵃ Wₐ
  MG3: float
    Soft supersymmetry breaking SU(3) gaugino (Gluino) mass: ℒ ⊃ -1/2 M₃ Gᵃ Gₐ
  Al: float  
    Soft supersymmetry breaking trilinear tau coupling: ℒ ⊃ Al Yl lR ⟨Hd, Q₃⟩
  At: float  
    Soft supersymmetry breaking trilinear top coupling: ℒ ⊃ At Yt tR ⟨Hu, Q₃⟩
  Ab: float  
    Soft supersymmetry breaking trilinear bottom coupling: ℒ ⊃ Ab Yb bR ⟨Hd, Q₃⟩
  sgn: float
    Sign of the Higgs Superpotential coupling μ: ℒ ⊃ μ Hu Hd
  MHu: float 
    Soft supersymmetry breaking mass of the up-type Higgs: ℒ ⊃ - MHu Hu^* Hu
  MHd: float 
    Soft supersymmetry breaking mass of the down-type Higgs: ℒ ⊃ - MHd Hd^* Hd
  Ml1: float 
    Soft supersymmetry breaking mass of the 1st and 2nd gen slepton doublet:
    ℒ ⊃ -ML₁ * (⟨L₁, L₁⟩ + ⟨L₂, L₂⟩)
  Ml3: float
    Soft supersymmetry breaking mass of the 3rd gen slepton doublet:
    ℒ ⊃ -ML₃ * ⟨L₃, L₃⟩
  Mr1: float 
    Soft supersymmetry breaking mass of the 1st and 2nd gen RH slepton:
    ℒ ⊃ -Ml₁ * (|lR₁|² + |lR₂|²)
  Mr3: float 
    Soft supersymmetry breaking mass of the 3rd gen RH slepton:
    ℒ ⊃ -Ml₃ * |lR₃|²
  Mq1: float 
    Soft supersymmetry breaking mass of the 1st and 2nd gen squark 
    doublet: ℒ ⊃ -MQ₁ * (⟨Q₁, Q₁⟩ + ⟨Q₂, Q₂⟩)
  Mq3: float
    Soft supersymmetry breaking mass of the 3rd gen squark doublet:
    ℒ ⊃ -MQ₃ * ⟨Q₃, Q₃⟩
  Mu1: float 
    Soft supersymmetry breaking mass of the 1st and 2nd gen RH up-squark:
    ℒ ⊃ -Mu₁ * (|uR₁|² + |uR₂|²)
  Mu3: float 
    Soft supersymmetry breaking mass of the 3rd gen RH up-squark:
    ℒ ⊃ -Mu₃ * |uR₃|²
  Md1: float 
    Soft supersymmetry breaking mass of the 1st and 2nd gen RH down-squark:
    ℒ ⊃ -Md₁ * (|dR₁|² + |dR₂|²)
  Md3: float
    Soft supersymmetry breaking mass of the 3rd gen RH down-squark:
    ℒ ⊃ -Md₃ * |dR₃|²
    )pbdoc",
      py::arg("tb"), py::arg("MG1"), py::arg("MG2"), py::arg("MG3"),
      py::arg("Al"), py::arg("At"), py::arg("Ab"), py::arg("sgn"),
      py::arg("MHu"), py::arg("MHd"), py::arg("Ml1"), py::arg("Ml3"),
      py::arg("Mr1"), py::arg("Mr3"), py::arg("Mq1"), py::arg("Mq3"),
      py::arg("Mu1"), py::arg("Mu3"), py::arg("Md1"), py::arg("Md3"));
}
