#include "pymicromegas.hpp"
#include <iomanip>
#include <sstream>

namespace py = pybind11;

static constexpr int DEFUALT_PRECISION = 5;

void define_ewsb_parameters(pybind11::class_<EwsbParameters> *ewsb) {
  // =====================
  // ---- Constructor ----
  // =====================

  ewsb->def(
      py::init<>(),
      "Construct an empty Parameters object with parameters defined at the "
      "electroweak scale");

  ewsb->def(py::init<double, double, double, double, double, double, double,
                     double, double, double, double, double, double, double,
                     double, double, double, double, double, double, double,
                     double, double, double>(),
            "Construct a Parameters object with parameters defined at the "
            "electroweak scale.",
            py::arg("mu"), py::arg("mg1"), py::arg("mg2"), py::arg("mg3"),
            py::arg("ml1"), py::arg("ml2"), py::arg("ml3"), py::arg("mr1"),
            py::arg("mr2"), py::arg("mr3"), py::arg("mq1"), py::arg("mq2"),
            py::arg("mq3"), py::arg("mu1"), py::arg("mu2"), py::arg("mu3"),
            py::arg("md1"), py::arg("md2"), py::arg("md3"), py::arg("mh3"),
            py::arg("tb"), py::arg("at"), py::arg("ab"), py::arg("al"));

  // ==============
  // ---- Repr ----
  // ==============

  ewsb->def("__repr__", [](const EwsbParameters &p) {
    std::stringstream stream;
    stream << std::scientific << std::setprecision(DEFUALT_PRECISION);
    stream << "EwsbParameters(\n";

    stream << "\tmg1 = " << p.get_mg1() << ", "
           << "\tmg2 = " << p.get_mg2() << ", "
           << "\tmg3 = " << p.get_mg3() << ", "
           << "\n";

    stream << "\tml1 = " << p.get_ml1() << ", "
           << "\tml2 = " << p.get_ml2() << ", "
           << "\tml3 = " << p.get_ml3() << ", "
           << "\n";

    stream << "\tmr1 = " << p.get_mr1() << ", "
           << "\tmr2 = " << p.get_mr2() << ", "
           << "\tmr3 = " << p.get_mr3() << ", "
           << "\n";

    stream << "\tmq1 = " << p.get_mq1() << ", "
           << "\tmq2 = " << p.get_mq2() << ", "
           << "\tmq3 = " << p.get_mq3() << ", "
           << "\n";

    stream << "\tmu1 = " << p.get_mu1() << ", "
           << "\tmu2 = " << p.get_mu2() << ", "
           << "\tmu3 = " << p.get_mu3() << ", "
           << "\n";

    stream << "\tmd1 = " << p.get_md1() << ", "
           << "\tmd2 = " << p.get_md2() << ", "
           << "\tmd3 = " << p.get_md3() << ", "
           << "\n";

    stream << "\tmh3 = " << p.get_mh3() << ", "
           << "\tmu = " << p.get_mu() << ", "
           << "\ttb = " << p.get_tb() << ", "
           << "\n";

    stream << "\tat = " << p.get_at() << ", "
           << "\tab = " << p.get_ab() << ", "
           << "\tal = " << p.get_al() << ", "
           << "\n";
    stream << ")";

    return stream.str();
  });

  // ========================
  // ---- Gaugino Masses ----
  // ========================

  ewsb->def_property("mg1", &EwsbParameters::get_mg1, &EwsbParameters::set_mg1,
                     "U(1)_Y gaugino (Bino) mass: ℒ ⊃ -1/2 M₁ B B");

  ewsb->def_property("mg2", &EwsbParameters::get_mg2, &EwsbParameters::set_mg2,
                     "SU(2) gaugino (Wino) mass: ℒ ⊃ -1/2 M₂ Wᵃ Wₐ.");

  ewsb->def_property("mg3", &EwsbParameters::get_mg3, &EwsbParameters::set_mg3,
                     "SU(3) gaugino (Gluino) mass: ℒ ⊃ -1/2 M₃ Gᵃ Gₐ.");

  // ============================
  // ---- LH Slepton Doublet ----
  // ============================

  ewsb->def_property(
      "ml1", &EwsbParameters::get_ml1, &EwsbParameters::set_ml1,
      "Mass of the 1st gen slepton doublet: ℒ ⊃ -ML₁ * ⟨L₁, L₁⟩");
  ewsb->def_property(
      "ml2", &EwsbParameters::get_ml2, &EwsbParameters::set_ml2,
      "Mass of the 2nd gen slepton doublet: ℒ ⊃ -ML₂ * ⟨L₂, L₂⟩");
  ewsb->def_property(
      "ml3", &EwsbParameters::get_ml3, &EwsbParameters::set_ml3,
      "Mass of the 3rd gen slepton doublet: ℒ ⊃ -ML₃ * ⟨L₃, L₃⟩");

  // ====================
  // ---- RH Slepton ----
  // ====================

  ewsb->def_property("mr1", &EwsbParameters::get_mr1, &EwsbParameters::set_mr1,
                     "Mass of the 1st gen RH slepton : ℒ ⊃ -Ml₁ * |lR₁|²");
  ewsb->def_property("mr2", &EwsbParameters::get_mr2, &EwsbParameters::set_mr2,
                     "Mass of the 2nd gen RH slepton : ℒ ⊃ -Ml₂ * |lR₂|²");
  ewsb->def_property("mr3", &EwsbParameters::get_mr3, &EwsbParameters::set_mr3,
                     "Mass of the 3rd gen RH slepton : ℒ ⊃ -Ml₃ * |lR₃|²");

  // ===========================
  // ---- LH Squark Doublet ----
  // ===========================

  ewsb->def_property("mq1", &EwsbParameters::get_mq1, &EwsbParameters::set_mq1,
                     "Mass of the 1st gen squark doublet: ℒ ⊃ -MQ₁ * ⟨Q₁, Q₁⟩");
  ewsb->def_property("mq2", &EwsbParameters::get_mq2, &EwsbParameters::set_mq2,
                     "Mass of the 2nd gen squark doublet: ℒ ⊃ -MQ₂ * ⟨Q₂, Q₂⟩");
  ewsb->def_property("mq3", &EwsbParameters::get_mq3, &EwsbParameters::set_mq3,
                     "Mass of the 3rd gen squark doublet: ℒ ⊃ -MQ₃ * ⟨Q₃, Q₃⟩");

  // ======================
  // ---- RH up-squark ----
  // ======================

  ewsb->def_property("mu1", &EwsbParameters::get_mu1, &EwsbParameters::set_mu1,
                     "Mass of the 1st gen RH up-squark : ℒ ⊃ -Mu₁ * |uR₁|²");
  ewsb->def_property("mu2", &EwsbParameters::get_mu2, &EwsbParameters::set_mu2,
                     "Mass of the 2nd gen RH up-squark : ℒ ⊃ -Mu₂ * |uR₂|²");
  ewsb->def_property("mu3", &EwsbParameters::get_mu3, &EwsbParameters::set_mu3,
                     "Mass of the 3rd gen RH up-squark : ℒ ⊃ -Mu₃ * |uR₃|²");

  // ========================
  // ---- RH down-squark ----
  // ========================

  ewsb->def_property("md1", &EwsbParameters::get_md1, &EwsbParameters::set_md1,
                     "Mass of the 1st gen RH down-squark : ℒ ⊃ -Md₁ * |dR₁|²");
  ewsb->def_property("md2", &EwsbParameters::get_md2, &EwsbParameters::set_md2,
                     "Mass of the 2nd gen RH down-squark : ℒ ⊃ -Md₂ * |dR₂|²");
  ewsb->def_property("md3", &EwsbParameters::get_md3, &EwsbParameters::set_md3,
                     "Mass of the 3rd gen RH down-squark : ℒ ⊃ -Md₃ * |dR₃|²");

  // ==========================
  // ---- Higgs Parameters ----
  // ==========================

  ewsb->def_property(
      "mu", &EwsbParameters::get_mu, &EwsbParameters::set_mu,
      "The coefficient of the Higgs superpotential term: ℒ ⊃ μ Hu Hd.");

  ewsb->def_property("mh3", &EwsbParameters::get_mh3, &EwsbParameters::set_mh3,
                     "Mass of the psuedo-scalar Higgs boson MA");
  ewsb->def_property("tb", &EwsbParameters::get_tb, &EwsbParameters::set_tb,
                     "Tangent beta: tan(β) = v₂ / v₁");

  // =============================
  // ---- Trilinear Couplings ----
  // =============================

  ewsb->def_property("at", &EwsbParameters::get_at, &EwsbParameters::set_at,
                     "Trilinear top coupling: ℒ ⊃ At Yt tR ⟨Hu, Q₃⟩");
  ewsb->def_property("ab", &EwsbParameters::get_ab, &EwsbParameters::set_ab,
                     "Trilinear bottom coupling: ℒ ⊃ Ab Yb bR ⟨Hd, Q₃⟩");
  ewsb->def_property("al", &EwsbParameters::get_al, &EwsbParameters::set_al,
                     "Trilinear tau coupling: ℒ ⊃ Al Yl lR ⟨Hd, Q₃⟩");
  ewsb->def_property("au", &EwsbParameters::get_au, &EwsbParameters::set_au,
                     "Trilinear 1st/2nd gen up-quark couplings: ℒ ⊃ Au [ Yu "
                     "uR ⟨Hu, Q₁⟩ + Yc cR ⟨Hu, Q₂⟩ ]");
  ewsb->def_property("ad", &EwsbParameters::get_ad, &EwsbParameters::set_ad,
                     "Trilinear 1st/2nd gen down-quark couplings: ℒ ⊃ Ad [ "
                     "Yd dR ⟨Hd, Q₁⟩ + Ys sR ⟨Hd, Q₂⟩ ]");

  ewsb->def_property("am", &EwsbParameters::get_am, &EwsbParameters::set_am);

  // ======================
  // ---- SM Couplings ----
  // ======================

  ewsb->def_property("alpha_s_mz", &EwsbParameters::get_alfsmz,
                     &EwsbParameters::set_alfsmz,
                     "Value of αQCD at Z-Boson mass.");

  ewsb->def_property("top_quark_pole_mass", &EwsbParameters::get_mtp,
                     &EwsbParameters::set_mtp, "Top-Quark pole mass.");

  ewsb->def_property("bottom_mass", &EwsbParameters::get_mbmb,
                     &EwsbParameters::set_mbmb,
                     "Bottom-quark mass evaluated at bottom-quark pole mass.");

  ewsb->def_property("charm_mass", &EwsbParameters::get_mcmc,
                     &EwsbParameters::set_mcmc,
                     "Charm-quark mass evaluated at charm-quark pole mass.");

  ewsb->def_property("electric_charge", &EwsbParameters::get_ee,
                     &EwsbParameters::set_ee,
                     "Value of the electromagnetic coupling 'e'.");

  ewsb->def_property("w_boson_mass", &EwsbParameters::get_mw,
                     &EwsbParameters::set_mw,
                     "Value of the W-Boson pole mass.");

  ewsb->def_property("z_boson_mass", &EwsbParameters::get_mz,
                     &EwsbParameters::set_mz,
                     "Value of the Z-Boson pole mass.");

  ewsb->def_property("qcd_scale", &EwsbParameters::get_q,
                     &EwsbParameters::set_q,
                     "QCD scale for running quark masses.");

  ewsb->def_property("tau_mass", &EwsbParameters::get_ml,
                     &EwsbParameters::set_ml, "Mass of the tau-lepton.");

  ewsb->def_property("light_quark_mass", &EwsbParameters::get_mq,
                     &EwsbParameters::set_mq, "Mass of the light-quarks.");
}
