#include "pymicromegas.hpp"

namespace py = pybind11;

void define_results(py::class_<MicromegasResults> *results) {

  results->def_property(
      "omega", [](const MicromegasResults &r) { return r.omega(); }, nullptr,
      "The dark matter relic density Ωh².");

  results->def_property(
      "xf", [](const MicromegasResults &r) { return r.xf(); }, nullptr,
      "The scaled dark matter freeze-out temperature: xᶠ = mᵡ/Tᶠ.");

  results->def_property(
      "bsgsm", [](const MicromegasResults &r) { return r.b_sg_sm(); }, nullptr,
      "The Standard Model value of the branching ratio for b → sγ");

  results->def_property(
      "bsgnlo", [](const MicromegasResults &r) { return r.b_sg_nlo(); },
      nullptr,
      "The value of the branching ratio for b → sγ, including contributions "
      "beyond leading-order.");

  results->def_property(
      "deltarho", [](const MicromegasResults &r) { return r.delta_rho(); },
      nullptr,
      "The ∆ρ parameter in the MSSM. It contains for example the stop/sbottom "
      "contributions, as well as the two-loop QCD corrections due to gluon "
      "exchange and the correction due to gluino exchange in the heavy gluino "
      "limit."

  );

  results->def_property(
      "bsmumu", [](const MicromegasResults &r) { return r.b_smumu(); }, nullptr,
      "The value of the branching ratio Bs → μ+μ− in the MSSM. It includes the "
      "loop contributions due to chargino, sneutrino, stop and Higgs exchange. "
      "The ∆mb effect relevant for high tan β is taken into account."

  );

  results->def_property(
      "bstaunu", [](const MicromegasResults &r) { return r.b_taunu(); },
      nullptr,
      "The ratio between the MSSM and SM branching fractions for  ̄B+ → τ +ντ "
      ".");

  results->def_property(
      "gmuon", [](const MicromegasResults &r) { return r.g_muon(); }, nullptr,
      "The value of the supersymmetric contribution to the anomalous magnetic "
      "moment of the muon.");

  results->def_property(
      "rl23", [](const MicromegasResults &r) { return r.rl23(); }, nullptr,
      "The ratio of the MSSM to SM value for Rl23 in K⁺ → μ⁺ + νμ due to a "
      "charged "
      "higgs contribution.");

  results->def_property(
      "dtaunu", [](const MicromegasResults &r) { return r.dtaunu(); }, nullptr,
      "The branching ratio for Ds⁺ → τ⁺ + ντ");

  results->def_property(
      "dmunu", [](const MicromegasResults &r) { return r.dmunu(); }, nullptr,
      "The branching ratio for Ds⁺ → μ⁺ + νμ");

  results->def_property(
      "masslimits", [](const MicromegasResults &r) { return r.masslimits(); },
      nullptr);

  results->def_property(
      "proton_si_amp",
      [](const MicromegasResults &r) { return r.proton_si_amp(); }, nullptr,
      "DM-proton spin-independent amplitude.");

  results->def_property(
      "proton_sd_amp",
      [](const MicromegasResults &r) { return r.proton_sd_amp(); }, nullptr,
      "DM-proton spin-dependent amplitude.");

  results->def_property(
      "neutron_si_amp",
      [](const MicromegasResults &r) { return r.neutron_si_amp(); }, nullptr,
      "DM-neutron spin-independent amplitude.");

  results->def_property(
      "neutron_sd_amp",
      [](const MicromegasResults &r) { return r.neutron_sd_amp(); }, nullptr,
      "DM-neutron spin-dependent amplitude.");

  results->def_property(
      "pval_xenon1T",
      [](const MicromegasResults &r) { return r.pval_xenon1T(); }, nullptr,
      "p-value from the XENON1T direct-detection experiment.");

  results->def_property(
      "pval_cresst", [](const MicromegasResults &r) { return r.pval_cresst(); },
      nullptr, "p-value from the CRESST-III direct-detection experiment.");

  results->def_property(
      "pval_darkside",
      [](const MicromegasResults &r) { return r.pval_darkside(); }, nullptr,
      "p-value from the DarkSide direct-detection experiment.");

  results->def_property(
      "pval_pico", [](const MicromegasResults &r) { return r.pval_pico(); },
      nullptr, "p-value from PICO direct-detection experiment.");

  // ==========================
  // ---- selectron Masses ----
  // ==========================

  results->def_property(
      "msel", [](const MicromegasResults &r) { return r.msel(); }, nullptr,
      "Left-handed selectron mass.");
  results->def_property(
      "mser", [](const MicromegasResults &r) { return r.mser(); }, nullptr,
      "Right-handed selectron mass.");
  results->def_property(
      "msne", [](const MicromegasResults &r) { return r.msne(); }, nullptr,
      "Electron-snuetrino mass.");

  // ======================
  // ---- smuon masses ----
  // ======================

  results->def_property(
      "msml", [](const MicromegasResults &r) { return r.msml(); }, nullptr,
      "Left-handed smuon mass.");
  results->def_property(
      "msmr", [](const MicromegasResults &r) { return r.msmr(); }, nullptr,
      "Right-handed smuon mass.");
  results->def_property(
      "msnm", [](const MicromegasResults &r) { return r.msnm(); }, nullptr,
      "Muon-snuetrino mass.");

  // =====================
  // ---- stau masses ----
  // =====================

  results->def_property(
      "msl1", [](const MicromegasResults &r) { return r.msl1(); }, nullptr,
      "Light stau mass.");
  results->def_property(
      "msl2", [](const MicromegasResults &r) { return r.msl2(); }, nullptr,
      "Heavy stau mass.");
  results->def_property(
      "msnl", [](const MicromegasResults &r) { return r.msnl(); }, nullptr,
      "Left-handed tau snuetrino mass.");

  // ==========================
  // ---- up squark masses ----
  // ==========================

  results->def_property(
      "msul", [](const MicromegasResults &r) { return r.msul(); }, nullptr,
      "Left-handed up-squark mass.");
  results->def_property(
      "msur", [](const MicromegasResults &r) { return r.msur(); }, nullptr,
      "Right-handed up-squark mass.");

  // ============================
  // ---- down squark masses ----
  // ============================

  results->def_property(
      "msdl", [](const MicromegasResults &r) { return r.msdl(); }, nullptr,
      "Left-handed down squark mass.");
  results->def_property(
      "msdr", [](const MicromegasResults &r) { return r.msdr(); }, nullptr,
      "Right-handed down squark mass.");

  // =============================
  // ---- charm squark masses ----
  // =============================

  results->def_property(
      "mscl", [](const MicromegasResults &r) { return r.mscl(); }, nullptr,
      "Left-handed charm squark mass.");
  results->def_property(
      "mscr", [](const MicromegasResults &r) { return r.mscr(); }, nullptr,
      "Right-handed charm squark mass.");

  // ===============================
  // ---- strange squark masses ----
  // ===============================

  results->def_property(
      "mssl", [](const MicromegasResults &r) { return r.mssl(); }, nullptr,
      "Left-handed strange squark mass.");
  results->def_property(
      "mssr", [](const MicromegasResults &r) { return r.mssr(); }, nullptr,
      "Right-handed strange squark mass.");

  // ===========================
  // ---- top squark masses ----
  // ===========================

  results->def_property(
      "mst1", [](const MicromegasResults &r) { return r.mst1(); }, nullptr,
      "Light stop mass.");
  results->def_property(
      "mst2", [](const MicromegasResults &r) { return r.mst2(); }, nullptr,
      "Heavy stop mass.");

  // ==============================
  // ---- bottom squark masses ----
  // ==============================

  results->def_property(
      "msb1", [](const MicromegasResults &r) { return r.msb1(); }, nullptr,
      "Light sbottom mass.");
  results->def_property(
      "msb2", [](const MicromegasResults &r) { return r.msb2(); }, nullptr,
      "Heavy sbottom mass.");

  // ===========================
  // ---- neutralino masses ----
  // ===========================

  results->def_property(
      "mneut1", [](const MicromegasResults &r) { return r.mneut1(); }, nullptr,
      "Lightest neutralino mass.");
  results->def_property(
      "mneut2", [](const MicromegasResults &r) { return r.mneut2(); }, nullptr,
      "Second neutralino mass.");
  results->def_property(
      "mneut3", [](const MicromegasResults &r) { return r.mneut3(); }, nullptr,
      "Third neutralino mass.");
  results->def_property(
      "mneut4", [](const MicromegasResults &r) { return r.mneut4(); }, nullptr,
      "Heaviest neutralino mass.");

  // =========================
  // ---- chargino masses ----
  // =========================

  results->def_property(
      "mchg1", [](const MicromegasResults &r) { return r.mchg1(); }, nullptr,
      "Light chargino mass.");
  results->def_property(
      "mchg2", [](const MicromegasResults &r) { return r.mchg2(); }, nullptr,
      "Heavy chargino mass.");

  // ======================
  // ---- Higgs masses ----
  // ======================

  results->def_property(
      "mhsm", [](const MicromegasResults &r) { return r.mhsm(); }, nullptr,
      "Higgs mass.");
  results->def_property(
      "mh", [](const MicromegasResults &r) { return r.mh(); }, nullptr,
      "Other higgs mass");
  results->def_property(
      "mhc", [](const MicromegasResults &r) { return r.mhc(); }, nullptr,
      "Charged Higgs mass");

  // ========================
  // ---- Gaugino masses ----
  // ========================

  results->def_property(
      "mg", [](const MicromegasResults &r) { return r.mg(); }, nullptr,
      "Gaugino mass");
}
