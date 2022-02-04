#include "execute.hpp"
#include "ewsb.hpp"
#include "results.hpp"
#include "settings.hpp"
#include "sugra.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

void execute(MicromegasResults *results, const MicromegasSettings &settings) {
  results->compute_relic_density(settings);
  results->compute_masses(settings);
  results->compute_gmuon(settings);
  results->compute_bsg(settings);
  results->compute_bsmumu(settings);
  results->compute_btaunu(settings);
  results->compute_deltarho(settings);
  results->compute_rl23(settings);
  results->compute_d_taunu_and_munu(settings);
  results->compute_masslimits(settings);
  results->compute_nucleon_amplitudes(settings);
  results->compute_direct_detection_pvalues(settings);
  results->compute_z_invisible(settings);
  results->compute_lsp_nlsp_lep(settings);
  results->compute_z_prime_limits(settings);
  results->compute_monojet(settings);
}

void execute(MicromegasResults *results, const MicromegasSettings &settings,
             const SugraParameters &sugra) {
  using micromegas::Micromegas;

  if (PyErr_CheckSignals() != 0) {
    throw py::error_already_set();
  }

  try {
    const double mhf = sugra.get_mhf();
    const double m0 = sugra.get_m0();
    const double a0 = sugra.get_a0();
    Micromegas::mssm_sugra(sugra.get_tb(), mhf, mhf, mhf, a0, a0, a0,
                           sugra.get_sgn(), m0, m0, m0, m0, m0, m0, m0, m0, m0,
                           m0, m0, m0);
    Micromegas::sort_odd_particles();
  } catch (const std::exception &e) {
    results->set_nans();
    if (settings.get_debug()) {
      py::print(e.what());
    }
  }

  execute(results, settings);
}

void execute(MicromegasResults *results, const MicromegasSettings &settings,
             const EwsbParameters &ewsb) {
  using micromegas::Micromegas;

  if (PyErr_CheckSignals() != 0) {
    throw py::error_already_set();
  }

  try {
    ewsb.assign_all();
    Micromegas::mssm_ewsb();
    micromegas::Micromegas::sort_odd_particles();
  } catch (const std::exception &e) {
    results->set_nans();
    if (settings.get_debug()) {
      py::print(e.what());
    }
  }

  execute(results, settings);
}

MicromegasResults execute(const MicromegasSettings &settings,
                          const SugraParameters &sugra) {
  MicromegasResults results(1);
  execute(&results, settings, sugra);
  return results;
}

MicromegasResults execute(const MicromegasSettings &settings,
                          const EwsbParameters &ewsb) {
  MicromegasResults results(1);
  execute(&results, settings, ewsb);
  return results;
}

MicromegasResults execute(const MicromegasSettings &settings,
                          const std::vector<SugraParameters> &sugras) {
  MicromegasResults results(sugras.size());
  for (const auto &sugra : sugras) { // NOLINT
    execute(&results, settings, sugra);
  }
  return results;
}

MicromegasResults execute(const MicromegasSettings &settings,
                          const std::vector<EwsbParameters> &ewsbs) {
  MicromegasResults results(ewsbs.size());
  for (const auto &ewsb : ewsbs) { // NOLINT
    execute(&results, settings, ewsb);
  }
  return results;
}
