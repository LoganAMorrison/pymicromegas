#include "results.hpp"
#include "micromegas.hpp"
#include <iostream>

void MicromegasResults::set_nans() {
  p_omega.push_back(std::nan("0"));
  p_xf.push_back(std::nan("0"));
  msne().push_back(std::nan("0"));
  msnm().push_back(std::nan("0"));
  msel().push_back(std::nan("0"));
  mser().push_back(std::nan("0"));
  msml().push_back(std::nan("0"));
  msmr().push_back(std::nan("0"));
  msdl().push_back(std::nan("0"));
  msdr().push_back(std::nan("0"));
  msul().push_back(std::nan("0"));
  msur().push_back(std::nan("0"));
  mssl().push_back(std::nan("0"));
  mssr().push_back(std::nan("0"));
  mscl().push_back(std::nan("0"));
  mscr().push_back(std::nan("0"));
  msnl().push_back(std::nan("0"));
  msl1().push_back(std::nan("0"));
  msl2().push_back(std::nan("0"));
  msb1().push_back(std::nan("0"));
  msb2().push_back(std::nan("0"));
  mst1().push_back(std::nan("0"));
  mst2().push_back(std::nan("0"));
  mg().push_back(std::nan("0"));
  mneut1().push_back(std::nan("0"));
  mneut2().push_back(std::nan("0"));
  mneut3().push_back(std::nan("0"));
  mneut4().push_back(std::nan("0"));
  mchg1().push_back(std::nan("0"));
  mchg2().push_back(std::nan("0"));
  mhsm().push_back(std::nan("0"));
  mh().push_back(std::nan("0"));
  mhc().push_back(std::nan("0"));

  p_monojet.push_back(std::nan("0"));
  p_z_prime_limits.push_back(std::nan("0"));
  p_lsp_nlsp_lep.push_back(std::nan("0"));
  p_z_invisible.push_back(std::nan("0"));
  p_pval_xenon1T.push_back(std::nan("0"));
  p_pval_cresst.push_back(std::nan("0"));
  p_pval_darkside.push_back(std::nan("0"));
  p_pval_pico.push_back(std::nan("0"));
  p_proton_si_amp.push_back(std::nan("0"));
  p_proton_sd_amp.push_back(std::nan("0"));
  p_neutron_si_amp.push_back(std::nan("0"));
  p_neutron_sd_amp.push_back(std::nan("0"));
  p_masslimits.push_back(std::nan("0"));
  p_dtaunu.push_back(std::nan("0"));
  p_dmunu.push_back(std::nan("0"));
  p_rl23.push_back(std::nan("0"));
  p_deltarho.push_back(std::nan("0"));
  p_btaunu.push_back(std::nan("0"));
  p_bsmumu.push_back(std::nan("0"));
  p_bsgnlo.push_back(std::nan("0"));
  p_bsgsm.push_back(std::nan("0"));
  p_gmuon.push_back(std::nan("0"));
}

void MicromegasResults::compute_relic_density(
    const MicromegasSettings &settings) {
  if (settings.get_relic_density()) {
    auto rd = micromegas::Micromegas::relic_density(settings.get_fast(),
                                                    settings.get_beps());
    p_omega.push_back(rd.first);
    p_xf.push_back(rd.second);
  } else {
    p_omega.push_back(std::nan("0"));
    p_xf.push_back(std::nan("0"));
  }
}

void MicromegasResults::compute_masses(const MicromegasSettings &settings) {
  using micromegas::Micromegas;
  if (settings.get_masses()) {
    msne().push_back(Micromegas::find_val("MSne"));
    msnm().push_back(Micromegas::find_val("MSnm"));
    msel().push_back(Micromegas::find_val("MSeL"));
    mser().push_back(Micromegas::find_val("MSeR"));
    msml().push_back(Micromegas::find_val("MSmL"));
    msmr().push_back(Micromegas::find_val("MSmR"));
    msdl().push_back(Micromegas::find_val("MSdL"));
    msdr().push_back(Micromegas::find_val("MSdR"));
    msul().push_back(Micromegas::find_val("MSuL"));
    msur().push_back(Micromegas::find_val("MSuR"));
    mssl().push_back(Micromegas::find_val("MSsL"));
    mssr().push_back(Micromegas::find_val("MSsR"));
    mscl().push_back(Micromegas::find_val("MScL"));
    mscr().push_back(Micromegas::find_val("MScR"));
    msnl().push_back(Micromegas::find_val("MSnl"));
    msl1().push_back(Micromegas::find_val("MSl1"));
    msl2().push_back(Micromegas::find_val("MSl2"));
    msb1().push_back(Micromegas::find_val("MSb1"));
    msb2().push_back(Micromegas::find_val("MSb2"));
    mst1().push_back(Micromegas::find_val("MSt1"));
    mst2().push_back(Micromegas::find_val("MSt2"));
    mg().push_back(Micromegas::find_val("MSG"));
    mneut1().push_back(Micromegas::find_val("MNE1"));
    mneut2().push_back(Micromegas::find_val("MNE2"));
    mneut3().push_back(Micromegas::find_val("MNE3"));
    mneut4().push_back(Micromegas::find_val("MNE4"));
    mchg1().push_back(Micromegas::find_val("MC1"));
    mchg2().push_back(Micromegas::find_val("MC2"));
    mhsm().push_back(Micromegas::find_val("Mh"));
    mh().push_back(Micromegas::find_val("MH"));
    mhc().push_back(Micromegas::find_val("MHc"));
  } else {
    msne().push_back(std::nan("0"));
    msnm().push_back(std::nan("0"));
    msel().push_back(std::nan("0"));
    mser().push_back(std::nan("0"));
    msml().push_back(std::nan("0"));
    msmr().push_back(std::nan("0"));
    msdl().push_back(std::nan("0"));
    msdr().push_back(std::nan("0"));
    msul().push_back(std::nan("0"));
    msur().push_back(std::nan("0"));
    mssl().push_back(std::nan("0"));
    mssr().push_back(std::nan("0"));
    mscl().push_back(std::nan("0"));
    mscr().push_back(std::nan("0"));
    msnl().push_back(std::nan("0"));
    msl1().push_back(std::nan("0"));
    msl2().push_back(std::nan("0"));
    msb1().push_back(std::nan("0"));
    msb2().push_back(std::nan("0"));
    mst1().push_back(std::nan("0"));
    mst2().push_back(std::nan("0"));
    mg().push_back(std::nan("0"));
    mneut1().push_back(std::nan("0"));
    mneut2().push_back(std::nan("0"));
    mneut3().push_back(std::nan("0"));
    mneut4().push_back(std::nan("0"));
    mchg1().push_back(std::nan("0"));
    mchg2().push_back(std::nan("0"));
    mhsm().push_back(std::nan("0"));
    mh().push_back(std::nan("0"));
    mhc().push_back(std::nan("0"));
  }
}

void MicromegasResults::compute_gmuon(const MicromegasSettings &settings) {
  if (settings.get_gmuon()) {
    auto res = micromegas::Micromegas::mssm_gmuon();
    p_gmuon.push_back(res);
  } else {
    p_gmuon.push_back(std::nan("0"));
  }
}

void MicromegasResults::compute_bsg(const MicromegasSettings &settings) {
  if (settings.get_bsg()) {
    auto res = micromegas::Micromegas::mssm_bsgnlo();
    p_bsgnlo.push_back(res.first);
    p_bsgsm.push_back(res.second);
  } else {
    p_bsgnlo.push_back(std::nan("0"));
    p_bsgsm.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_bsmumu(const MicromegasSettings &settings) {
  if (settings.get_bsmumu()) {
    auto res = micromegas::Micromegas::mssm_bsmumu();
    p_bsmumu.push_back(res);
  } else {
    p_bsmumu.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_btaunu(const MicromegasSettings &settings) {

  if (settings.get_btaunu()) {
    auto res = micromegas::Micromegas::mssm_btaunu();
    p_btaunu.push_back(res);
  } else {
    p_btaunu.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_deltarho(const MicromegasSettings &settings) {

  if (settings.get_deltarho()) {
    auto res = micromegas::Micromegas::mssm_deltarho();
    p_deltarho.push_back(res);
  } else {
    p_deltarho.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_rl23(const MicromegasSettings &settings) {

  if (settings.get_rl23()) {
    auto res = micromegas::Micromegas::mssm_rl23();
    p_rl23.push_back(res);
  } else {
    p_rl23.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_d_taunu_and_munu(
    const MicromegasSettings &settings) {

  if (settings.get_d_taunu_and_munu()) {
    auto res = micromegas::Micromegas::mssm_d_taunu_and_munu();
    p_dtaunu.push_back(res.first);
    p_dmunu.push_back(res.second);
  } else {
    p_dtaunu.push_back(std::nan("0"));
    p_dmunu.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_masslimits(const MicromegasSettings &settings) {

  if (settings.get_masslimits()) {
    auto res = micromegas::Micromegas::mssm_masslimits();
    p_masslimits.push_back(res);
  } else {
    p_masslimits.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_nucleon_amplitudes(
    const MicromegasSettings &settings) {

  if (settings.get_nucleon_amplitudes()) {
    auto res = micromegas::Micromegas::nucleon_amplitudes_cdm1();
    p_proton_si_amp.push_back(res.proton_si[0]);
    p_proton_sd_amp.push_back(res.proton_sd[0]);
    p_neutron_si_amp.push_back(res.neutron_si[0]);
    p_neutron_sd_amp.push_back(res.neutron_sd[0]);
  } else {
    p_proton_si_amp.push_back(std::nan("0"));
    p_proton_sd_amp.push_back(std::nan("0"));
    p_neutron_si_amp.push_back(std::nan("0"));
    p_neutron_sd_amp.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_direct_detection_pvalues(
    const MicromegasSettings &settings) {

  if (settings.get_direct_detection_pvalues()) {
    auto res = micromegas::Micromegas::direct_detection_pval_maxwell();
    p_pval_xenon1T.push_back(res.xenon1T);
    p_pval_cresst.push_back(res.cresst);
    p_pval_darkside.push_back(res.darkside);
    p_pval_pico.push_back(res.pico);
  } else {
    p_pval_xenon1T.push_back(std::nan("0"));
    p_pval_cresst.push_back(std::nan("0"));
    p_pval_darkside.push_back(std::nan("0"));
    p_pval_pico.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_z_invisible(
    const MicromegasSettings &settings) {

  if (settings.get_z_invisible()) {
    // auto res = micromegas::Micromegas::z_invisible();
    std::cout << "z_invisible not yet implemented.";
    p_z_invisible.push_back(std::nan("0"));
  } else {
    p_z_invisible.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_lsp_nlsp_lep(
    const MicromegasSettings &settings) {

  if (settings.get_lsp_nlsp_lep()) {
    // auto res = micromegas::Micromegas::lsp_nlsp_lep();
    std::cout << "lsp_nlsp_lep not yet implemented.";
    p_lsp_nlsp_lep.push_back(std::nan("0"));
  } else {
    p_lsp_nlsp_lep.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_z_prime_limits(
    const MicromegasSettings &settings) {
  if (settings.get_z_prime_limits()) {
    // auto res = micromegas::Micromegas::z_prime_limits();
    std::cout << "z_prime_limits not yet implemented.";
    p_z_prime_limits.push_back(std::nan("0"));
  } else {
    p_z_prime_limits.push_back(std::nan("0"));
  }
}
void MicromegasResults::compute_monojet(const MicromegasSettings &settings) {
  if (settings.get_monojet()) {
    auto res = micromegas::Micromegas::monojet();
    p_monojet.push_back(res);
  } else {
    p_monojet.push_back(std::nan("0"));
  }
}
