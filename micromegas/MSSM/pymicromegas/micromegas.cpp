#include "micromegas.hpp"
#include "../../include/micromegas.h"
#include "../../include/micromegas_aux.h"
#include "../lib/pmodel.h"
#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

namespace micromegas {

static std::vector<char> string_to_char_vec(const std::string &str) {
  return {str.begin(), str.end()};
}

/**
 * Set the value of VZdecay and VWdecay.
 *
 * @param flag If 0, off-shell Z/W decays are turned off. If 1, 3-body
 * final-states are included. If 2, Z's and W's are included in
 * coannihilations as well.
 */
void Micromegas::set_v_decay(int flag) {
  if (flag == 0 || flag == 1 || flag == 2) {
    VWdecay = flag;
    VZdecay = flag;
  } else {
    throw std::runtime_error("Invalid flag for VVdecay.");
  }
}

/**
 * Return the value of micromegas CDM1.
 */
std::string Micromegas::get_cdm1() { return CDM1 == nullptr ? "null" : CDM1; }

/**
 * Return the value of micromegas CDM2.
 */
std::string Micromegas::get_cdm2() { return CDM2 == nullptr ? "null" : CDM2; }

/**
 * Return the mass of micromegas first CDM.
 */
double Micromegas::get_mcdm1() {
  return CDM1 == nullptr ? std::nan("0") : Mcdm1;
}

/**
 * Return the mass of micromegas second CDM.
 */
double Micromegas::get_mcdm2() {
  return CDM2 == nullptr ? std::nan("0") : Mcdm2;
}

/**
 * Assign value val to parameter name. Returns a non-zero
 * value if it cannot recognize a parameter name.
 */
int Micromegas::assign_val(const std::string &name, double val) {
  return assignVal(name.c_str(), val);
}

/**
 * Assign value val to parameter name. Writes an error message if it cannot
 * recognize a parameter name.
 */
int Micromegas::assign_val_w(const std::string &name, double val) {
  return assignValW(name.c_str(), val);
}

/**
 * Assign value val to parameter name. Returns a non-zero
 * value if it cannot recognize a parameter name.
 */
int Micromegas::mssm_assign_val(MssmParameter param, double val) {
  switch (param) {
  case MssmParameter::Mu:
    return assignVal("mu", val);
  case MssmParameter::Mg1:
    return assignVal("MG1", val);
  case MssmParameter::Mg2:
    return assignVal("MG2", val);
  case MssmParameter::Mg3:
    return assignVal("MG3", val);
  case MssmParameter::Ml1:
    return assignVal("Ml1", val);
  case MssmParameter::Ml2:
    return assignVal("Ml2", val);
  case MssmParameter::Ml3:
    return assignVal("Ml3", val);
  case MssmParameter::Mr1:
    return assignVal("Mr1", val);
  case MssmParameter::Mr2:
    return assignVal("Mr2", val);
  case MssmParameter::Mr3:
    return assignVal("Mr3", val);
  case MssmParameter::Mq1:
    return assignVal("Mq1", val);
  case MssmParameter::Mq2:
    return assignVal("Mq2", val);
  case MssmParameter::Mq3:
    return assignVal("Mq3", val);
  case MssmParameter::Mu1:
    return assignVal("Mu1", val);
  case MssmParameter::Mu2:
    return assignVal("Mu2", val);
  case MssmParameter::Mu3:
    return assignVal("Mu3", val);
  case MssmParameter::Md1:
    return assignVal("Md1", val);
  case MssmParameter::Md2:
    return assignVal("Md2", val);
  case MssmParameter::Md3:
    return assignVal("Md3", val);
  case MssmParameter::Mh3:
    return assignVal("MH3", val);
  case MssmParameter::Tb:
    return assignVal("tb", val);
  case MssmParameter::At:
    return assignVal("At", val);
  case MssmParameter::Ab:
    return assignVal("Ab", val);
  case MssmParameter::Al:
    return assignVal("Al", val);
  case MssmParameter::Au:
    return assignVal("Au", val);
  case MssmParameter::Ad:
    return assignVal("Ad", val);
  case MssmParameter::Alfsmz:
    return assignVal("alfSMZ", val);
  case MssmParameter::Mz:
    return assignVal("MZ", val);
  case MssmParameter::Mw:
    return assignVal("MW", val);
  case MssmParameter::Mtp:
    return assignVal("Mtp", val);
  case MssmParameter::Mbmb:
    return assignVal("MbMb", val);
  case MssmParameter::Mcmc:
    return assignVal("McMc", val);
  case MssmParameter::Gg:
    return assignVal("GG", val);
  case MssmParameter::Q:
    return assignVal("Q", val);
  case MssmParameter::Ee:
    return assignVal("EE", val);
  case MssmParameter::Ml:
    return assignVal("Ml", val);
  case MssmParameter::Mq:
    return assignVal("Mq", val);
  case MssmParameter::Am:
    return assignVal("Am", val);
  default:
    throw std::runtime_error("Can't get here.");
  }
}

/**
 * Assign value val to parameter name. Writes an error message if it cannot
 * recognize a parameter name.
 */
int Micromegas::mssm_assign_val_w(MssmParameter param, double val) {
  switch (param) {
  case MssmParameter::Mu:
    return assignValW("mu", val);
  case MssmParameter::Mg1:
    return assignValW("MG1", val);
  case MssmParameter::Mg2:
    return assignValW("MG2", val);
  case MssmParameter::Mg3:
    return assignValW("MG3", val);
  case MssmParameter::Ml1:
    return assignValW("ML1", val);
  case MssmParameter::Ml2:
    return assignValW("ML2", val);
  case MssmParameter::Ml3:
    return assignValW("ML3", val);
  case MssmParameter::Mr1:
    return assignValW("MR1", val);
  case MssmParameter::Mr2:
    return assignValW("MR2", val);
  case MssmParameter::Mr3:
    return assignValW("MR3", val);
  case MssmParameter::Mq1:
    return assignValW("Mq1", val);
  case MssmParameter::Mq2:
    return assignValW("Mq2", val);
  case MssmParameter::Mq3:
    return assignValW("Mq3", val);
  case MssmParameter::Mu1:
    return assignValW("Mu1", val);
  case MssmParameter::Mu2:
    return assignValW("Mu2", val);
  case MssmParameter::Mu3:
    return assignValW("Mu3", val);
  case MssmParameter::Md1:
    return assignValW("Md1", val);
  case MssmParameter::Md2:
    return assignValW("Md2", val);
  case MssmParameter::Md3:
    return assignValW("Md3", val);
  case MssmParameter::Mh3:
    return assignValW("MH3", val);
  case MssmParameter::Tb:
    return assignValW("tb", val);
  case MssmParameter::At:
    return assignValW("At", val);
  case MssmParameter::Ab:
    return assignValW("Ab", val);
  case MssmParameter::Al:
    return assignValW("Al", val);
  case MssmParameter::Au:
    return assignValW("Au", val);
  case MssmParameter::Ad:
    return assignValW("Ad", val);
  case MssmParameter::Alfsmz:
    return assignValW("alfSMZ", val);
  case MssmParameter::Mz:
    return assignValW("MZ", val);
  case MssmParameter::Mw:
    return assignValW("MW", val);
  case MssmParameter::Mtp:
    return assignValW("Mtp", val);
  case MssmParameter::Mbmb:
    return assignValW("MbMb", val);
  case MssmParameter::Mcmc:
    return assignValW("McMc", val);
  case MssmParameter::Gg:
    return assignValW("GG", val);
  case MssmParameter::Q:
    return assignValW("Q", val);
  case MssmParameter::Ee:
    return assignValW("EE", val);
  case MssmParameter::Ml:
    return assignValW("Ml", val);
  case MssmParameter::Mq:
    return assignValW("Mq", val);
  case MssmParameter::Am:
    return assignValW("Am", val);
  default:
    throw std::runtime_error("Can't get here.");
  }
}

/**
 * Finds the value of variable name and assigns it to parameter val. It
 * returns a non-zero value if it cannot recognize a parameter name.
 */
double Micromegas::find_val(const std::string &name) {
  double val = 0.0;
  char *p = new char[name.size() + 1]; // NOLINT
  name.copy(p, name.size());
  p[name.size()] = '\0'; // NOLINT
  const int err = findVal(p, &val);
  delete[] p; // NOLINT
  if (err != 0) {
    return std::nan("0");
  }

  return val;
}

/**
 * Returns the value of variable name and writes an error message if it cannot
 * recognize a parameter name.
 */
double Micromegas::find_val_w(const std::string &fname) {
  char *p = new char[fname.size() + 1]; // NOLINT
  fname.copy(p, fname.size());
  p[fname.size()] = '\0'; // NOLINT
  const double result = findValW(p);
  delete[] p; // NOLINT
  return result;
}

/**
 * Finds the value of variable name and assigns it to parameter val. It
 * returns a non-zero value if it cannot recognize a parameter name.
 */
double Micromegas::mssm_find_val(MssmParameter param) {
  double val = 0.0;

  switch (param) {
  case MssmParameter::Mu: {
    char p[] = "mu";
    findVal(p, &val);
  } break;
  case MssmParameter::Mg1: {
    char p[] = "MG1"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mg2: {
    char p[] = "MG2"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mg3: {
    char p[] = "MG3"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Ml1: {
    char p[] = "Ml1"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Ml2: {
    char p[] = "Ml2"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Ml3: {
    char p[] = "Ml3"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mr1: {
    char p[] = "Mr1"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mr2: {
    char p[] = "Mr2"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mr3: {
    char p[] = "Mr3"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mq1: {
    char p[] = "Mq1"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mq2: {
    char p[] = "Mq2"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mq3: {
    char p[] = "Mq3"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mu1: {
    char p[] = "Mu1"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mu2: {
    char p[] = "Mu2"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mu3: {
    char p[] = "Mu3"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Md1: {
    char p[] = "Md1"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Md2: {
    char p[] = "Md2"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Md3: {
    char p[] = "Md3"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mh3: {
    char p[] = "MH3"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Tb: {
    char p[] = "tb";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::At: {
    char p[] = "At";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Ab: {
    char p[] = "Ab";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Al: {
    char p[] = "Al";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Au: {
    char p[] = "Au";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Ad: {
    char p[] = "Ad";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Alfsmz: {
    char p[] = "alfSMZ"; // NOLINT
    findVal(p, &val);    // NOLINT
  } break;
  case MssmParameter::Mz: {
    char p[] = "MZ";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mw: {
    char p[] = "MW";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mtp: {
    char p[] = "Mtp"; // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mbmb: {
    char p[] = "MbMb"; // NOLINT
    findVal(p, &val);  // NOLINT
  } break;
  case MssmParameter::Mcmc: {
    char p[] = "McMc"; // NOLINT
    findVal(p, &val);  // NOLINT
  } break;
  case MssmParameter::Gg: {
    char p[] = "GG";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Q: {
    char p[] = "Q";   // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Ee: {
    char p[] = "EE";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Ml: {
    char p[] = "Ml";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Mq: {
    char p[] = "Mq";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  case MssmParameter::Am: {
    char p[] = "Am";  // NOLINT
    findVal(p, &val); // NOLINT
  } break;
  default:
    break;
  }

  return val;
}

/**
 * Finds the value of variable name and assigns it to parameter val. It
 * returns a non-zero value if it cannot recognize a parameter name.
 */
double Micromegas::mssm_find_val_w(MssmParameter param) {
  switch (param) {
  case MssmParameter::Mu: {
    char p[] = "mu";
    return findValW(p);
  } break;
  case MssmParameter::Mg1: {
    char p[] = "MG1";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mg2: {
    char p[] = "MG2";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mg3: {
    char p[] = "MG3";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Ml1: {
    char p[] = "Ml1";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Ml2: {
    char p[] = "Ml2";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Ml3: {
    char p[] = "Ml3";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mr1: {
    char p[] = "Mr1";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mr2: {
    char p[] = "Mr2";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mr3: {
    char p[] = "Mr3";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mq1: {
    char p[] = "Mq1";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mq2: {
    char p[] = "Mq2";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mq3: {
    char p[] = "Mq3";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mu1: {
    char p[] = "Mu1";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mu2: {
    char p[] = "Mu2";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mu3: {
    char p[] = "Mu3";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Md1: {
    char p[] = "Md1";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Md2: {
    char p[] = "Md2";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Md3: {
    char p[] = "Md3";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mh3: {
    char p[] = "MH3";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Tb: {
    char p[] = "tb";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::At: {
    char p[] = "At";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Ab: {
    char p[] = "Ab";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Al: {
    char p[] = "Al";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Au: {
    char p[] = "Au";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Ad: {
    char p[] = "Ad";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Alfsmz: {
    char p[] = "alfSMZ"; // NOLINT
    return findValW(p);  // NOLINT
  } break;
  case MssmParameter::Mz: {
    char p[] = "MZ";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mw: {
    char p[] = "MW";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mtp: {
    char p[] = "Mtp";   // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mbmb: {
    char p[] = "MbMb";  // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mcmc: {
    char p[] = "McMc";  // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Gg: {
    char p[] = "GG";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Q: {
    char p[] = "Q";     // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Ee: {
    char p[] = "EE";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Ml: {
    char p[] = "Ml";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Mq: {
    char p[] = "Mq";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  case MssmParameter::Am: {
    char p[] = "Am";    // NOLINT
    return findValW(p); // NOLINT
  } break;
  default:
    throw std::runtime_error("Cannot get here");
  }
}

/**
 * Reads parameters from a file. The file should contain two columns with the
 * following format:
 *  name value
 *
 * Returns zero when the file has been read successfully, a negative value
 * when the file cannot be opened for reading and a positive value
 * corresponding to the line where a wrong file record was found.
 */
void Micromegas::read_var(const std::string &fname) {
  const int err = readVar(string_to_char_vec(fname).data());

  if (err == -1) {
    throw std::runtime_error("Cannot open the file " + fname + ".\n");
  }
  if (err > 0) {
    throw std::runtime_error("Wrong file contents at line " +
                             std::to_string(err) + " in file " + fname + ".\n");
  }
}

// ===========================
// ---- Utility Functions ----
// ===========================

/**
 * Sorts the odd particles with increasing masses.
 *
 * This routine fills the text parameters CDM1 and CDM2 with the names of the
 * lightest odd particle starting with one and two tildes respectively and
 * assigns the value of the mass of the lightest odd particle in each sector to
 * the global parameters Mcdm1 and Mcdm2. For models with only one DM candidate,
 * micrOMEGAs will set CDM2=NULL and Mcdm2=0. This routine returns a non zero
 * error code for a wrong set of parameters, for example parameters for which
 * some constraint cannot be calculated. The name of the corresponding
 * constraint is written in txt. This routine has to be called after a
 * reassignment of any input parameter.
 */
void Micromegas::sort_odd_particles() {
  static constexpr size_t NAME_LEN = 10;
  std::array<char, NAME_LEN> name{};
  const int err = sortOddParticles(name.data());

  if (err != 0) {
    throw std::runtime_error("Can't calculate " + std::string(name.data()));
  }
}

/**
 * The sortOddParticles command which must be used to recompute the particle
 * spectrum after changing the model parameters also clears the decay table.
 */
void Micromegas::clean_decay_table() { cleanDecayTable(); }

// ============================
// ---- Particle Functions ----
// ============================

/**
 * Returns the PDG code.
 */
int Micromegas::name_to_pdg(const std::string &name) {
  // TODO What does this return if unknown particle?
  return pNum(string_to_char_vec(name).data());
}
/**
 * Returns the numerical value of the particle mass.
 */
double Micromegas::name_to_mass(const std::string &name) {
  // TODO What does this return if unknown particle?
  return pMass(string_to_char_vec(name).data());
}

/**
 * Returns the name of the particle which PDG code is nPDG. If this particle
 * does not exist in the model the return value is NULL.
 */
std::string Micromegas::pdg_to_name(int pdg) { return pdg2name(pdg); }

// =========================
// ---- Pheno Functions ----
// =========================

/**
 * Returns the value of the supersymmetric contribution to the anomalous
 * magnetic moment of the muon.
 */
double Micromegas::mssm_gmuon() { return gmuon(); }

/**
 * Calculates the ∆ρ parameter in the MSSM. It contains for example the
 * stop/sbottom contributions, as well as the two-loop QCD corrections due to
 * gluon exchange and the correction due to gluino exchange in the heavy
 * gluino limit
 */
double Micromegas::mssm_deltarho() { return deltarho(); }

/**
 * Returns the value of the branching ratio for b → sγ, see Appendix A. We
 * have included some new contributions beyond the leading order that are
 * especially important for high tan β. SMbsg gives the SM contribution.
 */
std::pair<double, double> Micromegas::mssm_bsgnlo() {
  double smbsg = 0.0;
  const double nlo = bsgnlo(&smbsg);
  return {nlo, smbsg};
}

/**
 * Returns the value of the branching ratio Bs → μ+μ− in the MSSM. It includes
 * the loop contributions due to chargino, sneutrino, stop and Higgs exchange.
 * The ∆mb effect relevant for high tan β is taken into account.
 */
double Micromegas::mssm_bsmumu() { return bsmumu(); }

/**
 * Computes the ratio between the MSSM and SM branching fractions for  ̄B+ → τ
 * +ντ .
 */
double Micromegas::mssm_btaunu() { return btaunu(); }

/**
 * Computes the ratio of the MSSM to SM value for Rl23 in K+ → μν due to a
 * charged higgs contribution.
 */
double Micromegas::mssm_rl23() { return Rl23(); }

/**
 * Computes the branching ratio for D+s → τ +ντ . `dmunu` gives the branching
 * ratio for D+s →μ+νμ.
 */
std::pair<double, double> Micromegas::mssm_d_taunu_and_munu() {
  double muon = 0.0;
  double tau = dtaunu(&muon);
  return {tau, muon};
}

/**
 * Returns a positive value and prints a WARNING when the choice of parameters
 * conflicts with a direct accelerator limits on sparticle masses from LEP.
 * The constraint on the light Higgs mass from the LHC is included.
 */
double Micromegas::mssm_masslimits() { return masslimits(); }

// =================================
// ---- Relic Density Functions ----
// =================================

std::pair<double, double> Micromegas::relic_density(bool fast, double beps) {
  if (CDM1 == nullptr) {
    return {std::nan("0"), std::nan("0")};
  }

  double xf = 0.0;
  int err = 0;
  const int f = fast ? 1 : 0;
  const double omega = darkOmega(&xf, f, beps, &err);
  if (err != 0) {
    throw std::runtime_error(
        "Temperature where thermal equilibrium between the DM "
        "and SM sectors is too large");
  }
  return {omega, xf};
}

// =================================
// ---- Neutrino Flux Functions ----
// =================================

/**
 * Calculates the muon neutrino/anti-neutrino fluxes near the surface of the
 * Earth. Here f is the DM velocity distribution normalized such that
 *  ∫_{0}^{∞} v f (v)dv = 1.
 * The units are km/s for v and s2/km2 for f(v). For example, one can
 * use the same Maxwell function introduced for direct detection. This routine
 * implicitly depends on the WIMPSIM switch.
 */
std::pair<double, double> Micromegas::neutrino_flux(double (*vfv)(double),
                                                    int forSun) {
  double nu = 0.0;
  double nu_bar = 0.0;
  const int err = neutrinoFlux(vfv, forSun, &nu, &nu_bar);

  if (err != 0) {
    throw std::runtime_error("Error occurred in calculation of neutrino flux");
  }
  return {nu, nu_bar};
}

/**
 * calculates the exclusion confidence level for number of signal events
 * generated by given νμ and ̄νμ fluxes. The fluxes are assumed to be in [GeV
 * km2 Year]−1. This function uses the IC22BGdCos(cs) and IC22sigma(E) angular
 * distribution for background and signal as well as the event files
 * distributed by IceCube22 with φ < φcut = 8◦. The returned parameter B is a
 * Bayesian factor representing the ratio of likelihood functions for the
 * model with given fluxes and the model with null signal.
 */
double Micromegas::exclusion_confidence_icecube(double *nu, double *NU,
                                                double *B) {
  return exLevIC22(nu, NU, B);
}

// ====================================
// ---- Direct Detection Functions ----
// ====================================

/**
 * Returns
 *  F(v) = cnorm / (2pi v Rot^2)^(3/2) * exp( - v^2 / vRot^2) theta(vEsc - v)
 * which corresponds to the isothermal model. Here vRot is the orbital
 * velocity of stars in the Milky Way, it is also a global parameter of
 * micrOMEGAs. cnorm is the normalization factor.
 */
double Micromegas::maxwell(double v) { return Maxwell(v); }

/**
 * Returns
 *  F(v) = cnorm / (2pi v Rot^2)^(3/2) * exp( - v^2 / vRot^2) theta(vEsc - v)
 * which corresponds to the isothermal model. Here vRot is the orbital
 * velocity of stars in the Milky Way, it is also a global parameter of
 * micrOMEGAs. cnorm is the normalization factor.
 */
double Micromegas::shmpp(double v) { return SHMpp(v); }

/**
 * Calculates the amplitudes for CDM-nucleon elastic scattering at zero
 * momentum. pAsi(nAsi) are spin independent amplitudes for protons(neutrons)
 * whereas pAsd(nAsd) are the cor- responding spin dependent amplitudes. Each
 * of these four parameters is an array of dimension 2. The zeroth (first)
 * element of these arrays gives the χ-nucleon amplitudes whereas the second
 * element gives χ-nucleon amplitudes. Amplitudes (in GeV−2) are normalized
 * such that the total cross section for either χ or χ cross sections is:
 *  sigma_tot = 4 mx^2 mn^2 / (pi * (mx + mn)^2) * (|Asi|^2 + 3 |Asd|^2)
 */
DirectDetectionAmps Micromegas::nucleon_amplitudes(char *WINP) {
  DirectDetectionAmps results{};

  if (WINP == nullptr) {
    results.proton_si = {std::nan("0"), std::nan("0")};
    results.proton_sd = {std::nan("0"), std::nan("0")};
    results.neutron_si = {std::nan("0"), std::nan("0")};
    results.neutron_sd = {std::nan("0"), std::nan("0")};
    return results;
  }

  const int err = nucleonAmplitudes(
      WINP, results.proton_si.data(), results.proton_sd.data(),
      results.neutron_si.data(), results.neutron_sd.data());

  if (err != 0) {
    results.proton_si = {std::nan("0"), std::nan("0")};
    results.proton_sd = {std::nan("0"), std::nan("0")};
    results.neutron_si = {std::nan("0"), std::nan("0")};
    results.neutron_sd = {std::nan("0"), std::nan("0")};
  }
  return results;
}

DirectDetectionAmps Micromegas::nucleon_amplitudes_cdm1() {
  return nucleon_amplitudes(CDM1);
}

DirectDetectionAmps Micromegas::nucleon_amplitudes_cdm2() {
  return nucleon_amplitudes(CDM2);
}

/**
 * Convert string to direct-detection experiment name.
 */
static DirectDetectionExperiment
string_to_dd_experiment(const std::string &name) {
  if (name == "XENON1T_2018") {
    return DirectDetectionExperiment::Xenon1T2018;
  }

  if (name == "CRESST_2019") {
    return DirectDetectionExperiment::Cresst2019;
  }

  if (name == "DarkSide_2018") {
    return DirectDetectionExperiment::DarkSide2018;
  }

  if (name == "PICO_2019") {
    return DirectDetectionExperiment::Pico2019;
  }
  return DirectDetectionExperiment::Unknown;
}

/**
 * Returns the overall factor which should be applied to the cross sections,
 * σSIP , σSIN , σSDP , σSDN to reach the exclusion level α. All parameters
 * are the same as in DD_pvalCS.
 */
DirectDetectionResults
Micromegas::direct_detection_factor_cs(double pval, double (*vfv)(double),
                                       double cs_SI_P, double cs_SI_N,
                                       double cs_SD_P, double cs_SD_N) {
  if (CDM1 == nullptr) {
    return {std::nan("0"), std::nan("0"), std::nan("0"), std::nan("0")};
  }
  const char *expname = "";
  // clang-format off
  const double xenon = DD_factorCS(XENON1T_2018, pval, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  const double cresst = DD_factorCS(CRESST_2019, pval, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  const double darkside = DD_factorCS(DarkSide_2018, pval, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  const double pico = DD_factorCS(PICO_2019, pval, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  // clang-format off
  return {xenon, cresst, darkside, pico};
}

/**
 * Same as `direct_detection_factor_cs` but uses Maxwell distribution.
 */
DirectDetectionResults
Micromegas::direct_detection_factor_cs_maxwell(double pval, double cs_SI_P,
                                               double cs_SI_N, double cs_SD_P,
                                               double cs_SD_N) {
  return direct_detection_factor_cs(pval, Maxwell, cs_SI_P, cs_SI_N, cs_SD_P,
                                    cs_SD_N);
}

/**
 * Same as `direct_detection_factor_cs` but uses SHM++ distribution.
 */
DirectDetectionResults
Micromegas::direct_detection_factor_cs_shmpp(double pval, double cs_SI_P,
                                             double cs_SI_N, double cs_SD_P,
                                             double cs_SD_N) {
  return direct_detection_factor_cs(pval, SHMpp, cs_SI_P, cs_SI_N, cs_SD_P,
                                    cs_SD_N);
}

/**
 * calculates the value α = 1−C.L. for a model with DM-nucleon cross sections
 * σSIP , σSIN , σSDP , σSDN . Cross sections are specified in [pb] units. The
 * return value 0.1 corresponds to a 90% exclusion. The expCode parameter can
 * be any of the codes XENON1T_2018,DarkSide_2018, CRESST_2019,PICO_2019 or
 * their combination concatenated with the symbol |.
 */
DirectDetectionResults
Micromegas::direct_detection_pval_cs(double (*vfv)(double), double cs_SI_P,
                                     double cs_SI_N, double cs_SD_P,
                                     double cs_SD_N) {
  if (CDM1 == nullptr) {
    return {std::nan("0"), std::nan("0"), std::nan("0"), std::nan("0")};
  }
  const char *expname = "";
  // clang-format off
  const double xenon = DD_pvalCS(XENON1T_2018, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  const double cresst = DD_pvalCS(CRESST_2019, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  const double darkside = DD_pvalCS(DarkSide_2018, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  const double pico = DD_pvalCS(PICO_2019, vfv, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N, &expname);
  // clang-format on
  return {xenon, cresst, darkside, pico};
}
/**
 * Same as `direct_detection_pval_cs` but uses Maxwell distribution.
 */
DirectDetectionResults
Micromegas::direct_detection_pval_cs_maxwell(double cs_SI_P, double cs_SI_N,
                                             double cs_SD_P, double cs_SD_N) {
  return direct_detection_pval_cs(Maxwell, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N);
}

/**
 * Same as `direct_detection_pval_cs` but uses SMH++ distribution.
 */
DirectDetectionResults
Micromegas::direct_detection_pval_cs_shmpp(double cs_SI_P, double cs_SI_N,
                                           double cs_SD_P, double cs_SD_N) {
  return direct_detection_pval_cs(SHMpp, cs_SI_P, cs_SI_N, cs_SD_P, cs_SD_N);
}

/**
 * Similar to DD_factorCS but use the cross section calculated
 * from the DM model under consideration in micrOMEGAs. The necessary
 * corrections for a light mediator are implemented automatically, these
 * functions do not use dNdEFact.
 */
DirectDetectionResults
Micromegas::direct_detection_factor(double pval, double (*vfv)(double)) {
  if (CDM1 == nullptr) {
    return {std::nan("0"), std::nan("0"), std::nan("0"), std::nan("0")};
  }

  const char *expname = "";
  const double xenon = DD_factor(XENON1T_2018, pval, vfv, &expname);
  const double cresst = DD_factor(CRESST_2019, pval, vfv, &expname);
  const double darkside = DD_factor(DarkSide_2018, pval, vfv, &expname);
  const double pico = DD_factor(PICO_2019, pval, vfv, &expname);
  return {xenon, cresst, darkside, pico};
}

/**
 * Same as `direct_detection_factor` but uses Maxwell distribution.
 */
DirectDetectionResults
Micromegas::direct_detection_factor_maxwell(double pval) {
  return direct_detection_factor(pval, Maxwell);
}

/**
 * Same as `direct_detection_factor` but uses SMH++ distribution.
 */
DirectDetectionResults Micromegas::direct_detection_factor_shmpp(double pval) {
  return direct_detection_factor(pval, SHMpp);
}

/**
 * Similar to DD_pvalCS but use the cross section calculated
 * from the DM model under consideration in micrOMEGAs. The necessary
 * corrections for a light mediator are implemented automatically, these
 * functions do not use dNdEFact.
 */
DirectDetectionResults
Micromegas::direct_detection_pval(double (*vfv)(double)) {
  if (CDM1 == nullptr) {
    return {std::nan("0"), std::nan("0"), std::nan("0"), std::nan("0")};
  }
  const double xenon = DD_pval(XENON1T_2018, vfv, nullptr);
  const double cresst = DD_pval(CRESST_2019, vfv, nullptr);
  const double darkside = DD_pval(DarkSide_2018, vfv, nullptr);
  const double pico = DD_pval(PICO_2019, vfv, nullptr);
  return {xenon, cresst, darkside, pico};
}

/**
 * Same as `direct_detection_factor` but uses Maxwell distribution.
 */
DirectDetectionResults Micromegas::direct_detection_pval_maxwell() {
  return direct_detection_pval(Maxwell);
}

/**
 * Same as `direct_detection_factor` but uses SHM++ distribution.
 */
DirectDetectionResults Micromegas::direct_detection_pval_shmpp() {
  return direct_detection_pval(SHMpp);
}

/**
 * Returns the excluded 90% SI cross section in cm^2 from the XENON1T
 * experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.
 */
double Micromegas::xenon1T_90(double M) { return XENON1T_90(M); }

/**
 * Returns the excluded 90% SD proton cross section in cm^2 from the XENON1T
 * experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.
 */
double Micromegas::xenon1T_sdp_90(double M) { return XENON1T_SDp_90(M); }

/**
 * Returns the excluded 90% SD neutron cross section in cm^2 from the XENON1T
 * experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.
 */
double Micromegas::xenon1T_sdn_90(double M) { return XENON1T_SDn_90(M); }

/**
 * Returns the excluded 90% SI cross section in cm^2 from the DarkSide
 * experiment. For a DM mass outside 0.7 < MDM < 15 GeV, NaN is returned.
 */
double Micromegas::darkside50_90(double M) { return DS50_90(M); }

double Micromegas::darkside50_90_nob(double M) { return DS50_90_noB(M); }

/**
 * Returns the excluded 90% SI cross section in cm^2 from the CRESST III
 * experiment. For a DM mass outside 0.35 < MDM < 12 GeV, NaN is returned.
 */
double Micromegas::cresst3_90(double M) { return CRESST_III_90(M); }

/**
 * Returns the excluded 90% SD neutron cross section in cm^2 from the CRESST
 * III experiment. For a DM mass outside 0.35 < MDM < 12 GeV, NaN is returned.
 */
double Micromegas::cresst3_sdn_90(double M) { return CRESST_III_SDn_90(M); }

/**
 * Returns the excluded 90% SI cross section in cm^2 from the PICO-60
 * experiment. For a DM mass outside 3 < MDM < 10000 GeV, NaN is returned.
 */
double Micromegas::pico60_90(double M) { return PICO60_90(M); }

/**
 * Returns the excluded 90% SD proton cross section in cm^2 from the PICO-60
 * experiment. For a DM mass outside 3 < MDM < 10000 GeV, NaN is returned.
 */
double Micromegas::pico60_sdp_90(double M) { return PICO60_SDp_90(M); }

// ==============================
// ---- Constraint Functions ----
// ==============================

/**
 * Returns 1 and prints a WARNING if the invisible width of the Z boson of the
 * Standard Model is larger than 0.5 MeV and returns 0 if this
 * constraint is fulfilled. This routine can be used in any model with one DM
 * where the Z boson is uniquely defined by its PDG=23 and whether the neutral
 * LSP is its own antiparticle or not.
 */
bool Micromegas::z_invisible() {
  const int res = Zinvisible();
  return res == 1;
}

/**
 * Returns 0 if the scenario considered is not excluded by Z′ constraints, 1
 * if the point is excluded and 2 if both subroutines dealing with Z′
 * constraints cannot test the given scenario. The routine can be used for any
 * Z′ uniquely defined by the PDG code 32. Currently two types of searches
 * defined in different subroutines of Zprimelimits() are implemented.
 * The 3.2/fb Z′ search in the dilepton final state at √s = 13 TeV from ATLAS
 * is considered in the first subroutine Zprime_dilepton. If the scenario
 * considered is allowed or not tested by Zprime_dilepton, a second subroutine
 * called Zprime_dijet analyses the point using constraints from LHC dijet
 * searches at √s = 8 TeV and at √s = 13 TeV. This subroutine uses the
 * recasting performed in [88] for a combination of ATLAS and CMS searches.
 * Zprimelimits() returns 1 if MZ′ < 0.5 TeV and 2 for points for which the
 * narrow-width approximation is not valid, i.e. Γ/MZ′ > 0.3.
 */
std::pair<double, double> Micromegas::lsp_nlsp_lep() {
  double cs_out = 0.0;
  const int res = LspNlsp_LEP(&cs_out);
  return std::make_pair(res, cs_out);
}

/**
 * Checks the compatibility with the upper limit on the cross section for
 * the production of neutralinos σ(e+ e− →  ̃χ01  ̃χ0i ), i != 1, when the
 * heavier neutralino decays into quark pairs and the LSP,  ̃χ0i →  ̃χ01 + q +
 * q. The function returns σ × BR = ∑ i σ(e+e− → ̃χ01  ̃χ0i) × BR(χ0i →  ̃χ01 + q
 * + ̄q) in pb as well as a flag greater than one if σ × BR > 0.1(0.5) pb if
 * mNLSP > (<)100 GeV. This function can also be applied for non-SUSY models
 * which feature the same signature, in this case the function will compute
 * the cross section for production of the LSP and any other neutral particle
 * from the odd sector which can decay into the LSP and a Z boson.
 */
bool Micromegas::z_prime_limits() {
  const int res = Zprimelimits();
  return res == 1;
}

/**
 * Computes the cross section for p, p → pName1, pName2 + jet at √s = 8 TeV
 * where pName1, pName2 are the names of neutral outgoing particles and jet
 * includes light quarks (u,d,s) and gluons. The function returns the
 * resulting confidence level obtained with the CLs technique for each signal
 * region of the 8 TeV CMS mono-jet analysis and chooses the most
 * constraining one.
 */
double Micromegas::monojet() { return monoJet(); }

} // namespace micromegas
