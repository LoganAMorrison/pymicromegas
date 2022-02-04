#include "../../include/micromegas.h"
#include "../../include/micromegas_aux.h"
#include "../lib/pmodel.h"
#include "micromegas.hpp"

namespace micromegas {

/**
 * Calculates the masses of Higgs and supersymmetric particles in the MSSM
 * including one-loop corrections starting from weak scale input parameters
 * using `SoftSusy`.
 */
int Micromegas::mssm_ewsb() { return softSusyEwsbMSSM(); }

/**
 * Derives the parameters at the electroweak symmetry breaking scale assuming
 * thate all input parameters except `tb` and `signMu` are defined at the GUT
 * scale. Uses `SoftSusy`.
 */
int Micromegas::mssm_sugra(double tb, double gMG1, double gMG2, double gMG3,
                           double gAl, double gAt, double gAb, double sgn,
                           double gMHu, double gMHd, double gMl1, double gMl3,
                           double gMr1, double gMr3, double gMq1, double gMq3,
                           double gMu1, double gMu3, double gMd1, double gMd3) {
  return softSusySUGRA(tb, gMG1, gMG2, gMG3, gAl, gAt, gAb, sgn, gMHu, gMHd,
                       gMl1, gMl3, gMr1, gMr3, gMq1, gMq3, gMu1, gMu3, gMd1,
                       gMd3);
}
} // namespace micromegas
