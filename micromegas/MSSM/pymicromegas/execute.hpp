#ifndef PY_MICROMEGAS_EXECUTE_HPP
#define PY_MICROMEGAS_EXECUTE_HPP

#include "ewsb.hpp"
#include "results.hpp"
#include "settings.hpp"
#include "sugra.hpp"

void execute(MicromegasResults *results, const MicromegasSettings &settings);

void execute(MicromegasResults *results, const MicromegasSettings &settings,
             const SugraParameters &sugra);

void execute(MicromegasResults *results, const MicromegasSettings &settings,
             const EwsbParameters &ewsb);

MicromegasResults execute(const MicromegasSettings &settings,
                          const SugraParameters &sugra);

MicromegasResults execute(const MicromegasSettings &settings,
                          const EwsbParameters &ewsb);

MicromegasResults execute(const MicromegasSettings &settings,
                          const std::vector<SugraParameters> &sugras);

MicromegasResults execute(const MicromegasSettings &settings,
                          const std::vector<EwsbParameters> &ewsbs);

#endif // PY_MICROMEGAS_EXECUTE_HPP
