#ifndef PY_MICROMEGAS_HPP
#define PY_MICROMEGAS_HPP

#include "ewsb.hpp"
#include "micromegas.hpp"
#include "results.hpp"
#include "settings.hpp"
#include "sugra.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

void define_ewsb_parameters(pybind11::class_<EwsbParameters> *ewsb);
void define_sugra_parameters(pybind11::class_<SugraParameters> *sugra);
void define_results(pybind11::class_<MicromegasResults> *results);
void define_settings(pybind11::class_<MicromegasSettings> *settings);

#endif // PY_MICROMEGAS_HPP
