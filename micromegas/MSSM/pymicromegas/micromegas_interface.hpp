#ifndef PY_MICROMEGAS_HPP
#define PY_MICROMEGAS_HPP

#include "ewsb.hpp"
#include "sugra.hpp"
#include <pybind11/pybind11.h>

namespace pymicromegas {

// /**
//  * Defines the EWSB class for the input module.
//  */
// void define_ewsb_parameters(pybind11::class_<EwsbParameters> *ewsb);

// /**
//  * Defines the SUGRA class for the input module.
//  */
// void define_sugra_parameters(pybind11::class_<SugraParameters> *sugra);

// /**
//  * Defines the interface to pheno.
//  */
// void define_spheno(pybind11::module_ *m);

} // namespace pymicromegas

/**
 * Defines the all the function for getting and setting parameters.
 */
void define_micromegas_get_set(pybind11::module_ *ll);

/**
 * Defines the direct-detection functions.
 */
void define_micromegas_direct_detection(pybind11::module_ *ll);

/**
 * Defines the pheno functions.
 */
void define_micromegas_pheno(pybind11::module_ *ll);

/**
 * Defines the utility functions.
 */
void define_micromegas_utils(pybind11::module_ *ll);

#endif // PY_MICROMEGAS_HPP
