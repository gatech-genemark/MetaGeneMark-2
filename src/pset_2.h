// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: pset_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#pragma once

#ifndef PSET_2_H
#define PSET_2_H

#include <string>
#include <vector>
#include <map>

#include "model_2.h"
#include "settings_2.h"
#include "parameters_2.h"
#include "common_2.h"

// ----------------------------------------------------
class Pset
{
public:

    Pset();     // Never call this directly. It is only used by the copy constructor
	Pset( Settings & settings, Parameters & parameters, Logger * logger );
	~Pset(){};
    
    Pset deepCopy() const;

	std::vector< Model* > native;
	std::vector< Model* > first;
	std::vector< Model* > second;

	// check model files for the same setting of the genetic code

	unsigned int genetic_code;

	// different model files may have different minimum gene length
	// this is the shortest minium gene length

	unsigned int min_gene_length;

private:

	void LinkModelSetToModels( std::map< std::string, Model > const & source, std::vector< Model*> & target, std::string const label );
	void FillGapsInModelSet(std::vector< Model* > & vec);

	unsigned int GetMinGeneLength(void);
	unsigned int GetMinGeneLength( unsigned int min, const std::vector< Model* > & source );

	unsigned int GetGeneticCode(void);
	unsigned int GetGeneticCode(unsigned int code, std::vector< Model* > const & source);

	void SetStarts( std::vector< Model* > & source, std::vector< Model* > & target );
	void SetStops ( std::vector< Model* > & source, std::vector< Model* > & target );

	void UpdateToModel( std::vector< Model* > & model_set, double to_set_prob );
	void UpdateToModelByGC( std::vector< Model* > & model_set, double to_set_prob, std::vector<double> & prob_by_gc );

	Logger const * logger;
};
// ----------------------------------------------------
#endif // PSET_2_H

