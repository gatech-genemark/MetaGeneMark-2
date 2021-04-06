// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: pset_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <cassert>
#include <cmath>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

using std::string;
using std::vector;
using std::map;
using std::cerr;
using std::endl;
using std::cout;

#include "pset_2.h"
#include "exit.h"

Pset::Pset() {
    // empty constructor ; should only be used by copy constructor. Never call it directly
}
// ----------------------------------------------------
Pset::Pset( Settings & settings, Parameters & parameters, Logger * logger ) : logger(logger)
{
	Model* value = 0;

	native.assign(101,value);
	first.assign(101,value);
	second.assign(101,value);

	LinkModelSetToModels( parameters.model_set, native, "__NATIVE" );
	LinkModelSetToModels( parameters.model_set, first,  "__B_GC_" );
	LinkModelSetToModels( parameters.model_set, second, "__A_GC_" );

	FillGapsInModelSet(native);
	FillGapsInModelSet(first);
	FillGapsInModelSet(second);

	genetic_code    = GetGeneticCode();
	min_gene_length = GetMinGeneLength();
	
	SetStarts( native, first );
	SetStarts( native, second );

	SetStops( native, first );
	SetStops( native, second );

	UpdateToModel( native, parameters.to_native );
	UpdateToModel( first,  parameters.to_atypical_first  * parameters.to_mgm );
	UpdateToModel( second, parameters.to_atypical_second * parameters.to_mgm );
}
// ----------------------------------------------------
unsigned int Pset::GetGeneticCode(void)
{
	unsigned int code = 0;

	code = GetGeneticCode( code, native );
	code = GetGeneticCode( code, first );
	code = GetGeneticCode( code, second );

	return code;
}
// ----------------------------------------------------
unsigned int Pset::GetMinGeneLength(void)
{
	unsigned int min = -1;

	min = GetMinGeneLength( min, native );
	min = GetMinGeneLength( min, first );
	min = GetMinGeneLength( min, second );

	if ( min == -1 )
		Exit( "error, unexpected value found in 'Pset::GetMinGeneLength'" );

	return min;
}
// ----------------------------------------------------
unsigned int Pset::GetGeneticCode( unsigned int code, std::vector< Model* > const & source )
{
	vector< Model* >::const_iterator  itr     = source.begin();
	vector< Model* >::const_iterator  itr_end = source.end();

	for( ; itr != itr_end; ++itr )
	{
		if ( ! *itr ) continue;

		if ( !code )
			code = (*itr)->gcode;

		if ( code != (*itr)->gcode )
			Exit( "error, inconsistent genetic code was detected:", code, (*itr)->gcode );
	}

	return code;
}
// ----------------------------------------------------
unsigned int Pset::GetMinGeneLength( unsigned int min, const std::vector< Model* > & source )
{
	vector< Model* >::const_iterator  itr     = source.begin();
	vector< Model* >::const_iterator  itr_end = source.end();

	for( ; itr != itr_end; ++itr )
	{
		if ( ! *itr ) continue;
		
		if ( min > (*itr)->gene_min_length )
			min = (*itr)->gene_min_length;
	}

	return min;
}
// ----------------------------------------------------
void Pset::LinkModelSetToModels( std::map< std::string, Model > const & source, std::vector< Model* > & target , std::string const label )
{
	map< string, Model >::const_iterator  itr     = source.begin();
	map< string, Model >::const_iterator  itr_end = source.end();

	for( ; itr != itr_end; ++itr )
	{
		string name = itr->first;

		// find cases when label matches the name

		if ( !name.compare( 0, label.size(), label ) )
		{
			// if exact match, then source with target at [0]

			if ( !name.compare(label) )
			{
				if ( target[0] )
					Exit( "error, model duplication was detected for:", itr->first );

				target[0] = (Model*) & itr->second;
			}
			else
			{
				// if not match, then find GC value and link source with target at [GC]

				name.erase( 0, label.size() );

				int i;
				std::istringstream( name ) >> i;

				if ( i < 0 || i > 100 )
					Exit( "error, unexpected value found in label:", label, i );

				if ( target[i] )
					Exit( "error, model duplication was detected for:", itr->first );

				// link the model with corresponding GC bin

				target[i] = (Model*) & itr->second;
			}
		}
	}
}
// ----------------------------------------------------
void Pset::FillGapsInModelSet(std::vector< Model* > & vec)
{
	if  ( vec.size() != 101 )
		Exit( "error, unexpected size found in 'Pset::FillGapsInModelSet'" );

	// change "0" to the non-zero following the rule:
	// 0 0 0 0 A 0 0 0 0 B 0 0 0 0 C 0 0 0
	// A A A A A B B B B B C C C C C 0 0 0

	for( int i = 0; i < 101; ++i )
	{
		if ( vec[i] != 0 )
		{
			for( int j = i - 1; j >= 0 && vec[j] == 0; --j )
			{
					vec[j] = vec[i];
			}
		}
	}

	// fill the end
	// 0 0 0 0 A 0 0 0 0 B 0 0 0 0 C 0 0 0 0
	// A A A A A B B B B B C C C C C C C C C

	for( int i = 100; i >= 0; --i )
	{
		if ( vec[i] != 0 )
		{
			for( int j = i + 1; j < 101; ++j )
			{
				vec[j] = vec[i];
			}

			break;
		}
	}
}
// ----------------------------------------------------
void Pset::SetStarts( std::vector< Model* > & source, std::vector< Model* > & target )
{
	if ( source.at(0) && target.at(0) )
	{
		if ( source.size() != target.size() )
			Exit( "error, size mismatch in 'Pset::SetStarts'" );

		for( unsigned int i = 0; i < source.size(); ++i )
		{
			target[i]->pATG = source[i]->pATG;
			target[i]->pGTG = source[i]->pGTG;
			target[i]->pTTG = source[i]->pTTG;
		}
	}
}
// ----------------------------------------------------
void Pset::SetStops( std::vector< Model* > & source, std::vector< Model* > & target )
{
	if ( source.at(0) && target.at(0) )
	{
		if ( source.size() != target.size() )
			Exit( "error, size mismatch in 'Pset::SetStops'" );

		for( unsigned int i = 0; i < source.size(); ++i )
		{
			target[i]->pTAA = source[i]->pTAA;
			target[i]->pTAG = source[i]->pTAG;
			target[i]->pTGA = source[i]->pTGA;
		}
	}
}
// ----------------------------------------------------
void Pset::UpdateToModelByGC( std::vector< Model* > & model_set, double to_set_prob, std::vector<double> & prob_by_gc )
{
	if ( ! model_set.at(0) )
		return;

	if ( to_set_prob < 0 )
		Exit( "error, p < 0 'UpdateToModelByGC'" );

	if ( model_set.size() != prob_by_gc.size() )
		Exit( "error, size mismatch in 'UpdateToModelByGC'" );

	for( unsigned int i = 0; i < prob_by_gc.size(); ++i )
	{
		if ( prob_by_gc[i] < 0 || prob_by_gc[i] > 1 )
			Exit( "error, probability out of 0..1 range in 'UpdateToModelByGC'" );
	}

	for( unsigned int i = 0; i < model_set.size(); ++i )
	{
		model_set[i]->toModel = to_set_prob * (1 - prob_by_gc[i]);
	}
}
// ----------------------------------------------------
void Pset::UpdateToModel( std::vector< Model* > & model_set, double to_set_prob )
{
	if ( ! model_set.at(0) )
		return;

	if ( to_set_prob < 0 )
		Exit( "error, p < 0 'UpdateToModel'" );

	for( vector< Model* >::iterator itr = model_set.begin(); itr != model_set.end(); ++itr )
	{
		(*itr)->toModel = to_set_prob;
	}
}
// ----------------------------------------------------

void deep_copy_models_from_vector(const vector< Model* > & original, vector< Model *> &destination) {
    destination.resize(original.size());
    
    for (size_t i = 0; i < original.size(); i++) {
        if (original[i] != NULL) {
            destination[i] = original[i]; // new Model(*(original[i]));
        }
    }
}

Pset Pset::deepCopy() const {
    
    Pset newpset;
    
    deep_copy_models_from_vector(first, newpset.first);
    deep_copy_models_from_vector(native, newpset.native);
    deep_copy_models_from_vector(second, newpset.second);
    
    newpset.logger = logger;
    newpset.genetic_code = genetic_code;
    newpset.min_gene_length = min_gene_length;
    
    return newpset;
}
