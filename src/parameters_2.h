// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: parameters_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

// to do:  in 'Load' MGM key "__" is used;  keys "__B_GC" and "__A_GC" should be used;

#ifndef PARAMETERS_2_H
#define PARAMETERS_2_H

#include <string>
#include <map>
#include <vector>

#include "logger.h"
#include "model_2.h"
#include "settings_2.h"
#include "parameter_parser_2.h"
#include "multi_shift_site.hpp"

// ----------------------------------------------------
class Parameters
{
public:
	Parameters( Settings const & settings, Logger * const logger );
	~Parameters(){};

	std::map< std::string, Model > model_set;

	double to_atypical_first;
	double to_atypical_second;
	double to_native;
	double to_mgm;

	// not used yet; untested consept: native should have close GC to native
	std::vector<double> to_native_by_gc;

	void Initialize( std::map< std::string, Model > & target );

private:

	void Load(void);
	
	void LoadSection( ParameterParser & parser, std::string & buffer, std::string const label );
	void LoadSite( ParameterParser & parser, parameter_map & par, std::string & buffer, std::string const label, Site * ptr, bool with_dur );
    void LoadMultiShiftSite( ParameterParser & parser, parameter_map & par, std::string & buffer, std::string const label, MultiShiftSite * ptr, bool with_dur );

	std::string Summary(void);

	// settings

	std::string  native_filename;
	std::string  mgm_filename;
	std::string  tis_filename;

	Logger * const logger;
};
// ----------------------------------------------------
#endif // PARAMETERS_2_H

