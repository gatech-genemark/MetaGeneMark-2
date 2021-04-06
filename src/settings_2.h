// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Project: GeneMark.hmm-2 (no introns)
// Code was tested with boost 1.48 <boost/program_options.hpp>
// File: settings_2.h
// ====================================================

#ifndef GMHMMP2_SETTINGS_2_H
#define GMHMMP2_SETTINGS_2_H

#include <string>
#include <set>

#include "common_2.h"
#include "logger.h"

// ----------------------------------------------------
class InputOptions
{
public:
	InputOptions();
	~InputOptions(){}

	void CheckAndSet(void);

	std::string   filename;
	std::string   format;
	std::string   shape;
	std::string   strand;
	std::string   infofile;
	std::string   evidence;
	bool          softmasked;
	bool          parse_defline;
	unsigned int  incomplete_at_gaps;
	std::string   ignore_records;
	std::string   run_on_records;

	SHAPE_TYPE    shape_t;
	FILE_FORMAT   format_t;
	STRAND_TYPE   strand_t;

private:
	SHAPE_TYPE    StringToShape( std::string const & str );
	FILE_FORMAT   StringToFileFormat( std::string const & str );
	STRAND_TYPE   StringToStrand( std::string const & str );
};
// ----------------------------------------------------
class OutputOptions
{
public:
	OutputOptions();
	~OutputOptions(){}

	void CheckAndSet(void);

	std::string  filename;
	std::string  format;
	std::string  aa_filename;
	std::string  nt_filename;
	bool         add_aa;
	bool         add_nt;
	bool         gid_per_contig;
	int          gid_start;
	std::string  gid_label;
	bool         standard_partial;   // not implemented

	FILE_FORMAT  format_t;

private:
	FILE_FORMAT  StringToFileFormat( std::string const & str );
};
// ----------------------------------------------------
class HmmOptions
{
public:
	HmmOptions();
	~HmmOptions(){}

	void CheckAndSet(void);

	std::string  native_filename;
	std::string  mgm_filename;
	std::string  tis_filename;
	unsigned int use_this_gc_mgm;   // not implemented
	
	double       bac_prob;
	double       arc_prob;
	double       native_prob;
	double       mgm_prob;

	unsigned int file_load_size;
	double       noncoding_state_initiation_and_termination_probability;
	double       overlap_penalty;
	bool         allow_gene_overlap;   // not implemented
	bool         allow_start_models;   // not implemented
	bool         allow_incomplete_gene_duration;   // not implemented
	bool         best_start_before_dp;
	std::string  genetic_code_label;   // not implemented

	// "defaulted" is "true" if value was not modified due to command line options

	bool  bac_prob_defaulted;
	bool  arc_prob_defaulted;
	bool  nat_prob_defaulted;
	bool  mgm_prob_defaulted;

	unsigned int genetic_code;   // not implemented

	double       delta; 

private:
	unsigned int StringToGeneticCode(std::string const & str);
	void CheckUpdateTransitionProb(void);
	void CheckModelSettings(void);
};
// ----------------------------------------------------
class MGM2Options
{
public:
	MGM2Options();
	~MGM2Options() {}

	void CheckAndSet(void);

	bool         mgm;
};
// ----------------------------------------------------
class GMS2Options
{
public:
	GMS2Options();
	~GMS2Options(){}

	void CheckAndSet(void);

	bool         gms;
	std::string  species_name;
};
// ----------------------------------------------------
class Settings
{
public:
	Settings( int ac, char* av[], Logger * const logger, std::string const & version );
	~Settings(){}

	InputOptions  in;
	OutputOptions out;
	HmmOptions    hmm;
	GMS2Options   gms2;
	MGM2Options   mgm2;

	bool          verbose;
	int unsigned  debug;
	std::string   logfile;
	bool          progress_bar;

	std::string   config;
	std::string   user_config;

private:

	Logger * const logger;
	std::string const version;

	void ParseCmd( int ac, char* av[] );
	void CheckAndSet(void);
	void InitializeLoger( Logger * const ptr, std::string const & name, unsigned int const debug_level, bool const verbose );

	std::string Summary(void);

	void CheckFileNamesForCollision(void);
	void AddNameToSet( std::string const & name, std::set<std::string> & unique_names );
	std::string UsageHeader( char* name );

	// static for call with ->notifier(&function)

	static void CheckFileReadable(std::string const & name);
	static void CheckProbabilityRange(double const p);
};
// ----------------------------------------------------
#endif // GMHMMP2_SETTINGS_2_H

