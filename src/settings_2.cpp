// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Project: GeneMark.hmm-2 (no introns)
// Code was tested with boost 1.48 <boost/program_options.hpp>
// File: settings_2.cpp
// ====================================================

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <cmath>

using std::string;
using std::set;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::exception;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "settings_2.h"
#include "exit.h"

// ====================================================
InputOptions::InputOptions()
{
	shape_t  = SHAPE_TYPES::NON;
	format_t = FILE_FORMATS::NON;
	strand_t = STRAND_TYPES::NON;
}
// ----------------------------------------------------
FILE_FORMAT InputOptions::StringToFileFormat( std::string const & str )
{
	FILE_FORMAT f = FILE_FORMATS::NON;

	if	    ( !str.compare(TAGS::AUTO) )     f = FILE_FORMATS::AUTO;
	else if ( !str.compare(TAGS::FASTA) )    f = FILE_FORMATS::FASTA;
	else if ( !str.compare(TAGS::PLAIN) )    f = FILE_FORMATS::PLAIN;
	else if ( !str.compare(TAGS::GENBANK) )  f = FILE_FORMATS::GENBANK;
	else if ( !str.compare(TAGS::EMBL) )     f = FILE_FORMATS::EMBL;
	else
		Exit( "error, unexpected input sequence file format specified:", str );

	return f;
}
// ----------------------------------------------------
STRAND_TYPE InputOptions::StringToStrand( std::string const & str )
{
	STRAND_TYPE s = STRAND_TYPES::NON;

	if      ( !str.compare(TAGS::BOTH) )     s = STRAND_TYPES::BOTH;
	else if ( !str.compare(TAGS::DIRECT) )   s = STRAND_TYPES::DIRECT;
	else if ( !str.compare(TAGS::REVERSE) )  s = STRAND_TYPES::REVERSE;
	else if ( !str.compare(TAGS::DOT) )      s = STRAND_TYPES::BOTH;
	else if ( !str.compare(TAGS::PLUS) )     s = STRAND_TYPES::DIRECT;
	else if ( !str.compare(TAGS::MINUS) )    s = STRAND_TYPES::REVERSE;
	else
		Exit( "error, unexpected input strand specified:", str );

	return s;
}
// ----------------------------------------------------
SHAPE_TYPE InputOptions::StringToShape( std::string const & str )
{
	SHAPE_TYPE s = SHAPE_TYPES::NON;

	if      ( !str.compare(TAGS::INCOMPLETE) )  s = SHAPE_TYPES::INCOMPLETE;
	else if ( !str.compare(TAGS::LINEAR) )      s = SHAPE_TYPES::LINEAR;
	else if ( !str.compare(TAGS::CIRCULAR) )    s = SHAPE_TYPES::CIRCULAR;
	else
		Exit( "error, unexpected input shape specified:", str );

	return s;
}
// ----------------------------------------------------
void InputOptions::CheckAndSet(void)
{
	format_t = StringToFileFormat( format );
	strand_t = StringToStrand( strand );
	shape_t  = StringToShape( shape );
}

// ====================================================
OutputOptions::OutputOptions()
{
	format_t = FILE_FORMATS::NON;
}
// ----------------------------------------------------
FILE_FORMAT OutputOptions::StringToFileFormat(std::string const & str)
{
	FILE_FORMAT f = FILE_FORMATS::NON;

	if      ( !str.compare(TAGS::LST) )   f = FILE_FORMATS::LST;
	else if ( !str.compare(TAGS::GFF) )   f = FILE_FORMATS::GFF;
	else if ( !str.compare(TAGS::GTF) )   f = FILE_FORMATS::GTF;
	else if ( !str.compare(TAGS::GFF3) )  f = FILE_FORMATS::GFF3;
	else if ( !str.compare(TAGS::TRAIN) ) f = FILE_FORMATS::TRAIN;
	else if ( !str.compare(TAGS::EXT) )   f = FILE_FORMATS::EXT;
	else
		Exit( "error, unexpected output file format specified:", str );

	return f;
}
// ----------------------------------------------------
void OutputOptions::CheckAndSet(void)
{
	format_t = StringToFileFormat( format );
}

// ====================================================
HmmOptions::HmmOptions()
{
	bac_prob_defaulted = true;
	arc_prob_defaulted = true;
	nat_prob_defaulted = true;
	mgm_prob_defaulted = true;

	genetic_code = -1;
}
// ----------------------------------------------------
void HmmOptions::CheckModelSettings(void)
{
	if ( native_filename.empty() && mgm_filename.empty() )
		Exit( "error, at least one of the model options must be specified: -m --mod and/or -M  --Meta" );
}
// ----------------------------------------------------
void HmmOptions::CheckUpdateTransitionProb(void)
{
	if (!bac_prob_defaulted && arc_prob_defaulted)
	{
		arc_prob = 1 - bac_prob;
		arc_prob_defaulted = false;
	}
	else if (bac_prob_defaulted && !arc_prob_defaulted)
	{
		bac_prob = 1 - arc_prob;
		bac_prob_defaulted = false;
	}

	if (std::abs(bac_prob + arc_prob - 1) > 0.01)
		Exit("error on cmd: contradiction detected between --bac_prob and --arc_prob");


	if (!mgm_prob_defaulted || !nat_prob_defaulted)
	{
		if (!mgm_prob_defaulted && nat_prob_defaulted)
		{
			native_prob = 1 - mgm_prob;
			nat_prob_defaulted = false;
		}
		else if (mgm_prob_defaulted && !nat_prob_defaulted)
		{
			mgm_prob = 1 - native_prob;
			mgm_prob_defaulted = false;
		}

		if (std::abs(mgm_prob + native_prob - 1) > 0.01)
			Exit("error on cmd: contradiction detected between --mgm_prob and --native_prob");
	}
}
// ----------------------------------------------------
unsigned int HmmOptions::StringToGeneticCode(std::string const & str)
{
	unsigned int gcode = -1;

	if (!str.compare(TAGS::NO))         gcode = -1;
	else if (!str.compare(TAGS::AUTO))  gcode = 0;
	else if (!str.compare("11"))        gcode = 11;
	else if (!str.compare("4"))         gcode = 4;
	else if (!str.compare("25"))        gcode = 25;
	else if (!str.compare("1"))         gcode = 1;
	else if (!str.compare("15"))        gcode = 15;
	else
		Exit("error, unexpected genetic code specified:", str);

	return gcode;
}
// ----------------------------------------------------
void HmmOptions::CheckAndSet(void)
{
	CheckModelSettings();
	CheckUpdateTransitionProb();
	genetic_code = StringToGeneticCode(genetic_code_label);
}

// ====================================================
MGM2Options::MGM2Options()
{
	;
}
// ----------------------------------------------------
void MGM2Options::CheckAndSet(void)
{
	;
}

// ====================================================
GMS2Options::GMS2Options()
{
	;
}
// ----------------------------------------------------
void GMS2Options::CheckAndSet(void)
{
	;
}

// ====================================================
Settings::Settings(int ac, char* av[], Logger * const logger, std::string const & version) : logger(logger), version(version)
{
	ParseCmd(ac, av);
	CheckAndSet();
	InitializeLoger(logger, logfile, debug, verbose);

	if (logger->verbose)
	{
		logger->Print(0, "# Starting GeneMark.hmm ...");
		logger->Print(0, "# Parsing command line arguments ... done");
	}

	if (logger->debug)
		logger->Print(10, Summary());
}
// ----------------------------------------------------
void Settings::CheckAndSet(void)
{
	CheckFileNamesForCollision();

	in.CheckAndSet();
	out.CheckAndSet();
	hmm.CheckAndSet();

	if (gms2.gms && mgm2.mgm)
	{
		Exit("error, both --mgm2 and --gms2 options can not specified");
	}
	else
	{
		if (mgm2.mgm)
			mgm2.CheckAndSet();
		if (gms2.gms)
			gms2.CheckAndSet();
	}
}
// ----------------------------------------------------
void Settings::InitializeLoger( Logger * const ptr, std::string const & name, unsigned int const debug_level, bool const verbose )
{
	if( !name.empty() )
		ptr->Open(name);

	ptr->debug = debug_level;
	ptr->verbose = ( debug_level ) ? true : verbose;
}
// ----------------------------------------------------
void Settings::ParseCmd( int ac, char* av[] )
{
	// bool variables, which are 'TRUE' if present on command line and 'FALSE' otherwise are processed by bool_switch().

	po::options_description input_op_short( "Input" );

	input_op_short.add_options()
	( "seq,s",  po::value<string>( &in.filename )->required()->notifier( &CheckFileReadable ),
		"Name of file with sequence/s" )
	;

	po::options_description input_op_long = input_op_short;

	input_op_long.add_options()
	( "seq_format",          po::value<string>( &in.format )->default_value( TAGS::AUTO ),
		"Format of file with input sequence. Supported formats are: 'fasta', 'plain', 'genbank', 'embl' and 'auto'. In 'auto' mode algorithms will try to determine format from input file itself." )

	( "softmask",            po::bool_switch( &in.softmasked )->default_value( false, "false" ),
		"Mask lowercase letters in input sequence" )

	( "incomplete_at_gaps",  po::value<unsigned int>( &in.incomplete_at_gaps )->default_value( 40 ),
		"Allow prediction of incomplete genes at boundaries of internal gaps longer than the specified length. Gaps are defined as strings of unknown letter 'NN...N'. To disable such predictions set this value to zero" )

	( "strand",              po::value<string>( &in.strand )->default_value( TAGS::BOTH ),
		"Strand to predict genes on. Supported: 'direct', 'reverse', 'plus', 'minus', '.' or 'both'" )

	( "shape",               po::value<string>( &in.shape )->default_value( TAGS::INCOMPLETE ),
		"Input sequence organization. Supported: 'incomplete', 'circular' or 'linear'. With 'incomplete' incomplete genes are allowed at both sequence boundaries. With 'circular' a gene can span right-left boundary of circular genome. With 'linear' incomplete genes are not allowed at sequence boundaries of linear genomes." )

	( "defline_parse",       po::bool_switch( &in.parse_defline )->default_value( false, "false" ),
		"Parse sequence shape/gcode from FASTA defline. Shape must be specified as [key=value], where the value is one of the allowed shapes. Defaults to 'incomplete' if value was not provided at defline." )

	( "infofile",            po::value<string>( &in.infofile )->notifier( &CheckFileReadable ),
		"Parse sequence shape/gcode from file. File must have at least two columns with the name of the sequence in the first column and the shape in the second. Defaults to 'incomplete' if value was not provided in the file." )
	
	( "ignore", po::value<string>(&in.ignore_records)->notifier(&CheckFileReadable),
		"Ignore FASTA records listed in this file. File must have at least two columns with the name of the sequence in the first column.")

	( "run_on", po::value<string>(&in.run_on_records)->notifier(&CheckFileReadable),
		"Run algorithm only on FASTA records listed in this file. File must have at least two columns with the name of the sequence in the first column.")
	;

	po::options_description output_op_short( "Output" );

	output_op_short.add_options()
	( "out,o",     po::value<string>( &out.filename )->required(),
		"Output file" )

	( "format,f",  po::value<string>( &out.format )->default_value( TAGS::LST ),
		"Output format; supported: lst, gff, gtf, gff3, train and ext" )
	;

	po::options_description output_op_long = output_op_short;

	output_op_long.add_options()
	( "AA",  po::value<string>( &out.aa_filename ),
		"Output protein sequences of predicted genes to this file" )

	( "NT",  po::value<string>( &out.nt_filename ),
		"Output nucleotide sequences of predicted genes to this file" )

	( "aa",  po::bool_switch( &out.add_aa )->default_value( false, "false" ),
		"Add protein sequences of predicted genes to output file" )

	( "nt",  po::bool_switch( &out.add_nt )->default_value( false, "false" ),
		"Add nucleotide sequences of predicted genes to output file" )

	( "gid_per_contig", po::bool_switch(&out.gid_per_contig)->default_value(false, "false"),
		"Construct gene ID from contig ID and gene order in the contig. By default use gene order in the file as ID.")
	
	( "gid_start", po::value<int>(&out.gid_start)->default_value(1, "1"),
		"Start gene ID with this number.")

	( "gid_label", po::value<string>(&out.gid_label)->default_value("", "no label"),
		"Add this label as prefix to gene IDs.")

	("standard", po::bool_switch(&out.standard_partial)->default_value(false, "false"),
		"Format output files according to strict standard specification")
	;

	po::options_description plus_op_short( "Plus" );
	
	plus_op_short.add_options()
	( "evi,e",  po::value<string>( &in.evidence )->notifier( &CheckFileReadable ),
		"File with evidence data in GFF format for PLUS prediction algorithm" )
	;

	po::options_description plus_op_long = plus_op_short;

	po::options_description algorithm_op_short( "Algorithm parameters" );

	algorithm_op_short.add_options()
	( "mod,m",   po::value<string>( &hmm.native_filename )->notifier( &CheckFileReadable ),
		"File with species specific parameters" )

	( "Meta,M",  po::value<string>( &hmm.mgm_filename )->notifier( &CheckFileReadable ),
		"File with MetaGeneMark parameters" )
	;

	po::options_description algorithm_op_long = algorithm_op_short;

	algorithm_op_long.add_options()
	( "tis_par",   po::value<string>( &hmm.tis_filename )->notifier( &CheckFileReadable ),
		"Use TIS parameters from this file" )
	;

	po::options_description other_op_short( "Other options" );

	other_op_short.add_options()
	( "help,h",
		"Full help message" )

	( "verbose,v", po::bool_switch( &verbose )->default_value( false, "false" ),
		"Verbose" )
	;

	po::options_description other_op_long = other_op_short;

	other_op_long.add_options()
	( "pbar",     po::bool_switch( &progress_bar )->default_value( false, "false" ),
		"Show progress bar" )

	( "usr_cfg",  po::value<string>( &user_config )->notifier( &CheckFileReadable ),
		"Read program option values from this file" )
	;

	po::options_description developer_op( "Developer options" );

	developer_op.add_options()
	( "all",
		"Full developer help message" )

	( "debug",    po::value<unsigned int>( &debug )->default_value( 0, "no" ),
			"Debug level; zero - no debug" )

	( "logfile",  po::value<string>( &logfile ),
		"Logger" )

	// configuration file is placed in the same folder as the executable
	// not clear how to find this path portable way in C++
	// ->default_value("gmhmmp2.cfg")

	( "cfg",  po::value<string>( &config )->notifier(&CheckFileReadable),
		"Configuration file" )

	// Algo

	( "file_load_size",    po::value<unsigned int>( &hmm.file_load_size )->default_value( 20000000, "20,000,000" ),
		"Load sequence files below this size into memory at once" )

	( "ini_non_prob",     po::value<double>( &hmm.noncoding_state_initiation_and_termination_probability )->default_value( 0.5 )->notifier( &CheckProbabilityRange ),
		"Initiation and termination probability of noncoding state" )

	( "bac_prob",         po::value<double>( &hmm.bac_prob )->default_value( 0.5 )->notifier(&CheckProbabilityRange),
		"Transition probability to Atypical bacteria" )

	( "arc_prob",         po::value<double>( &hmm.arc_prob )->default_value( 0.5 )->notifier(&CheckProbabilityRange),
		"Transition probability to Atypical archaea" )

	( "mgm_prob",         po::value<double>( &hmm.mgm_prob )->default_value( 1 )->notifier(&CheckProbabilityRange),
		"Transition probability to Atypical state" )

	( "native_prob",      po::value<double>( &hmm.native_prob )->default_value( 1 )->notifier(&CheckProbabilityRange),
		"Transition probability to Native state" )

	( "overlap_penalty",  po::value<double>( &hmm.overlap_penalty )->default_value( 0.5 )->notifier( &CheckProbabilityRange ),
		"Overlap penalty" )

	( "allow_overlap",    po::value<bool>( &hmm.allow_gene_overlap )->default_value( true, "true" ),
		"Allow gene overlap; supported 'on' and 'off' or 'true' and 'false'. Turn gene overlap OFF in intronless eukaryotes." )

	( "allow_tis",        po::value<bool>( &hmm.allow_start_models )->default_value( true, "true" ),
		"Allow TIS models; 'on' and 'off' or 'true' and 'false'" )
		
	( "incomplete_dur",   po::value<bool>( &hmm.allow_incomplete_gene_duration )->default_value( true, "true" ),
		"Allow incomplete duration for incomplete states" )

	( "best_start_before_dp",  po::value<bool>( &hmm.best_start_before_dp )->default_value( false, "false" ),
		"Select best ORF start before DP" )

	( "use_this_gc_mgm",       po::value<unsigned int>( &hmm.use_this_gc_mgm)->default_value( -1 , "use all" ),
		"Use only this GC model from MGM set" )
			
	("mgm2", po::value<bool>(&mgm2.mgm)->default_value(false, "false"),
		"Run in MetaGeneMark-2 mode, under construction")

	("gms2", po::value<bool>(&gms2.gms)->default_value(false, "false"),
		"Run in GeneMarkS-2 mode, under construction")

	( "gcode,g",  po::value<std::string>( &hmm.genetic_code_label )->default_value( TAGS::NO ),
		"Genetic code of input sequence. Supported: 11, 4, 25, 15 and 1. Value 'auto' is reserved for auto-detection." )

	( "species",  po::value<string>( &gms2.species_name )->default_value( "undefined" ),
		"The name of the species to output in the species specific parameter file" )

	("delta", po::value<double>(&hmm.delta)->default_value(0),
		"gene prediction threshold in logood space.")
	;

	// join specifications

	po::options_description short_op("");
	short_op.add( input_op_short ).add(output_op_short).add(plus_op_short).add(algorithm_op_short).add(other_op_short);

	po::options_description long_op("");
	long_op.add( input_op_long ).add(output_op_long).add(plus_op_long).add(algorithm_op_long).add(other_op_long);

	po::options_description all_op("");
	all_op.add( input_op_long ).add(output_op_long).add(plus_op_long).add(algorithm_op_long).add(other_op_long).add(developer_op);

	// --

	po::variables_map vm;

	try
	{
		// parse command line

		po::store( po::parse_command_line(ac, av, all_op), vm );

		// previous call has higher priority than next: cmd parameters overwrite parameters from user_cfg

		// error in boost documentation description "boost parse_config_file"
		// this works on WINDOWS: po::store( po::parse_config_file<char>( user_config.c_str(), long_op), vm );
		// fails on LINUX
		// pass file stream instead of file name; works on both tested platforms

		if ( !user_config.empty() )
		{
			ifstream file( user_config.c_str() );
			CheckFileReadable( user_config );
			po::store( po::parse_config_file( file , all_op ), vm );
		}

		// previous calls have higher priority than next: 'cmd' and 'user_cfg' overwrites 'config'

		if ( !config.empty() )
		{
			ifstream file( config.c_str() );
			CheckFileReadable( config );
			po::store( po::parse_config_file( file , all_op ), vm );
		}

		// Usage and EXIT

		if      ( ac == 1 )          Exit( UsageHeader(av[0]), short_op );
		else if ( vm.count("help") ) Exit( UsageHeader(av[0]), long_op );
		else if ( vm.count("all") )  Exit( UsageHeader(av[0]), all_op );

		// set and check values

		hmm.bac_prob_defaulted = vm["bac_prob"].defaulted();
		hmm.arc_prob_defaulted = vm["arc_prob"].defaulted();
		hmm.nat_prob_defaulted = vm["native_prob"].defaulted();
		hmm.mgm_prob_defaulted = vm["mgm_prob"].defaulted();

		po::notify(vm);
	}
	catch( exception & e )
	{
		Exit( "error:", e.what() );
	}
}
// ----------------------------------------------------
string Settings::UsageHeader( char* name )
{
	stringstream stst;

	stst << endl;
	stst << "GeneMark.hmm-2 version " << version << endl;
	stst << "Usage: " << name << " --option value" << endl;

	return stst.str();
}
// ----------------------------------------------------
string Settings::Summary(void)
{
	stringstream s;

	s << "# START Settings" << endl;

	s << "version : " << version << endl;
	s << "logger : " << logger << endl;

	s << "verbose : " << verbose << endl;
	s << "debug : " << debug << endl;
	s << "logfile : " << logfile << endl;
	s << "progress_bar : " << progress_bar << endl;
	s << "config : " << config << endl;
	s << "user_config : " << user_config << endl;

	s << "in.filename : " << in.filename << endl;
	s << "in.format : " << in.format << endl;
	s << "in.shape : " << in.shape << endl;
	s << "in.incomplete_at_gaps : " << in.incomplete_at_gaps << endl;
	s << "in.softmasked : " << in.softmasked << endl;
	s << "in.strand : " << in.strand << endl;
	s << "in.parse_defline : " << in.parse_defline << endl;
	s << "in.infofile : " << in.infofile << endl;
	s << "in.evidence : " << in.evidence << endl;
	s << "in.format_t : " << in.format_t << endl;
	s << "in.strand_t : " << in.strand_t << endl;
	s << "in.shape_t : " << in.shape_t << endl;
	s << "in.ignore_records : " << in.ignore_records << endl;
	s << "in.run_on_records : " << in.run_on_records << endl;

	s << "gms2.species_name : " << gms2.gms << endl;
	s << "gms2.species_name : " << gms2.species_name << endl;
	
	s << "mgm2.mgm : " << mgm2.mgm << endl;

	s << "out.filename : " << out.filename << endl;
	s << "out.format : " << out.format << endl;
	s << "out.add_aa : " << out.add_aa << endl;
	s << "out.add_nt : " << out.add_nt << endl;
	s << "out.aa_filename : " << out.aa_filename << endl;
	s << "out.nt_filename : " << out.nt_filename << endl;
	s << "out.format_t : " << out.format_t << endl;
	s << "out.gid_per_contig : " << out.gid_per_contig << endl;
	s << "out.gid_start : " << out.gid_start << endl;
	s << "out.gid_label : " << out.gid_label << endl;

	s << "hmm.native_filename : " << hmm.native_filename << endl;
	s << "hmm.mgm_filename : " << hmm.mgm_filename << endl;
	s << "hmm.tis_filename : " << hmm.tis_filename << endl;
	s << "hmm.bac_prob : " << hmm.bac_prob << endl;
	s << "hmm.arc_prob : " << hmm.arc_prob << endl;
	s << "hmm.native_prob : " << hmm.native_prob << endl;
	s << "hmm.mgm_prob : " << hmm.mgm_prob << endl;
	s << "hmm.file_load_size : " << hmm.file_load_size << endl;
	s << "hmm.noncoding_state_initiation_and_termination_probability : " << hmm.noncoding_state_initiation_and_termination_probability << endl;
	s << "hmm.overlap_penalty : " << hmm.overlap_penalty << endl;
	s << "hmm.allow_gene_overlap : " << hmm.allow_gene_overlap << endl;
	s << "hmm.allow_start_models : " << hmm.allow_start_models << endl;
	s << "hmm.allow_incomplete_gene_duration : " << hmm.allow_incomplete_gene_duration << endl;
	s << "hmm.best_start_before_dp : " << hmm.best_start_before_dp<< endl;
	s << "hmm.bac_prob_defaulted : " << hmm.bac_prob_defaulted << endl;
	s << "hmm.arc_prob_defaulted : " << hmm.arc_prob_defaulted << endl;
	s << "hmm.nat_prob_defaulted : " << hmm.nat_prob_defaulted << endl;
	s << "hmm.mgm_prob_defaulted : " << hmm.mgm_prob_defaulted << endl;
	s << "hmm.use_this_gc_mgm : " << hmm.use_this_gc_mgm << endl;
	s << "hmm.genetic_code_label : " << hmm.genetic_code_label << endl;
	s << "hmm.genetic_code : " << hmm.genetic_code << endl;

	s << "# END Settings" << endl;

	return s.str();
}
// ----------------------------------------------------
void Settings::CheckFileNamesForCollision(void)
{
	set<string> unique_names;

	// add input file names to set

	unique_names.insert( in.filename );
	unique_names.insert( in.infofile );
	unique_names.insert( in.evidence );
	unique_names.insert( hmm.native_filename );
	unique_names.insert( hmm.mgm_filename );
	unique_names.insert( hmm.tis_filename );
	unique_names.insert( user_config );
	unique_names.insert( config );

	// check that output files are unique

	AddNameToSet( out.filename,    unique_names );
	AddNameToSet( out.aa_filename, unique_names );
	AddNameToSet( out.nt_filename, unique_names );
	AddNameToSet( logfile,         unique_names );
}
// ----------------------------------------------------
void Settings::AddNameToSet( std::string const & name, std::set<std::string> & unique_names )
{
	if ( name.empty() ) return;

	if ( unique_names.find( name ) == unique_names.end() )
		unique_names.insert( name );
	else
		Exit( "error, file name collision detected for:", name );
}
// ----------------------------------------------------
void Settings::CheckFileReadable(std::string const & name)
{
	if ( name.empty() )
		Exit("error, file name is empty in Settings::CheckFileReadable");

	ifstream file( name.c_str() );

	if ( ! file.is_open() )
		Exit( "error on open file:", name );
}
// ----------------------------------------------------
void Settings::CheckProbabilityRange(double const p)
{
	if ( p < 0  || p > 1 )
		Exit( "error, parameter value is out of allowed range 0..1:", p );
}
// ----------------------------------------------------

