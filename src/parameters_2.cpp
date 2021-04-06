// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: parameters_2.cpp
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

using std::string;
using std::map;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::pair;

#include "parameters_2.h"
#include "parameter_parser_2.h"
#include "exit.h"

// ----------------------------------------------------
Parameters::Parameters( Settings const & settings, Logger * const logger ) : logger(logger)
{
	if ( logger->verbose ) logger->Print( 0, "# Loading parameters ..." );

	native_filename = settings.hmm.native_filename;
	mgm_filename    = settings.hmm.mgm_filename;
	tis_filename    = settings.hmm.tis_filename;

	// the same parameter can be specified in three locations
	// location processing order: defaults < from model file < from command line

	to_native = 0.5;
	to_mgm = 0.5;
	to_atypical_first  = 0.5;
	to_atypical_second = 0.5;

	// load from file
	Load();

	// load from command line
	if ( ! settings.hmm.bac_prob_defaulted )  to_atypical_first  = settings.hmm.bac_prob;
	if ( ! settings.hmm.arc_prob_defaulted )  to_atypical_second = settings.hmm.arc_prob;
	if ( ! settings.hmm.nat_prob_defaulted )  to_native          = settings.hmm.native_prob;
	if ( ! settings.hmm.mgm_prob_defaulted )  to_mgm             = settings.hmm.mgm_prob;

	if ( logger->debug )   logger->Print( 10, Summary() );
	if ( logger->verbose ) logger->Print( 0, "# Loading parameters ... done" );
}
// ----------------------------------------------------
void Parameters::Load(void)
{
	ParameterParser  parser;

	if ( !mgm_filename.empty() )
	{
		parser.PutFileToBuffer( mgm_filename );

		// "__B_GC" and "__A_GC" are expected

		LoadSection( parser, parser.buffer, "__" );
	}

	if ( !native_filename.empty() )
	{
		parser.PutFileToBuffer( native_filename );

		LoadSection( parser, parser.buffer, "__NATIVE" );
	}

	if ( !tis_filename.empty() )
	{
		parser.PutFileToBuffer( tis_filename );

		LoadSection( parser, parser.buffer, "__START" );
	}

	if ( !model_set.size() )
		Exit( "error, model set is empty" );
}
// ----------------------------------------------------
void Parameters::LoadSection( ParameterParser & parser, std::string & buffer,  const std::string label )
{
	// for logger
	unsigned int set_size = model_set.size();

	// Parse 'buffer' into sections with "__" separator

	parameter_map  section;

	parser.LoadFromString( section, buffer, "__" );

	parameter_map::iterator itr     = section.begin();
	parameter_map::iterator itr_end = section.end();

	for( ; itr != itr_end; ++itr )
	{
		string section_name = itr->first;

		// skip sections whithout partial match to the label

		if ( section_name.compare( 0, label.size(), label ) ) 
			continue;

		if (logger->debug) logger->Print( 10, "Loading section: ", section_name  );

		// put empty model into map
	
		if ( model_set.find( section_name ) == model_set.end() )
		{
			model_set.insert( pair<string, Model>( section_name, Model(logger) ) );
		}
		else
			Exit( "error, section duplication found in parameter file:", section_name );

		// Parse model into variables using '$' separator

		parameter_map  par;

		parser.LoadFromSubString( par, buffer, "$", itr->second.first, itr->second.second );

		map< string, Model >::iterator current = model_set.find( section_name );

		// Fill in empty model

		// =====================================================================================================
		{
			if (parser.IsKey(par, "$NAME"))            current->second.model_name         = parser.asString ( par, "$NAME", buffer );
			if (parser.IsKey(par, "$GCODE"))           current->second.gcode              = parser.asPInt   ( par, "$GCODE", buffer );
			if (parser.IsKey(par, "$GENE_MIN_LENGTH")) current->second.gene_min_length    = parser.asPInt   ( par, "$GENE_MIN_LENGTH", buffer );
			if (parser.IsKey(par, "$COD_ORDER"))       current->second.order_cod          = parser.asPInt   ( par, "$COD_ORDER", buffer );
			if (parser.IsKey(par, "$NON_ORDER"))       current->second.order_non          = parser.asPInt   ( par, "$NON_ORDER", buffer );
			if (parser.IsKey(par, "$COD_P_N")) current->second.probability_N_in_coding    = parser.asPDouble( par, "$COD_P_N", buffer );
			if (parser.IsKey(par, "$NON_P_N")) current->second.probability_N_in_noncoding = parser.asPDouble( par, "$NON_P_N", buffer );
			if (parser.IsKey(par, "$NON_DURATION_DECAY")) current->second.noncoding_duration_decay = parser.asPDouble( par, "$NON_DURATION_DECAY", buffer );
			if (parser.IsKey(par, "$COD_DURATION_DECAY")) current->second.coding_duration_decay    = parser.asPDouble( par, "$COD_DURATION_DECAY", buffer );

			if (parser.IsKey(par, "$BUILD"))   current->second.build = parser.asString(par, "$BUILD", buffer);

			if ( parser.IsKey( par, "$ATG" ) )  current->second.pATG = parser.asPDouble( par, "$ATG", buffer );
			if ( parser.IsKey( par, "$GTG" ) )  current->second.pGTG = parser.asPDouble( par, "$GTG", buffer );
			if ( parser.IsKey( par, "$TTG" ) )  current->second.pTTG = parser.asPDouble( par, "$TTG", buffer );
			if ( parser.IsKey( par, "$TAA" ) )  current->second.pTAA = parser.asPDouble( par, "$TAA", buffer );
			if ( parser.IsKey( par, "$TAG" ) )  current->second.pTAG = parser.asPDouble( par, "$TAG", buffer );
			if ( parser.IsKey( par, "$TGA" ) )  current->second.pTGA = parser.asPDouble( par, "$TGA", buffer );
            
            // MGM PARAMETERS
            if ( parser.IsKey( par, "$ATG_A" ) )  current->second.pATG_A = parser.asPDouble( par, "$ATG_A", buffer );
            if ( parser.IsKey( par, "$GTG_A" ) )  current->second.pGTG_A = parser.asPDouble( par, "$GTG_A", buffer );
            if ( parser.IsKey( par, "$TTG_A" ) )  current->second.pTTG_A = parser.asPDouble( par, "$TTG_A", buffer );
            
            if ( parser.IsKey( par, "$ATG_B" ) )  current->second.pATG_B = parser.asPDouble( par, "$ATG_B", buffer );
            if ( parser.IsKey( par, "$GTG_B" ) )  current->second.pGTG_B = parser.asPDouble( par, "$GTG_B", buffer );
            if ( parser.IsKey( par, "$TTG_B" ) )  current->second.pTTG_B = parser.asPDouble( par, "$TTG_B", buffer );
            
            if ( parser.IsKey( par, "$ATG_C" ) )  current->second.pATG_C = parser.asPDouble( par, "$ATG_C", buffer );
            if ( parser.IsKey( par, "$GTG_C" ) )  current->second.pGTG_C = parser.asPDouble( par, "$GTG_C", buffer );
            if ( parser.IsKey( par, "$TTG_C" ) )  current->second.pTTG_C = parser.asPDouble( par, "$TTG_C", buffer );
            
            if ( parser.IsKey( par, "$ATG_D" ) )  current->second.pATG_D = parser.asPDouble( par, "$ATG_D", buffer );
            if ( parser.IsKey( par, "$GTG_D" ) )  current->second.pGTG_D = parser.asPDouble( par, "$GTG_D", buffer );
            if ( parser.IsKey( par, "$TTG_D" ) )  current->second.pTTG_D = parser.asPDouble( par, "$TTG_D", buffer );
            
            if ( parser.IsKey( par, "$ATG_X" ) )  current->second.pATG_X = parser.asPDouble( par, "$ATG_X", buffer );
            if ( parser.IsKey( par, "$GTG_X" ) )  current->second.pGTG_X = parser.asPDouble( par, "$GTG_X", buffer );
            if ( parser.IsKey( par, "$TTG_X" ) )  current->second.pTTG_X = parser.asPDouble( par, "$TTG_X", buffer );
            
            // MGM STOP PARAMETERS
            if ( parser.IsKey( par, "$TAA_A" ) )  current->second.pTAA_A = parser.asPDouble( par, "$TAA_A", buffer );
            if ( parser.IsKey( par, "$TAG_A" ) )  current->second.pTAG_A = parser.asPDouble( par, "$TAG_A", buffer );
            if ( parser.IsKey( par, "$TGA_A" ) )  current->second.pTGA_A = parser.asPDouble( par, "$TGA_A", buffer );
            
            if ( parser.IsKey( par, "$TAA_B" ) )  current->second.pTAA_B = parser.asPDouble( par, "$TAA_B", buffer );
            if ( parser.IsKey( par, "$TAG_B" ) )  current->second.pTAG_B = parser.asPDouble( par, "$TAG_B", buffer );
            if ( parser.IsKey( par, "$TGA_B" ) )  current->second.pTGA_B = parser.asPDouble( par, "$TGA_B", buffer );
            
            if ( parser.IsKey( par, "$TAA_C" ) )  current->second.pTAA_C = parser.asPDouble( par, "$TAA_C", buffer );
            if ( parser.IsKey( par, "$TAG_C" ) )  current->second.pTAG_C = parser.asPDouble( par, "$TAG_C", buffer );
            if ( parser.IsKey( par, "$TGA_C" ) )  current->second.pTGA_C = parser.asPDouble( par, "$TGA_C", buffer );
            
            if ( parser.IsKey( par, "$TAA_D" ) )  current->second.pTAA_D = parser.asPDouble( par, "$TAA_D", buffer );
            if ( parser.IsKey( par, "$TAG_D" ) )  current->second.pTAG_D = parser.asPDouble( par, "$TAG_D", buffer );
            if ( parser.IsKey( par, "$TGA_D" ) )  current->second.pTGA_D = parser.asPDouble( par, "$TGA_D", buffer );
            
            if ( parser.IsKey( par, "$TAA_X" ) )  current->second.pTAA_X = parser.asPDouble( par, "$TAA_X", buffer );
            if ( parser.IsKey( par, "$TAG_X" ) )  current->second.pTAG_X = parser.asPDouble( par, "$TAG_X", buffer );
            if ( parser.IsKey( par, "$TGA_X" ) )  current->second.pTGA_X = parser.asPDouble( par, "$TGA_X", buffer );
		
			if (parser.IsKey(par, "$NON_MAT") || parser.IsKey(par, "$COD_MAT"))
			{
				current->second.ReserveSpace();

				parser.asVectorOfDoublesWithLabel(par, "$NON_MAT", buffer, current->second.non);
				parser.as3VectorOfDoublesWithLabel(par, "$COD_MAT", buffer, current->second.cod1, current->second.cod2, current->second.cod3);
			}
		}
		// =====================================================================================================
		{
			LoadSite( parser, par, buffer, "$RBS",          &current->second.RBS,                  true );
            LoadSite( parser, par, buffer, "$PROMOTER",     &current->second.Promoter,             true );
            LoadSite( parser, par, buffer, "$SC",           &current->second.StartContent,         false );
            LoadSite( parser, par, buffer, "$SC_RBS",       &current->second.StartContentRBS,      false );
            LoadSite( parser, par, buffer, "$SC_PROMOTER",  &current->second.StartContentPromoter, false );
            
            
            // MGM PARAMETERS
            LoadMultiShiftSite( parser, par, buffer, "$RBS_A",  &current->second.RBS_A,  true );
            LoadMultiShiftSite( parser, par, buffer, "$RBS_B",  &current->second.RBS_B,  true );
            LoadMultiShiftSite( parser, par, buffer, "$RBS_C",  &current->second.RBS_C,  true );
            LoadMultiShiftSite( parser, par, buffer, "$RBS_D",  &current->second.RBS_D,  true );
            LoadMultiShiftSite( parser, par, buffer, "$RBS_X",  &current->second.RBS_X,  true );
            
            LoadMultiShiftSite( parser, par, buffer, "$PROMOTER_C",     &current->second.PROMOTER_C,             true );
            LoadMultiShiftSite( parser, par, buffer, "$PROMOTER_D",     &current->second.PROMOTER_D,             true );
            
            LoadSite( parser, par, buffer, "$SC_RBS_A",       &current->second.SC_RBS_A,      false );
            LoadSite( parser, par, buffer, "$SC_RBS_B",       &current->second.SC_RBS_B,      false );
            LoadSite( parser, par, buffer, "$SC_RBS_C",       &current->second.SC_RBS_C,      false );
            LoadSite( parser, par, buffer, "$SC_RBS_D",       &current->second.SC_RBS_D,      false );
            LoadSite( parser, par, buffer, "$SC_RBS_X",       &current->second.SC_RBS_X,      false );
            LoadSite( parser, par, buffer, "$SC_PROMOTER_C",  &current->second.SC_PROMOTER_C, false );
            LoadSite( parser, par, buffer, "$SC_PROMOTER_D",  &current->second.SC_PROMOTER_D, false );
		}
		// =====================================================================================================
		{
			if ( parser.IsKey( par, "$TO_ATYPICAL_FIRST_BACTERIA" ) )  to_atypical_first  = parser.asPDouble( par, "$TO_ATYPICAL_FIRST_BACTERIA", buffer );
			if ( parser.IsKey( par, "$TO_ATYPICAL_SECOND_ARCHAEA" ) )  to_atypical_second = parser.asPDouble( par, "$TO_ATYPICAL_SECOND_ARCHAEA", buffer );
			if ( parser.IsKey( par, "$TO_NATIVE" ) )                   to_native          = parser.asPDouble( par, "$TO_NATIVE", buffer );
			if ( parser.IsKey( par, "$TO_MGM" ) )                      to_mgm             = parser.asPDouble( par, "$TO_MGM", buffer );

			if ( parser.IsKey( par, "$GENE_GC_DIST" ) && parser.asBool(par, "$GENE_GC_DIST", buffer ) )
			{
				to_native_by_gc.assign( 101, 0 );

				parser.asVectorOfDoublesWithPos( par, "$GENE_GC_DIST_MAT", buffer, to_native_by_gc );
			}
		}
	}

	if (logger->debug) logger->Print( 10, "Sections loaded:", label,  model_set.size() - set_size  );
}
// ----------------------------------------------------
void Parameters::LoadSite( ParameterParser & parser, parameter_map & par, std::string & buffer, std::string const label, Site * ptr, bool with_dur )
{
	if ( parser.IsKey( par, label ) && parser.asBool( par, label, buffer ) )
	{
		ptr->order    = parser.asPInt( par, label + "_ORDER",   buffer );
		ptr->width    = parser.asPInt( par, label + "_WIDTH",   buffer );
		ptr->margin   = parser.asInt ( par, label + "_MARGIN",  buffer );

		if ( with_dur )
			ptr->max_dur  = parser.asPInt( par, label + "_MAX_DUR", buffer ) + 1;

		ptr->ReserveSpace();

		parser.asMatrixOfDoublesWithLabel( par, label + "_MAT", buffer, ptr->matrix );

		if ( with_dur)
			parser.asVectorOfDoublesWithPos( par, label + "_POS_DISTR", buffer, ptr->duration );

		ptr->is_valid = true;
	}
}

void Parameters::LoadMultiShiftSite( ParameterParser & parser, parameter_map & par, std::string & buffer, std::string const label, MultiShiftSite * ptr_ms, bool with_dur )
{
    int max_allowed_shift = 5;
    
    
    
    for (int shift = 0; shift < max_allowed_shift; shift++) {
        
        stringstream iss;
        iss << label;
        iss << "_" << shift;
        
        // get new label with shift by adding _SHIFT to end
        std::string label_ws = iss.str();
        
        if ( parser.IsKey( par, label_ws ) && parser.asBool( par, label_ws, buffer ) )
        {
            Site *ptr = new Site();        // this memory is freed in destructor of MultiShiftSite
            
            ptr->order    = parser.asPInt( par, label_ws + "_ORDER",   buffer );
            ptr->width    = parser.asPInt( par, label_ws + "_WIDTH",   buffer );
            ptr->margin   = parser.asInt ( par, label_ws + "_MARGIN",  buffer );

            if ( with_dur )
                ptr->max_dur  = parser.asPInt( par, label_ws + "_MAX_DUR", buffer ) + 1;

            ptr->ReserveSpace();

            parser.asMatrixOfDoublesWithLabel( par, label_ws + "_MAT", buffer, ptr->matrix );

            if ( with_dur)
                parser.asVectorOfDoublesWithPos( par, label_ws + "_POS_DISTR", buffer, ptr->duration );

            ptr->is_valid = true;
            
            float shift_prior = parser.asPDouble(par, label_ws + "_SHIFT",   buffer );
            
            ptr_ms->is_valid = true;
            ptr_ms->add_site_with_shift(ptr, shift_prior);
        }
    }
}
// ----------------------------------------------------
string Parameters::Summary(void)
{
	stringstream s;

	s << "# START Parameters:" << endl;
	s << "parameter sets loaded into 'model_set' : " << model_set.size() << endl;

	s << "to_native: " << to_native << endl;
	s << "to_mgm: " << to_mgm << endl;
	s << "to_atypical_first: " << to_atypical_first << endl;
	s << "to_atypical_second: " << to_atypical_second << endl;

	s << "# STAR KL of model" << endl;

	map< string, Model >::iterator  itr     = model_set.begin();
	map< string, Model >::iterator  itr_end = model_set.end();

	for( ; itr != itr_end; ++itr )
	{
		s << itr->first << " " << itr->second.KL() << endl;
	}

	s << "# END KL of model" << endl;
	s << "# END Parameters " << endl;

	return s.str();
}
// ----------------------------------------------------
void Parameters::Initialize( std::map< std::string, Model > & target )
{
	map< string, Model >::iterator  itr     = target.begin();
	map< string, Model >::iterator  itr_end = target.end();

	for( ; itr != itr_end; ++itr )
	{
		itr->second.Initialize();
	}
}
// ----------------------------------------------------
