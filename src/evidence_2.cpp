// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: evidence_2.cpp
// Code was tested with boost 1.48 <boost/algorithm/string.hpp>
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <vector>
#include <string>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <set>

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::set;
using std::fill;
using std::pair;

#include <boost/algorithm/string.hpp>

#include "evidence_2.h"
#include "exit.h"

// ----------------------------------------------------
Evidence::Evidence( Settings const & settings, Logger * const logger ) : logger(logger)
{
	filename = settings.in.evidence;

	if ( filename.empty() )  return;

	SetEvidenceLabels();
	ParseFile( filename, data );
	
	if( logger->verbose ) logger->Print( "# Loading evidence ... done" );
	if( logger->debug )   logger->Print( Status(logger->debug) );
}
// ----------------------------------------------------
void Evidence::SetEvidenceLabels(void)
{
	str_to_evidence.clear();

	// all letters to upper case

	string str;

	str.assign(TAGS::COD_EVIDENCE);  boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::COD;
	str.assign(TAGS::NON_EVIDENCE);  boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::NON;
	str.assign(TAGS::MASK_EVIDENCE); boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::MASK;
	str.assign(TAGS::RDNA_EVIDENCE); boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::RDNA;
	str.assign(TAGS::TRNA_EVIDENCE); boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::TRNA;
	str.assign(TAGS::NRNA_EVIDENCE); boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::NRNA;
	str.assign(TAGS::FSH_EVIDENCE);  boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::FSH;
	str.assign(TAGS::CDS_EVIDENCE);  boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::CDS;

	str.assign(TAGS::C_EVIDENCE); boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::COD;
	str.assign(TAGS::N_EVIDENCE); boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::NON;
	str.assign(TAGS::M_EVIDENCE); boost::to_upper(str); str_to_evidence[str] = EVIDENCE_TYPES::MASK;

}
// ----------------------------------------------------
void Evidence::ParseFile( std::string const & name, EvidenceMap & target )
{
	target.clear();

	ifstream gff_file( name.c_str() );

	if ( !gff_file.is_open() )
		Exit( "error on open file:", name );

	string line;

	while( getline( gff_file, line ) )
	{
		string::size_type pos = line.find_first_not_of( " \t" );

		if ( pos == string::npos ) continue;
		if ( line[pos] == '#' )    continue;

		INTERVAL_EVI interval;
		string seq_name;

		if ( ParseGFFLine( line, seq_name, interval ) )
		{
			if ( interval.key == EVIDENCE_TYPES::COD || interval.key == EVIDENCE_TYPES::CDS )
			{
				if ( IsBadCodInterval(interval) )
				{
					if (logger->verbose ) logger->Print( 0, "warning, line from GFF file was ignored:", line );
					
					continue;
				}
			}

			target[seq_name].push_back( interval );
		}
		else
		{
			if (logger->debug) logger->Print( 10, "warning, line from GFF file was ignored:", line );
		}
	}
}
// ----------------------------------------------------
bool Evidence::IsBadCodInterval( INTERVAL_EVI const & rec )
{
	if ( (rec.right - rec.left + 1)%3 != 0 )
	{
		if (logger->verbose ) logger->Print( 0, "warning, length is not 3x" );

		return true;
	}

	if ( rec.strand != STRAND_TYPES::DIRECT && rec.strand != STRAND_TYPES::REVERSE )
	{
		if (logger->verbose ) logger->Print( 0, "warning, direct or reverse strand is expected" );

		return true;
	}

	return false;
}
// ----------------------------------------------------
std::string Evidence::Status(unsigned int const debug_level)
{
	stringstream s;

	s << "# START Evidence:" << endl;
	s << "data.size() " << data.size() << endl;

	int gff_lines_accepted = 0;

	for( EvidenceMapItr itr = data.begin(); itr != data.end(); ++itr )
			gff_lines_accepted = itr->second.size();

	s << "gff_lines_accepted " << gff_lines_accepted << endl;

	if (debug_level > 2)
	{
		s << "seqname number_of_intervals" << endl;

		for( EvidenceMapItr itr = data.begin(); itr != data.end(); ++itr )
		{
			s << itr->first << " " << itr->second.size()  << endl;
		}

		// print all accepted GFF lines

		std::map< EVIDENCE_TYPE, std::string > evidence_to_str;

		evidence_to_str[EVIDENCE_TYPES::COD]  = TAGS::COD_EVIDENCE;
		evidence_to_str[EVIDENCE_TYPES::NON]  = TAGS::NON_EVIDENCE;
		evidence_to_str[EVIDENCE_TYPES::MASK] = TAGS::MASK_EVIDENCE;
		evidence_to_str[EVIDENCE_TYPES::RDNA] = TAGS::RDNA_EVIDENCE;
		evidence_to_str[EVIDENCE_TYPES::TRNA] = TAGS::TRNA_EVIDENCE;
		evidence_to_str[EVIDENCE_TYPES::NRNA] = TAGS::NRNA_EVIDENCE;
		evidence_to_str[EVIDENCE_TYPES::FSH]  = TAGS::FSH_EVIDENCE;
		evidence_to_str[EVIDENCE_TYPES::CDS]  = TAGS::CDS_EVIDENCE;

		for( map< string, vector<INTERVAL_EVI> >::iterator itr = data.begin(); itr != data.end(); ++itr )
		{
			for( vector<INTERVAL_EVI>::iterator i = itr->second.begin();  i != itr->second.end(); ++i )
			{
				//       1                  2                      3                           4                  5                   6                   7                   8               9
				s << itr->first << "\t" << "." << "\t" << evidence_to_str[i->key] << "\t" << i->left << "\t" << i->right << "\t" << i->score << "\t" << i->strand << "\t" << "."  << "\t" << i->attr << endl;
			}
		}
	}

	s << "# END Evidence " << endl;

	return s.str();
}
// ----------------------------------------------------
bool Evidence::ParseGFFLine( std::string const & source, std::string & seq_name, INTERVAL_EVI & target )
{
	// check for 9 TAB separated fields

	int count_tabs = 0;

	for ( unsigned int i = 0; i < source.size(); ++i )
	{
		if ( source[i] == '\t' )
			++count_tabs;
	}

	if ( count_tabs < 8 )
		Exit( "error, unexpected number of TABs in GFF line:", count_tabs, source );

	// 1 column: seqname

	size_t pos_L = 0;
	size_t pos_R = source.find_first_of( '\t', pos_L );
	{
		seq_name.assign( source.c_str() + pos_L, pos_R - pos_L );

		if ( seq_name.empty() )
			Exit( "error, empty 'seqname' field in GFF line found:", source );
	}

	// 2 column: source
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( '\t', pos_L );
	
	// 3 column: feature
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( '\t', pos_L );
	{
		string feature( source.c_str() + pos_L, pos_R - pos_L );

		// check for allowed action keys
	
		boost::to_upper( feature );

		map< std::string, EVIDENCE_TYPE >::iterator itr = str_to_evidence.find( feature );

		if ( itr != str_to_evidence.end() )
			target.key = itr->second;
		else
		{
			if (logger->debug ) logger->Print( 10, "warning, feature was ignored", feature );
			return false;
		}
	}
	
	// 4 column: start
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( '\t', pos_L );
	{
		target.left = atoi( source.c_str() + pos_L );
	}
	
	// 5 column: end
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( '\t', pos_L );
	{
		target.right = atoi( source.c_str() + pos_L );

		// check coordinates
		if ( target.left < 1 || target.left > target.right  )
			Exit( "error in GFF 'start' 'end' position values:", target.left, target.right );
	}

	// 6 column: score
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( '\t', pos_L );
	{
		if ( source.at(pos_L) == '.' )
			target.score = VALUE_ZERO;
		else
			target.score = atof( source.c_str() + pos_L );
	}
	
	// 7 column: strand
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( '\t', pos_L );
	{
		if ( pos_R - pos_L != 1 )
			Exit( "error, unexpectd label found in GFF strand field:", source.substr( pos_L, pos_R - pos_L ) );

		if ( source.at(pos_L) == '.' ) 
			target.strand = STRAND_TYPES::BOTH;
		else if ( source.at(pos_L) == '+' )
			target.strand = STRAND_TYPES::DIRECT;
		else if ( source.at(pos_L) == '-' )
			target.strand = STRAND_TYPES::REVERSE;
		else
			Exit( "error, unexpectd label found in GFF strand field:", source.substr( pos_L, pos_R - pos_L ) );
	}

	// 8 column: frame
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( '\t', pos_L );
	
	// 9 column: attribute
	
	pos_L = source.find_first_not_of( '\t', pos_R );
	pos_R = source.find_first_of( " \t", pos_L );

	target.attr = source.substr( pos_L, pos_R - pos_L );

	return true;
}
// ----------------------------------------------------
void Evidence::SoftMaskToIntervals( FastaRecord const & rec )
{
	string::size_type L;
	string::size_type R;

	L = rec.data.find_first_of( "acgturyswkmbdhvn", 0 );

	while ( L != string::npos )
	{
		R = rec.data.find_first_not_of("acgturyswkmbdhvn", L );

		INTERVAL_EVI interval;

		// GFF positions are "1" based

		interval.left = static_cast<int>( L ) + 1;

		if ( R != string::npos)
			interval.right = static_cast<int>( R );
		else
			interval.right = static_cast<int>( rec.data.size() );

		interval.strand = STRAND_TYPES::BOTH;
		interval.score  = VALUE_ZERO;
		interval.key    = EVIDENCE_TYPES::MASK;

		data[ rec.name ].push_back( interval );

		L = rec.data.find_first_of( "acgturyswkmbdhvn", R );
	}
}
// ----------------------------------------------------
void Evidence::SoftMaskToEvidence( FastaVector const & fasta )
{
	for ( vector< FastaRecord >::const_iterator itr = fasta.begin(); itr != fasta.end(); ++itr )
	{
		SoftMaskToIntervals( *itr );
	}
}
// ----------------------------------------------------
void Evidence::SyncSeqNameWithSizeOne(FastaVector const & fasta, EvidenceMap & evi)
{
	typedef map< string, vector<INTERVAL_EVI> >::iterator evi_i;

	// rename in evi when name differs
	// rename in map by: insert - swap - delete

	if ( fasta.at(0).name.compare( evi.begin()->first ) ) 
	{
		string new_name = fasta.at(0).name;
		string old_name = evi.begin()->first;

		vector<INTERVAL_EVI> tmp;
		data[ new_name] = tmp;

		evi_i  itr_to_new = evi.find(new_name);
		evi_i  itr_to_old = evi.find(old_name);

		std::swap( itr_to_old->second, itr_to_new->second );

		evi.erase( itr_to_old );

		if ( logger->verbose ) logger->Print( 0, "warning, 'seqname' in GFF was renamed from-to:", old_name, new_name );
	}
}
// ----------------------------------------------------
void Evidence::VerifySeqNames( FastaVector const & fasta, EvidenceMap const & evi)
{
	typedef vector< FastaRecord >::const_iterator  fasta_i;
	typedef map< string, vector<INTERVAL_EVI> >::const_iterator  evi_i;

	set<string> names;

	for ( fasta_i fasta_itr = fasta.begin(); fasta_itr != fasta.end(); ++fasta_itr )
	{
		names.insert(fasta_itr->name);
	}

	for ( evi_i  evi_itr = evi.begin(); evi_itr != evi.end(); ++evi_itr )
	{
		if ( names.find( evi_itr->first ) == names.end() )
		{
			logger->Print( 0, "warning, 'seqname' from GFF not found in input sequence:", evi_itr->first );
		}
	}
}
// ----------------------------------------------------
void Evidence::CheckIntervalsRightBound( FastaVector const & fasta)
{
	for ( vector< FastaRecord >::const_iterator itr = fasta.begin(); itr != fasta.end(); ++itr )
	{
		map< string, vector<INTERVAL_EVI> >::iterator  gff = data.find( itr->name );

		if ( gff == data.end() ) 
			continue;

		int size = itr->data.size();

		for( vector<INTERVAL_EVI>::iterator interval = gff->second.begin(); interval != gff->second.end();  )
		{
			if ( interval->right > size )
			{
				gff->second.erase( interval );

				logger->Print( 0, "warning, interval is out of sequence:", interval->right, itr->name );
			}
			else
				++interval;
		}
	}
}
// ----------------------------------------------------
void Evidence::SyncWithSequence( FastaVector const & fasta )
{
	if ( data.size() == 0 )
		return;

	if ( data.size() == 1 && fasta.size() == 1 )
		SyncSeqNameWithSizeOne( fasta, data );
	else
		VerifySeqNames( fasta, data );

	CheckIntervalsRightBound( fasta );
}
// ----------------------------------------------------
void Evidence::HardMaskSequence( FastaVector & fasta )
{
	typedef vector< FastaRecord >::iterator  fasta_i;
	typedef map< string, vector<INTERVAL_EVI> >::const_iterator  gff_i;
	typedef vector< INTERVAL_EVI >::const_iterator  interval_i;

	for( fasta_i  fasta_itr = fasta.begin(); fasta_itr != fasta.end(); ++fasta_itr )
	{
		gff_i  gff_itr = data.find( fasta_itr->name );
		
		if ( gff_itr != data.end() )
		{
			for( interval_i  interval_itr = gff_itr->second.begin(); interval_itr != gff_itr->second.end(); ++interval_itr )
			{
				if ( interval_itr->key == EVIDENCE_TYPES::MASK )
				{
					fill( fasta_itr->data.begin() + interval_itr->left - 1, fasta_itr->data.begin() + interval_itr->right, NT::N );
				}
			}
		}
	}
}
// ----------------------------------------------------
