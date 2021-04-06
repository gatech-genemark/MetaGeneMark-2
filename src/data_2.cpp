// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: data_2.cpp
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::stringstream;
using std::map;

#include "data_2.h"
#include "common_2.h"
#include "index_2.h"
#include "exit.h"

// ----------------------------------------------------
Data::Data( Settings & settings, Pset & pset, Logger * const logger ) : logger(logger)
{
	incomplete_at_gaps = settings.in.incomplete_at_gaps;
	softmask = settings.in.softmasked;
	strand = settings.in.strand_t;
	gcode = pset.genetic_code;
	min_gene_length = pset.min_gene_length;

	SetNucDataAlphabet( softmask, nuc_data_alphabet );
	SetCodonToFunctionTable( gcode, strand, stop_and_start_codon_alphabet );
}
// ----------------------------------------------------
void Data::Set(std::string const & source)
{
	LoadNT(source);
	CountCumulativeGC(nt);
	SetFlags(nt);

	if ( incomplete_at_gaps > 0 )
		SetFlagsAtGaps(incomplete_at_gaps);

	if (logger->debug) logger->Print( 3, PrintNT() );
	if (logger->debug) logger->Print( 2, PrintFlags() );
}
// ----------------------------------------------------
void Data::LoadNT(std::string const & source)
{
	nt.assign( source.size(), 0 );

	for ( unsigned int i = 0; i < source.size(); ++i )
	{
		nt[i] = nuc_data_alphabet[ source[i] ];
	}
}
// ----------------------------------------------------
void Data::SetNucDataAlphabet( bool const mask, std::vector<char> & target )
{
	target.assign( 128, NT::N );

	target[ 'A' ] = NT::A;
	target[ 'C' ] = NT::C;
	target[ 'G' ] = NT::G;
	target[ 'T' ] = NT::T;

	if ( !mask )
	{
		target[ 'a' ] = NT::A;
		target[ 'c' ] = NT::C;
		target[ 'g' ] = NT::G;
		target[ 't' ] = NT::T;
	}
}
// ----------------------------------------------------
void Data::CountCumulativeGC(std::vector<unsigned char> const & source)
{
	gc.assign( source.size(), 0 );

	std::vector<char> nt_to_gc( 64, 1 );

	nt_to_gc[NT::G] = 2;
	nt_to_gc[NT::C] = 2;

	nt_to_gc[NT::A] = 0;
	nt_to_gc[NT::T] = 0;

	int sum = 0;

	for( unsigned int i = 0; i < source.size(); ++i )
	{
		sum += nt_to_gc[ source[i] ];
		gc[i] = sum;
	}
}
// ----------------------------------------------------
// L & R - C-style inclusive indeces of vector<int>
//
// X1, X1+X2, X1+X2+X3, X1+X2+X3+X4, X1+X2+X3+X4+X5, ...
//      L-1      L                           R
//  gc( L, R ) = [R] - [L-1] = (X1+X2+X3+X4+X5)-(X1+X2) = X3+X4+X5
// 
// Letters 'N' makes GC calculation slow
// For speeding up the process (on expence of accuracy of calculation): GC is 2, AT is 0 and N is 1
// ----------------------------------------------------
unsigned int Data::GCasINT( unsigned int const L, unsigned int const R )
{
	if ( L > R )
		Exit( "error in 'Data::GCasINT' interval range L..R:", L, R );

	if (R > gc.size())
		Exit( "error in 'Data::GCasINT' interval boundary R:", L, gc.size() );

	double x;

	if ( L >= 1 )
	{
		x = ( gc[R] - gc[L-1] )/( 2.0*(R-L+1) );
	}
	else  // L == 0
	{
		if ( L != R )
			x = gc[R]/( 2.0*(R+1) );
		else 
			x = gc[L]/2.0;
	}

	return (unsigned int)( 100.0 * x + 0.5 );
}
// ----------------------------------------------------
void Data::SetFlags(std::vector<unsigned char> const & source)
{
	flag.assign( source.size(), 0 );

	struct Frame
	{
		unsigned int frame;
		Frame* next;
	};
	
	Frame fr1;
	Frame fr2;
	Frame fr3;

	fr1.frame = frame1;
	fr2.frame = frame2; 
	fr3.frame = frame3;

	fr1.next = &fr2;
	fr2.next = &fr3;
	fr3.next = &fr1;
	
	Frame* fr;
	fr = &fr1;

	// triplet index to mark start and stop codons

	Index codon(2);

	unsigned int idx;
	unsigned int const npos = -1;

	for( unsigned int i = 0; i < source.size(); ++i )
	{
		idx = codon.Get(source[i]);

		if ( idx != npos )
		{
			flag[i] = stop_and_start_codon_alphabet[ idx ];
		}

		flag[i] |= fr->frame;

		fr = fr->next;
	}
}
// ----------------------------------------------------
void Data::MergeGaps(std::vector<INTERVAL> & d, int max )
{
	// do nothing on empty input
	if (d.size() == 0)
		return;

	// check from ini to first gap

	if (d.begin()->left < max)
	{
		// clean labels here

		for (int k = 0; k < d.begin()->left; ++k)
		{
			flag[k] = 0;
		}
	}

	// check from last gap to term

	if (nt.size() - d.back().right < max)
	{
		// clean labels here

		for (int k = d.back().right; k < nt.size(); ++k)
		{
			flag[k] = 0;
		}
	}

	// more than two gaps
	if (d.size() > 1)
	{
		int i = 0;
		int j = 1;

		while (j < d.size())
		{
			if (d[j].left - d[i].right < max)
			{
				// clean flags between the gaps

				for (int k = d[i].right; k < d[j].left; ++k)
				{
					flag[k] = 0;
				}

				// merge gaps:
				//    * extend the right gap bounday
				//    * zero the left gap

				d[j].left = d[i].left;

				d[i].left = 0;
				d[i].right = 0;
			}

			++i;
			++j;
		}
	}
}
// ----------------------------------------------------
void Data::SetFlagsAtGaps(unsigned int len)
{
	// "gap" is a continuous stretch of symbols "NT::N" longer than "threhold"
	// minimum "gap" length is currently hard coded

	const int min_gap = 6 + 6;

	if ( len < min_gap)
	{
		len = min_gap + 1;
		logger->Print( 0, "warning, minimum gap length was reset from - to ", len, min_gap );
	}

	unsigned int size = nt.size();
	unsigned int L = 0;
	unsigned int R = 0;
	vector<INTERVAL> gaps;

	// return, if "len" is more than sequence length
	if ( len > size )
		return;

	while ( L < size )
	{
		// find gap start
		while (L < size && nt[L] != NT::N)
			++L;

		// no more "N"
		if (L == size)
			break;

		// find gap end
		R = L;
		while ( R < size && nt[R] == NT::N )
			++R;

		// record
		if ( R - L >= len )
		{
			INTERVAL interval;

			interval.left = L;
			interval.right = R;

			gaps.push_back(interval);
		}

		// look for next from last position
		L = R;
	}	

	// gaps not found

	if (gaps.size() == 0)
		return;

	// merge closely located gaps

	MergeGaps(gaps, min_gene_length );

	// set gap flags
	for( vector<INTERVAL>::iterator itr = gaps.begin(); itr < gaps.end(); ++itr )
	{
		L = itr->left;
		R = itr->right;

		// one from the two merged gaps is reset to zero, so skip such gaps zero gaps
		if (R == 0)
			continue;
		
		if (L > 3)
		{
			if (!(flag[L - 3] & dirEnd))
			{
				flag[L + 0] |= isGap;
				flag[L + 0] |= dirEnd;
			}

			if (!(flag[L - 2] & dirEnd))
			{
				flag[L + 1] |= isGap;
				flag[L + 1] |= dirEnd;
			}

			if (!(flag[L - 1] & dirEnd))
			{
				flag[L + 2] |= isGap;
				flag[L + 2] |= dirEnd;
			}

			flag[L + 3] |= isGap;
			flag[L + 3] |= revStart;
			flag[L + 4] |= isGap;
			flag[L + 4] |= revStart;
			flag[L + 5] |= isGap;
			flag[L + 5] |= revStart;
		}
		
		if (R < nt.size() - 4)
		{
			if (!(flag[R + 2] & revEnd))
			{
				flag[R - 1] |= isGap;
				flag[R - 1] |= revEnd;
			}

			if (!(flag[R + 3] & revEnd))
			{
				flag[R + 0] |= isGap;
				flag[R + 0] |= revEnd;
			}

			if (!(flag[R + 4] & revEnd))
			{
				flag[R + 1] |= isGap;
				flag[R + 1] |= revEnd;
			}

			flag[R - 2] |= isGap;
			flag[R - 2] |= dirStart;
			flag[R - 3] |= isGap;
			flag[R - 3] |= dirStart;
			flag[R - 4] |= isGap;
			flag[R - 4] |= dirStart;
		}
	}
}
// ----------------------------------------------------
unsigned int Data::CountFlags( std::vector<unsigned int> const & source )
{
	unsigned int count = 0;

	for( unsigned int i = 0; i < source.size(); ++i )
	{
		if ( source[i] != 0 )
			++count;
	}

	return count;
}
// ----------------------------------------------------
std::string Data::PrintNT(void)
{
	stringstream s;

	s << "# nucleotide sequence from current record " << endl;

	for ( unsigned int i = 0; i < nt.size(); ++i )
	{	
		switch( nt[i] )
		{
			case NT::A:
				s << "A";
				break;
			case NT::C:
				s << "C";
				break;
			case NT::G:
				s << "G";
				break;
			case NT::T:
				s << "T";
				break;
			case NT::N:
				s << "N";
				break;
		}
	}
	
	return s.str();
}
// ----------------------------------------------------
std::string Data::PrintFlags(void)
{
	stringstream s;
	int count_flags = 0; 

	s << "# flags from current record " << endl;

	for ( unsigned int i = 0; i < flag.size(); ++i )
	{	
		if ( ! (flag[i] & emptyMark) )
			continue;
		
		++count_flags;

		s << i << " ";
			
		if ( flag[i] & dirEnd )   s << "d_end ";
		if ( flag[i] & revEnd )   s << "r_end ";
		if ( flag[i] & dirStart ) s << "d_sta ";
		if ( flag[i] & revStart ) s << "r_sta ";
		if ( flag[i] & isGap )    s << "gap ";
		if ( flag[i] & frame1 )   s << "fr1 ";
		if ( flag[i] & frame2 )   s << "fr2 ";
		if ( flag[i] & frame3 )   s << "fr3 ";
		if ( flag[i] & isATG )    s << "atg ";
		if ( flag[i] & isGTG )    s << "gtg ";
		if ( flag[i] & isTTG )    s << "ttg ";
		if ( flag[i] & isTAA )    s << "taa ";
		if ( flag[i] & isTAG )    s << "tag ";
		if ( flag[i] & isTGA )    s << "tga ";
		if ( flag[i] & iniCod )   s << "ini_c ";
		if ( flag[i] & iniNon )   s << "ini_n ";
		if ( flag[i] & terCod )   s << "ter_c ";
		if ( flag[i] & terNon )   s << "ter_n ";

		s << endl;
	}

	s << "total non zero flags: " << count_flags << endl;

	return s.str();
}
// ----------------------------------------------------
void Data::SetCodonToFunctionTable( unsigned int const gcode, STRAND_TYPE const st, std::vector<unsigned int> & target )
{
	target.assign( 64, 0 );

	if ( st == STRAND_TYPES::BOTH || st == STRAND_TYPES::DIRECT )
	{
		switch(gcode)
		{
			case 1:
				target[ NT::A<<4|NT::T<<2|NT::G ] = dirStart + isATG;
				target[ NT::T<<4|NT::A<<2|NT::A ] = dirEnd + isTAA;
				target[ NT::T<<4|NT::G<<2|NT::A ] = dirEnd + isTGA;
				target[ NT::T<<4|NT::A<<2|NT::G ] = dirEnd + isTAG;
				break;
			case 4:
			case 25:
				target[ NT::A<<4|NT::T<<2|NT::G ] = dirStart + isATG;
				target[ NT::G<<4|NT::T<<2|NT::G ] = dirStart + isGTG;
				target[ NT::T<<4|NT::T<<2|NT::G ] = dirStart + isTTG;
				target[ NT::T<<4|NT::A<<2|NT::A ] = dirEnd + isTAA;
				target[ NT::T<<4|NT::A<<2|NT::G ] = dirEnd + isTAG;
				break;
			case 11:
				target[ NT::A<<4|NT::T<<2|NT::G ] = dirStart + isATG;
				target[ NT::G<<4|NT::T<<2|NT::G ] = dirStart + isGTG;
				target[ NT::T<<4|NT::T<<2|NT::G ] = dirStart + isTTG;
				target[ NT::T<<4|NT::A<<2|NT::A ] = dirEnd + isTAA;
				target[ NT::T<<4|NT::G<<2|NT::A ] = dirEnd + isTGA;
				target[ NT::T<<4|NT::A<<2|NT::G ] = dirEnd + isTAG;
				break;
			case 15:
				target[ NT::A<<4|NT::T<<2|NT::G ] = dirStart + isATG;
				target[ NT::G<<4|NT::T<<2|NT::G ] = dirStart + isGTG;
				target[ NT::T<<4|NT::T<<2|NT::G ] = dirStart + isTTG;
				target[ NT::T<<4|NT::A<<2|NT::A ] = dirEnd + isTAA;
				target[ NT::T<<4|NT::G<<2|NT::A ] = dirEnd + isTGA;
				break;
			case 101:
				target[ NT::A<<4|NT::T<<2|NT::G ] = dirStart + isATG;
				target[ NT::G<<4|NT::T<<2|NT::G ] = dirStart + isGTG;
				target[ NT::T<<4|NT::T<<2|NT::G ] = dirStart + isTTG;
				target[ NT::T<<4|NT::G<<2|NT::A ] = dirEnd + isTGA;
				target[ NT::T<<4|NT::A<<2|NT::G ] = dirEnd + isTAG;
				break;
			default:
				Exit( "error, unsupported genetic code found:", gcode );
		}
	}
	
	if ( st == STRAND_TYPES::BOTH || st == STRAND_TYPES::REVERSE )
	{
		switch(gcode)
		{
			case 1:
				target[ NT::C<<4|NT::A<<2|NT::T ] = revStart + isATG;
				target[ NT::T<<4|NT::T<<2|NT::A ] = revEnd + isTAA;
				target[ NT::T<<4|NT::C<<2|NT::A ] = revEnd + isTGA;
				target[ NT::C<<4|NT::T<<2|NT::A ] = revEnd + isTAG;
				break;
			case 4:
			case 25:
				target[ NT::C<<4|NT::A<<2|NT::T ] = revStart + isATG;
				target[ NT::C<<4|NT::A<<2|NT::C ] = revStart + isGTG;
				target[ NT::C<<4|NT::A<<2|NT::A ] = revStart + isTTG;
				target[ NT::T<<4|NT::T<<2|NT::A ] = revEnd + isTAA;
				target[ NT::C<<4|NT::T<<2|NT::A ] = revEnd + isTAG;
				break;
			case 11:
				target[ NT::C<<4|NT::A<<2|NT::T ] = revStart + isATG;
				target[ NT::C<<4|NT::A<<2|NT::C ] = revStart + isGTG;
				target[ NT::C<<4|NT::A<<2|NT::A ] = revStart + isTTG;
				target[ NT::T<<4|NT::T<<2|NT::A ] = revEnd + isTAA;
				target[ NT::T<<4|NT::C<<2|NT::A ] = revEnd + isTGA;
				target[ NT::C<<4|NT::T<<2|NT::A ] = revEnd + isTAG;
				break;
			case 15:
				target[ NT::C<<4|NT::A<<2|NT::T ] = revStart + isATG;
				target[ NT::C<<4|NT::A<<2|NT::C ] = revStart + isGTG;
				target[ NT::C<<4|NT::A<<2|NT::A ] = revStart + isTTG;
				target[ NT::T<<4|NT::T<<2|NT::A ] = revEnd + isTAA;
				target[ NT::T<<4|NT::C<<2|NT::A ] = revEnd + isTGA;
				break;
			case 101:
				target[ NT::C<<4|NT::A<<2|NT::T ] = revStart + isATG;
				target[ NT::C<<4|NT::A<<2|NT::C ] = revStart + isGTG;
				target[ NT::C<<4|NT::A<<2|NT::A ] = revStart + isTTG;
				target[ NT::T<<4|NT::C<<2|NT::A ] = revEnd + isTGA;
				target[ NT::C<<4|NT::T<<2|NT::A ] = revEnd + isTAG;
				break;
			default:
				Exit( "error, unsupported genetic code found:", gcode );
		}
	}
}
// ----------------------------------------------------
void Data::ApplyEvidence( std::map< std::string, std::vector<INTERVAL_EVI> > const & evi, std::string const & name )
{
	evi_dir_orf.clear();
	evi_rev_orf.clear();

	map< string, vector<INTERVAL_EVI> >::const_iterator  evi_itr = evi.find(name);

	if ( evi_itr == evi.end() )
		return;

	// evidence for sequence exists

	vector<INTERVAL_EVI>::const_iterator  itr     = evi_itr->second.begin();
	vector<INTERVAL_EVI>::const_iterator  itr_end = evi_itr->second.end();

	for( ; itr != itr_end; ++itr )
	{
		if      ( itr->key == EVIDENCE_TYPES::MASK )   SetMASK( itr->left, itr->right );
		else if ( itr->key == EVIDENCE_TYPES::NON )    SetNON( itr->left, itr->right );
		else if ( itr->key == EVIDENCE_TYPES::RDNA )   SetNoFullOverlap( itr->left, itr->right );
		else if ( itr->key == EVIDENCE_TYPES::TRNA )   SetNoFullOverlap( itr->left, itr->right );
		else if ( itr->key == EVIDENCE_TYPES::NRNA )   SetNoFullOverlap( itr->left, itr->right );
		else if ( itr->key == EVIDENCE_TYPES::COD )    SetCOD( *itr );
		else if ( itr->key == EVIDENCE_TYPES::CDS )    SetCDS( *itr );
	}
}
// ----------------------------------------------------
void Data::SetCOD( INTERVAL_EVI const & d )
{
	int L = d.left;
	int R = d.right;

	// move to index ; left and right inclusive
	--L;
	--R;

	if ( d.strand == STRAND_TYPES::DIRECT )
	{
		// find stop codon
		int stop;
		for( stop = R; stop < static_cast<int>( flag.size() ); stop += 3 )
		{
			if ( flag[stop] & dirEnd ) break;
		}

		// boundary
		if ( stop >= static_cast<int>( flag.size() ) )
			stop -= 3;

		// check for internal stops
		for( int i = stop - 3; i >= L+2; i -= 3 )
		{
			if ( flag[i] & dirEnd )
			{
				if ( logger->verbose )
					logger->Print( 0, "warning, internal stop codon was found in CDS evidence, record ignored:", d.left, d.right );
				return;
			}
		}

		// find minimum start
		int start = L + 2;

		// is boundary a start ?
		if ( flag[start] & dirStart )
			;
		else
		{
			// search for first upstream start
			for( int i = L - 1; i >= 2; i -= 3 )
			{
				if ( flag[i] & dirEnd )
				{
					if ( logger->verbose )
						logger->Print( 0, "warning, valid start is missing for CDS evidence, record ignored:", d.left, d.right );
					return;
				}

				if ( flag[i] & dirStart )
				{
					start = i;
					break;
				}
			}
		}

		// boundary
		if ( start < 2 )
			start += 3;

		// remove all internal starts
		for( int i = start + 3; i < stop; i += 3 )
		{
			if (flag[i] & dirStart ) flag[i] ^= dirStart;
		}

		// check if stop interval exists 

		map< int, INTERVAL_EVI >::iterator itr = evi_dir_orf.find( stop ); 

		if ( itr == evi_dir_orf.end() )
		{
			INTERVAL_EVI current_evi = d;

			current_evi.left = start - 2 + 1;
			current_evi.right = stop + 1;

			evi_dir_orf.insert( std::pair<int, INTERVAL_EVI>( stop, current_evi ) );
		}
		else
		{
			if ( itr->second.left > start - 2 + 1 )
				itr->second.left = start - 2 + 1;

			itr->second.attr += ( " " + d.attr );
		}
	}
	else if ( d.strand == STRAND_TYPES::REVERSE )
	{
		// find stop codon
		int stop;
		for( stop = L + 2; stop >= 2; stop -= 3 )
		{
			if ( flag[stop] & revEnd ) break;
		}

		// boundary
		if ( stop < 2 )
			stop += 3;

		// check for internal stops
		for( int i = stop + 3; i <= R; i += 3 )
		{
			if ( flag[i] & revEnd )
			{
				if ( logger->verbose )
					logger->Print( 0, "warning, internal stop codon was found in CDS evidence, record ignored:", d.left, d.right );
				return;
			}
		}

		// find minimum start
		int start = R;

		// is boundary a start
		if ( flag[start] & revStart )
			;
		else
		{
			// search for first upstream start
			for( int i = R + 3; i < static_cast<int>( flag.size() ); i += 3)
			{
				if ( flag[i] & revEnd )
				{
					if ( logger->verbose )
						logger->Print( 0, "warning, valid start is missing for CDS evidence, record ignored:", d.left, d.right );
					return;
				}

				if ( flag[i] & revStart )
				{
					start = i;
					break;
				}
			}
		}

		if ( start >= static_cast<int>( flag.size() ) )
			start -= 3;

		// remove all internal starts

		for( int i = start - 3; i > stop; i -= 3)
		{
			if (flag[i] & revStart ) flag[i] ^= revStart;
		}

		// check if interval exists
		std::map<int, INTERVAL_EVI >::iterator itr = evi_rev_orf.find( stop - 2 ); 

		if ( itr == evi_rev_orf.end() )
		{
			INTERVAL_EVI current_evi = d;

			current_evi.left = stop - 2 + 1;
			current_evi.right = start + 1;

			evi_rev_orf.insert( std::pair<int, INTERVAL_EVI>( stop - 2, current_evi ) );
		}
		else
		{
			if ( itr->second.right < start + 1 )
				itr->second.right = start + 1;

			itr->second.attr += ( " " + d.attr );
		}
	}
}
// ----------------------------------------------------
void Data::SetCDS( INTERVAL_EVI const & d )
{
	;
}
// ----------------------------------------------------
void Data::SetNON( int L, int R )
{
	// move to index ; left and right inclusive
	--L;
	--R;

	// remove starts inside the L..R

	int end = R + 2;
	if ( static_cast<int>( flag.size() ) - 1 < R + 2 )
		end = flag.size() - 1;

	for( int i = L; i <= end; ++i )
	{
		if( flag[i] & dirStart ) flag[i] ^= dirStart;
		if( flag[i] & revStart ) flag[i] ^= revStart;
	}

	// remove dir starts on Left side, till stop found

	for( int j = 0; j < 3; ++j )
	{
		for( int i = L - 1 - j; i >= 0; i -= 3 )
		{
			if( flag[i] & dirEnd ) 
				break;

			if( flag[i] & dirStart )
				flag[i] ^= dirStart;
		}
	}

	// remove rev starts on Right side, till stop found

	for( int j = 0; j < 3; ++j )
	{
		for( int i = R + 3 + j; i < static_cast<int>( flag.size() ); i += 3 )
		{
			if( flag[i] & revEnd ) 
				break;

			if( flag[i] & revStart )
				flag[i] ^= revStart;
		}
	}
}
// ----------------------------------------------------
void Data::SetMASK( int L, int R )
{
	// move to index ; left and right inclusive
	--L;
	--R;

	// remove starts inside the L..R

	int end = R + 2;
	if ( static_cast<int>( flag.size() ) - 1 < R + 2 )
		end = flag.size() - 1;

	for( int i = L; i <= end; ++i )
	{
		if( flag[i] & dirStart ) flag[i] ^= dirStart;
		if( flag[i] & revStart ) flag[i] ^= revStart;
	}
}
// ----------------------------------------------------
void Data::SetNoFullOverlap( int L, int R )
{
	// move to index ; left and right inclusive
	--L;
	--R;

	// remove starts of fully overlapping ORFs

	int end = R + 2;
	if ( static_cast<int>( flag.size() ) - 1 < R + 2 )
		end = flag.size() - 1;

	for( int i = L; i <= end; ++i )
	{
		if( flag[i] & dirEnd )
		{
			for( int j = i - 3; j >= L; j -= 3 )
			{
				if( flag[j] & dirEnd ) 
					break;

				if ( flag[j] & dirStart )
					flag[j] ^= dirStart;
			}
		}

		if( flag[i] & revEnd )
		{
			for( int j = i + 3; j <= end; j += 3 )
			{
				if( flag[j] & revEnd )
					break;

				if ( flag[j] & revStart )
					flag[j] ^= revStart;
			}
		}
	}
}
// ----------------------------------------------------

