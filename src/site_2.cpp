// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2016
// File: site.cpp
// Project: GeneMark.hmm-2 (no introns)
// Version:
//   0.1
// Tested with boost 1.48
// ====================================================

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <string>
#include <vector>
#include <iostream>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

#include "site_2.h"
#include "common_2.h"
#include "logger.h"
#include "index_2.h"

// ----------------------------------------------------
Site::Site()
{
	is_valid = false;

	width = 0;
	order = 0;
	margin = 0;
	max_dur = 0;
	height = 0;

	// unknown letters are penelized toward the noncoding states
	unknow_letter_logodd = log(0.4/0.6);

	max_dur = 0;
	norm_for_duration = 0;

	verbose = 0;
	debug = 0;
	log_here = 0;
}
// ----------------------------------------------------
void Site::ReserveSpace(void)
{
	// zero order background is used in current implementation
	
	non.assign( 4, 0 );

	height = 4 << (2*order);

	// two arrays are allocated for development state ; one array is enogth for release code	

	matrix.assign( height, vector<double>(width, 0) );
	logodd.assign( height, vector<double>(width, 0) );

	// duration may be empty : if max_dur == 0

	duration.assign( max_dur, 0 );
}
// ----------------------------------------------------
void Site::NormalizeArr( std::vector< std::vector<double> > & arr, string const message )
{
	unsigned int h = arr.size();
	unsigned int w = 0;

	for( unsigned int i = 0; i < h; ++i )
	{
		if ( w == 0 )
			w = arr[i].size();
		else
			assert( w == arr[i].size() );
	}

	for( unsigned int j = 0; j < w; ++j )
	{
		double sum = 0;

		for( unsigned int i = 0; i < h; ++i )
			sum += matrix[i][j];

		// lower orders are part of high order matrix, therefore "sum" is x4

		for( unsigned int k = j; k < order; ++k )
			sum /= 4.0;

		if ( sum < 0.99 || sum > 1.01  )
		{
			if (debug) cerr << "warning, probability array norm is not of expected range 0.99-1.01 : " << message << " " << sum << endl;
			if (debug) log_here->Print( 1, "warning, probability array norm is not of expected range 0.99-1.01 : " , message );
		}
	}
}
// ----------------------------------------------------
void Site::CheckForZeros(std::vector<double> const & arr)
{
	for( vector<double>::const_iterator itr = arr.begin(); itr != arr.end(); ++itr )
	{
		if ( *itr == 0 )
		{
			cerr << "error, unexpected zero probability value was found" << endl;
			exit(1);
		}
	}
}
// ----------------------------------------------------
void Site::Initialize( std::string const message, std::vector<double> const & arr )
{
	NormalizeArr( matrix, message );
	AbsoluteToConditionalByLastPosition(matrix);
	non = arr;
	CheckForZeros(non);
	CalcLogOdd();
	ClearSitePositions();
}
// ----------------------------------------------------
void Site::InitializeDuration( std::string const message, double const p )
{
	norm_for_duration = -1.0 / p;

	double sum = 0;

	for( vector<double>::iterator itr = duration.begin(); itr != duration.end(); ++itr )
	{
		sum += *itr;
	}

	if ( sum < 0.99 || sum > 1.01  )
	{
		if (debug) cerr << "warning, probability array norm is not of expected range 0.99-1.01 : " << message << " " << sum << endl;
		if (log_here) log_here->Print( 1, "warning, probability array norm is not of expected range 0.99-1.01 : " , message );
	}

	assert(sum);

	for( vector<double>::iterator itr = duration.begin(); itr != duration.end(); ++itr )
	{
		*itr /= sum;
	}

	duration_logodd = duration;

	for ( unsigned int i = 0; i < duration_logodd.size(); ++i )
	{
		if ( duration_logodd[i] == 0 )
		{
			duration_logodd[i] = LOG_ZERO;
		}
		else
		{
			duration_logodd[i] = log( duration_logodd[i] ) - (width + i)*norm_for_duration;
//            duration_logodd[i] = log( duration_logodd[i] ) - log(-norm_for_duration);
		}
	}
}
// ----------------------------------------------------
void Site::CalcLogOdd(void)
{
	unsigned int non_order = 0;
	if ( non.size() == 64 )
	{
		non_order = 2;
		AbsoluteToConditionalByLastPosition(non);
	}
	else if ( non.size() == 4 )
		non_order = 0;
	else
	{
		cerr << "error, this order of noncoding matrix is not supported yet : " << non.size() << endl;
		exit(1);
	}

	for( unsigned int i = 0; i < height; ++i )
	{
		for( unsigned int j = 0; j < width; ++j )
		{
			if ( matrix[i][j] )
			{
				if ( non_order == 0 )
					logodd[i][j] = log( matrix[i][j] / non[ i & 3 ] );
				else
					logodd[i][j] = log( matrix[i][j] / non[ i ] );
			}
			else
			{
				logodd[i][j] = LOG_ZERO;
			}
		}
	}
}
// ----------------------------------------------------
void Site::AbsoluteToConditionalByLastPosition( std::vector< std::vector<double> > & arr )
{
	for( unsigned int i = 0; i < height; i += 4 )
	{
		for( unsigned int j = 0; j < width; ++j )
		{
			double sum = arr[i+NT::A][j] + arr[i+NT::C][j] + arr[i+NT::G][j] + arr[i+NT::T][j];

			if ( sum != 0 )
			{
				arr[i+NT::A][j] /= sum;
				arr[i+NT::C][j] /= sum;
				arr[i+NT::G][j] /= sum;
				arr[i+NT::T][j] /= sum;
			}
			else
			{
				arr[i+NT::A][j] = 0;
				arr[i+NT::C][j] = 0;
				arr[i+NT::G][j] = 0;
				arr[i+NT::T][j] = 0;
			}
		}
	}
}
// ----------------------------------------------------
void Site::AbsoluteToConditionalByLastPosition( std::vector<double> & arr )
{
	for( unsigned int i = 0; i < arr.size(); i += 4 )
	{
		double sum = arr[i+NT::A] + arr[i+NT::C] + arr[i+NT::G] + arr[i+NT::T];

		if ( sum != 0 )
		{
			arr[i+NT::A] /= sum;
			arr[i+NT::C] /= sum;
			arr[i+NT::G] /= sum;
			arr[i+NT::T] /= sum;
		}
		else
		{
			arr[i+NT::A] = 0;
			arr[i+NT::C] = 0;
			arr[i+NT::G] = 0;
			arr[i+NT::T] = 0;
		}
	}
}
// ----------------------------------------------------
void Site::ClearSitePositions(void)
{
	if ( margin < 0 && abs(margin) < width + 3 )
	{
		unsigned int L = margin + width;
		unsigned int R = L + 3;

		for( unsigned int i = 0; i < height; ++i )
		{
			for( unsigned int j = L; j < R; ++j )
			{
				logodd[i][j] = 0;
			}
		}
	}
}
// ----------------------------------------------------
// vector of vectors
//                  position index
//                   0 1 2 3 ... j ... K   - width
//  letter index 0
//  letter index 1
//  letter index 2
//           ... i ...
//  letter index N
//           |
//         hight
//
//                logodd[i][j]   logodd[ letter index ][ position index ]
//
// ----------------------------------------------------
double Site::GetDir( std::vector<unsigned char> const & nt, unsigned int pos )
{
	double x = 0;

	// check that matrix is not outside of the sequence boundaries

	if ( width + margin <= pos )
	{
		pos = pos - margin - width;

		Index idx(order);

		for( unsigned int j = 0; j < width; ++j, ++pos )
		{
			unsigned int letter_index = idx.Get( nt[pos] );

			if ( letter_index != idx.npos )
			{
				x += logodd[ letter_index ][j];
			}
			else
			{
				x += unknow_letter_logodd;
			}
		}
	}
	else
	{
		int start_pos = margin + width - pos;

		x += ( start_pos * unknow_letter_logodd );

		Index idx(order);

		for( unsigned int j = start_pos; j < width; ++j )
		{
			unsigned int letter_index = idx.Get( nt[j] );

			if ( letter_index != idx.npos )
			{
				x += logodd[ letter_index ][j];
			}
			else
			{
				x += unknow_letter_logodd;
			}
		}	
	}

	return x;
}
// ----------------------------------------------------
double Site::GetRev( std::vector<unsigned char> const & nt, unsigned int pos )
{
	double x = 0;

	if ( width + margin + pos < nt.size() )
	{
		pos = pos + margin + width;

		Index idx(order);

		for( unsigned int j = 0; j < width; ++j, --pos )
		{
			unsigned int letter_index = idx.GetRC( nt[pos] );

			if ( letter_index != idx.npos )
			{
				x += logodd[ letter_index ][j];
			}
			else
			{
				x += unknow_letter_logodd;
			}
		}
	}
	else
	{
		int start_pos = margin + width + pos + 1 - nt.size();

		x += ( start_pos * unknow_letter_logodd );

		Index idx(order);

		for( unsigned int j = start_pos; j < width; ++j )
		{
			unsigned int letter_index = idx.GetRC( nt[nt.size() - width + j] );

			if ( letter_index != idx.npos)
			{
				x += logodd[ letter_index ][j];
			}
			else
			{
				x += unknow_letter_logodd;
			}
		}	
	}

	return x;
}
// ----------------------------------------------------
double Site::GetDirWithDur( std::vector<unsigned char> const & nt, unsigned int pos )
{
	double best = LOG_ZERO;
	unsigned int best_pos = -1;

	unsigned int max_pos = 0;

	if ( width + max_dur <= pos )
		max_pos = max_dur;
	else if ( pos >= width )
		max_pos = pos - width + 1;

	Index idx(order);

	for( unsigned int k = 0; k < max_pos; ++k )
	{
		idx.Reset();
		double x = 0;

		unsigned int current = pos - width - k;

		for( unsigned int j = 0; j < width; ++j, ++current )
		{
			unsigned int letter_index = idx.Get( nt[current] );

			if ( letter_index != idx.npos )
				x += logodd[ letter_index ][j];
			else
				x += unknow_letter_logodd;
		}
//
//        if (pos == 336) {
//            std::cout << k << " " << x << " " << duration_logodd[k] << std::endl;
//        }
		x += duration_logodd[k];

		if ( x > best )
		{
			best = x;
			best_pos = k;
		}
        
        
	}

	return best;
}
// ----------------------------------------------------
double Site::GetRevWithDur( std::vector<unsigned char> const & nt, unsigned int pos )
{
	double best = LOG_ZERO;
	unsigned int best_pos = 0;

	unsigned int max_pos = 0;

	if ( pos + width + max_dur < nt.size() )
		max_pos = max_dur;
	else if ( pos + width < nt.size() )
		max_pos = nt.size() - pos - width;

	Index idx(order);

	for( unsigned int k = 0; k < max_pos; ++k )
	{
		idx.Reset();
		double x = 0;

		unsigned int current = pos + width + k;

		for( unsigned int j = 0; j < width; ++j, --current )
		{
			unsigned int letter_index = idx.GetRC( nt[current] );

			if ( letter_index != idx.npos)
				x += logodd[ letter_index ][j];
			else
				x += unknow_letter_logodd;
		}

		x += duration_logodd[k];

		if ( x > best )
		{
			best = x;
			best_pos = k;
		}
	}

	return best;
}
// ----------------------------------------------------
SiteInfo Site::GetDirWithDurFullInfo( std::vector<unsigned char> const & nt, unsigned int pos )
{
	SiteInfo si;

	double best = LOG_ZERO;
	unsigned int best_pos = -1;

	unsigned int max_pos = 0;

	if ( width + max_dur <= pos )
		max_pos = max_dur;
	else if ( pos >= width )
		max_pos = pos - width + 1;

	Index idx(order);

	for( unsigned int k = 0; k < max_pos; ++k )
	{
		idx.Reset();
		double x = 0;

		unsigned int current = pos - width - k;

		for( unsigned int j = 0; j < width; ++j, ++current )
		{
			unsigned int letter_index = idx.Get( nt[current] );

			if ( letter_index != idx.npos )
				x += logodd[ letter_index ][j];
			else
				x += unknow_letter_logodd;
		}

		x += duration_logodd[k];

		if ( x > best )
		{
			best = x;
			best_pos = k;
		}
	}

	if ( best != LOG_ZERO )
	{
		for( unsigned int i = 0; i < width; ++i )
		{
			si.seq.push_back( LETTER[ nt[ pos - width - best_pos + i ] ] );
		}
	}

	si.score = best;
	si.spacer = best_pos;

	return si;
};
// ----------------------------------------------------
SiteInfo Site::GetRevWithDurFullInfo( std::vector<unsigned char> const & nt, unsigned int pos )
{
	SiteInfo si;

	double best = LOG_ZERO;
	unsigned int best_pos = 0;

	unsigned int max_pos = 0;

	if ( pos + width + max_dur < nt.size() )
		max_pos = max_dur;
	else if ( pos + width < nt.size() )
		max_pos = nt.size() - pos - width;

	Index idx(order);

	for( unsigned int k = 0; k < max_pos; ++k )
	{
		idx.Reset();
		double x = 0;

		unsigned int current = pos + width + k;

		for( unsigned int j = 0; j < width; ++j, --current )
		{
			unsigned int letter_index = idx.GetRC( nt[current] );

			if ( letter_index != idx.npos )
				x += logodd[ letter_index ][j];
			else
				x += unknow_letter_logodd;
		}

		x += duration_logodd[k];

		if ( x > best )
		{
			best = x;
			best_pos = k;
		}
	}

	if ( best != LOG_ZERO )
	{
		for( unsigned int i = 0; i < width; ++i )
		{
			si.seq.push_back( RC_LETTER[ nt[ pos + best_pos + width - i ] ] );
		}
	}

	si.score = best;
	si.spacer = best_pos;

	return si;
}
// ----------------------------------------------------
double Site::Get( std::vector<unsigned char> const & nt, unsigned int pos, unsigned int status )
{
	if ( status & dirStart )  return GetDir( nt, pos );
	if ( status & revStart )  return GetRev( nt, pos );

	return LOG_ZERO;
}
// ----------------------------------------------------
double Site::GetWithDur( std::vector<unsigned char> const & nt, unsigned int pos, unsigned int status )
{
	if ( status & dirStart )  return GetDirWithDur( nt, pos );
	if ( status & revStart )  return GetRevWithDur( nt, pos );

	return LOG_ZERO;
}
// ----------------------------------------------------
double Site::GetMaxDurationScore() const {
    double max = LOG_ZERO;
    for (size_t i = 0; i < duration_logodd.size(); i++) {
        if (duration_logodd[i] > max)
            max = duration_logodd[i];
    }
    return max;
}

