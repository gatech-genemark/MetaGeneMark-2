// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Project: GeneMark.hmm-2 (no introns)
// File: model_2.cpp
// ====================================================

#include <cstdlib>
#include <cmath>
#include <cassert>

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <sstream>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::stringstream;
using std::distance;

#include "model_2.h"
#include "common_2.h"
#include "exit.h"

// ----------------------------------------------------
void Model::Initialize(void)
{
	if (logger->debug)
	{
		VerifyProbabilityArray( cod1, order_cod, "cod1" );
		VerifyProbabilityArray( cod2, order_cod, "cod2" );
		VerifyProbabilityArray( cod3, order_cod, "cod3" );
		VerifyProbabilityArray( non,  order_non, "non" );
	}

	NormalizeArray( cod1, "cod1" );
	NormalizeArray( cod2, "cod2" );
	NormalizeArray( cod3, "cod3" );
	NormalizeArray( non,  "non" );

	// P[i]=0 is not allowed in non-coding arrays

	CheckForNoZeroValues( non, order_non, "non" );

	CalculateLogOdd();
	CalculateLogProb();

	CalculateStartStopLoggodd();
	
	InitiateLogoddDuration( RESERVED_DURATION );
	InitiateLogoddIncompleDuration( RESERVED_INCOMPLETE_DURATION );

	toModel = LogRatio( toModel, 1 );

	InitializeSite( &RBS, true, 0, "RBS" );
	InitializeSite( &StartContentRBS, false, 2, "StartContentRBS" );
	InitializeSite( &Promoter, true, 0, "Promoter" );
	InitializeSite( &StartContentPromoter, false, 2, "StartContentPromoter" );
	InitializeSite( &StartContent, false, 2, "StartContent" );
    
    InitializeMultiShiftSite( &RBS_A, true, 0, "RBS_A" );
    InitializeMultiShiftSite( &RBS_B, true, 0, "RBS_B" );
    InitializeMultiShiftSite( &RBS_C, true, 0, "RBS_C" );
    InitializeMultiShiftSite( &RBS_D, true, 0, "RBS_D" );
    InitializeMultiShiftSite( &RBS_X, true, 0, "RBS_X" );
    
    InitializeSite( &SC_RBS_A, false, 2, "SC_RBS_A" );
    InitializeSite( &SC_RBS_B, false, 2, "SC_RBS_B" );
    InitializeSite( &SC_RBS_C, false, 2, "SC_RBS_C" );
    InitializeSite( &SC_RBS_D, false, 2, "SC_RBS_D" );
    InitializeSite( &SC_RBS_X, false, 2, "SC_RBS_X" );
    
//    InitializeSite( &PROMOTER_C, true, 0, "Promoter_A" );
//    InitializeSite( &PROMOTER_D, true, 0, "Promoter_B" );
    InitializeMultiShiftSite( &PROMOTER_C, true, 0, "PROMOTER_C" );
    InitializeMultiShiftSite( &PROMOTER_D, true, 0, "PROMOTER_D" );
    
    
    InitializeSite( &SC_PROMOTER_C, false, 2, "SC_PROMOTER_C" );
    InitializeSite( &SC_PROMOTER_D, false, 2, "SC_PROMOTER_D" );

	if (logger->debug > 10)
	{
		logger->Print( PrintLogOdds() );
		logger->Print( PrintLogOdds_abs() );
		logger->Print( PrintDuration() );
	}
}
// ----------------------------------------------------
Model::Model( Logger * const logger ) : logger(logger)
{
	gcode = 0;
	gene_min_length = 0;
	order_cod = 0;
	order_non = 0;

	noncoding_duration_decay = 0;
	coding_duration_decay = 0;

	probability_N_in_noncoding = 0;
	probability_N_in_coding = 0;

	pATG = 0;
	pGTG = 0;
	pTTG = 0;

	pTAA = 0;
	pTAG = 0;
	pTGA = 0;
    
    // mgm starts/stops
    pATG_A = 0; pATG_B = 0; pATG_C = 0; pATG_D = 0; pATG_X = 0;
    pGTG_A = 0; pGTG_B = 0; pGTG_C = 0; pGTG_D = 0; pGTG_X = 0;
    pTTG_A = 0; pTTG_B = 0; pTTG_C = 0; pTTG_D = 0; pTTG_X = 0;
    pTAA_A = 0; pTAA_B = 0; pTAA_C = 0; pTAA_D = 0; pTAA_X = 0;
    pTGA_A = 0; pTGA_B = 0; pTGA_C = 0; pTGA_D = 0; pTGA_X = 0;
    pTAG_A = 0; pTAG_B = 0; pTAG_C = 0; pTAG_D = 0; pTAG_X = 0;

	logodd_ATG = 0;
	logodd_GTG = 0;
	logodd_TTG = 0;

	logodd_TAA = 0;
	logodd_TAG = 0;
	logodd_TGA = 0;

	toModel = 0;

	ORF_start_marging = 3;

	dur_c = 0;
	dur_r = 0;
}
// ----------------------------------------------------
void Model::ReserveSpace(void)
{
	// 4^(order+1) ; 4^(0+1) = 4 ; 4^(1+1) = 16
	//  bin 0000 0100 << 0 = 0000 0100 =  4
	//  bin 0000 0100 << 2 = 0001 0000 = 16

	unsigned int cod_array_size = ( 4<<(2*order_cod) );
	unsigned int non_array_size = ( 4<<(2*order_non) );
	
	cod1.assign( cod_array_size, DEFVALUES::ERROR );
	cod2.assign( cod_array_size, DEFVALUES::ERROR );
	cod3.assign( cod_array_size, DEFVALUES::ERROR );
	non.assign ( non_array_size, DEFVALUES::ERROR );

	// reserve and initialise for size + 1 ; use "+1" location for probabilities of k-mers with letters 'N' in  

	unsigned int logodd_array_size = ( order_cod > order_non ) ? cod_array_size : non_array_size;

	// logodds

	logodd_1_abs.assign( logodd_array_size + 1, DEFVALUES::ZERO );
	logodd_2_abs.assign( logodd_array_size + 1, DEFVALUES::ZERO );
	logodd_3_abs.assign( logodd_array_size + 1, DEFVALUES::ZERO );

	logodd_1.assign( logodd_array_size + 1, DEFVALUES::ZERO );
	logodd_2.assign( logodd_array_size + 1, DEFVALUES::ZERO );
	logodd_3.assign( logodd_array_size + 1, DEFVALUES::ZERO );

	// log normalized probabilities

	logP_1_abs_n.assign( cod_array_size + 1, DEFVALUES::ZERO );
	logP_2_abs_n.assign( cod_array_size + 1, DEFVALUES::ZERO );
	logP_3_abs_n.assign( cod_array_size + 1, DEFVALUES::ZERO );
	
	logP_1_n.assign( cod_array_size + 1, DEFVALUES::ZERO );
	logP_2_n.assign( cod_array_size + 1, DEFVALUES::ZERO );
	logP_3_n.assign( cod_array_size + 1, DEFVALUES::ZERO );

	logP_N_abs_n.assign( non_array_size + 1, DEFVALUES::ZERO );
	logP_N_n.assign    ( non_array_size + 1, DEFVALUES::ZERO );
}
// ----------------------------------------------------
void Model::CalculateLogOdd(void)
{
	vector<double> non_conditional( (4<<(2*order_non)), DEFVALUES::ERROR );
	
	AbsoluteToConditionalByLastPosition( non, non_conditional );

	if ( order_cod == order_non )
	{
		CalculateAbsoluteCodVsNonRatio( cod1, non, logodd_1_abs );
		CalculateAbsoluteCodVsNonRatio( cod2, non, logodd_2_abs );
		CalculateAbsoluteCodVsNonRatio( cod3, non, logodd_3_abs );
	}
	else if ( order_cod > order_non )
	{
		CalculateAbsoluteCodVsNonRatio( cod1, non, non_conditional, logodd_1_abs );
		CalculateAbsoluteCodVsNonRatio( cod2, non, non_conditional, logodd_2_abs );
		CalculateAbsoluteCodVsNonRatio( cod3, non, non_conditional, logodd_3_abs );
	}
	else
		assert(0); // todo

	AbsoluteToConditionalByLastPosition( cod1, logodd_1 );
	AbsoluteToConditionalByLastPosition( cod2, logodd_2 );
	AbsoluteToConditionalByLastPosition( cod3, logodd_3 );

	DivideByConditionalNon( non_conditional, logodd_1 );
	DivideByConditionalNon( non_conditional, logodd_2 );
	DivideByConditionalNon( non_conditional, logodd_3 );

	Log(logodd_1_abs);
	Log(logodd_2_abs);
	Log(logodd_3_abs);

	Log(logodd_1);
	Log(logodd_2);
	Log(logodd_3);

	SetLogOddForN( logodd_1_abs, LogRatio( probability_N_in_coding, probability_N_in_noncoding) );
	SetLogOddForN( logodd_2_abs, LogRatio( probability_N_in_coding, probability_N_in_noncoding) );
	SetLogOddForN( logodd_3_abs, LogRatio( probability_N_in_coding, probability_N_in_noncoding) );

	SetLogOddForN( logodd_1, LogRatio( probability_N_in_coding, probability_N_in_noncoding) );
	SetLogOddForN( logodd_2, LogRatio( probability_N_in_coding, probability_N_in_noncoding) );
	SetLogOddForN( logodd_3, LogRatio( probability_N_in_coding, probability_N_in_noncoding) );
}
// ----------------------------------------------------
void Model::CalculateLogProb(void)
{
	vector<double> norm( 4, 0.25 );

	if ( order_cod == 0 )
	{
		CalculateAbsoluteCodVsNonRatio( cod1, norm, logP_1_abs_n );
		CalculateAbsoluteCodVsNonRatio( cod2, norm, logP_2_abs_n );
		CalculateAbsoluteCodVsNonRatio( cod3, norm, logP_3_abs_n );
	}
	else // order_cod > 0
	{
		CalculateAbsoluteCodVsNonRatio( cod1, norm, norm, logP_1_abs_n );
		CalculateAbsoluteCodVsNonRatio( cod2, norm, norm, logP_2_abs_n );
		CalculateAbsoluteCodVsNonRatio( cod3, norm, norm, logP_3_abs_n );
	}

	if ( order_non == 0 )
		CalculateAbsoluteCodVsNonRatio( non, norm, logP_N_abs_n );
	else // order_non > 0
		CalculateAbsoluteCodVsNonRatio( non, norm, norm, logP_N_abs_n );

	AbsoluteToConditionalByLastPosition( cod1, logP_1_n );
	AbsoluteToConditionalByLastPosition( cod2, logP_2_n );
	AbsoluteToConditionalByLastPosition( cod3, logP_3_n );
	AbsoluteToConditionalByLastPosition( non,  logP_N_n );

	DivideByConditionalNon( norm, logP_1_n );
	DivideByConditionalNon( norm, logP_2_n );
	DivideByConditionalNon( norm, logP_3_n );
	DivideByConditionalNon( norm, logP_N_n );

	Log(logP_1_abs_n);
	Log(logP_2_abs_n);
	Log(logP_3_abs_n);
	Log(logP_N_abs_n);

	Log(logP_1_n);
	Log(logP_2_n);
	Log(logP_3_n);
	Log(logP_N_n);

	SetLogOddForN( logP_1_abs_n, LogRatio(probability_N_in_coding, 0.25) );
	SetLogOddForN( logP_2_abs_n, LogRatio(probability_N_in_coding, 0.25) );
	SetLogOddForN( logP_3_abs_n, LogRatio(probability_N_in_coding, 0.25) );
	SetLogOddForN( logP_N_abs_n, LogRatio(probability_N_in_noncoding, 0.25) );

	SetLogOddForN( logP_1_n, LogRatio(probability_N_in_coding, 0.25) );
	SetLogOddForN( logP_2_n, LogRatio(probability_N_in_coding, 0.25) );
	SetLogOddForN( logP_3_n, LogRatio(probability_N_in_coding, 0.25) );
	SetLogOddForN( logP_N_n, LogRatio(probability_N_in_noncoding, 0.25) );
}
// ----------------------------------------------------
void Model::CalculateStartStopLoggodd(void)
{
	vector<double> &non_2 = this->non_2;
	vector<double> non_0;
	
	non_2.reserve(64);
	non_0.reserve(4);

	if (order_non == 2)
	{
		non_2 = non;
	}
	else if ( order_non > 2 )
	{
		non_2 = ReduceOrderAbs( non, order_non, 2);
	}
	else
	{
		if ( order_non == 0 )
			non_0 = non;
		else
			non_0 = ReduceOrderAbs( non, order_non, 0);

		for( int i = 0; i < 4; ++i )
		{
			for( int j = 0; j < 4; ++j )
			{
				for( int k = 0; k < 4; ++k )
				{
					non_2.push_back( non_0[i]*non_0[j]*non_0[k] );
				}
			}
		}
	}

	logodd_ATG = LogRatio( pATG, non_2[NT::A<<4|NT::T<<2|NT::G] );
	logodd_GTG = LogRatio( pGTG, non_2[NT::G<<4|NT::T<<2|NT::G] );
	logodd_TTG = LogRatio( pTTG, non_2[NT::T<<4|NT::T<<2|NT::G] );

	logodd_TAA = LogRatio( pTAA, non_2[NT::T<<4|NT::A<<2|NT::A] );
	logodd_TAG = LogRatio( pTAG, non_2[NT::T<<4|NT::A<<2|NT::G] );
	logodd_TGA = LogRatio( pTGA, non_2[NT::T<<4|NT::G<<2|NT::A] );
}
// ----------------------------------------------------
void Model::VerifyProbabilityArray( std::vector<double> const & arr, unsigned int const order, std::string const & message )
{
	vector<double>::const_iterator  itr     = arr.begin();
	vector<double>::const_iterator  itr_end = arr.end();

	for( ; itr != itr_end; ++itr )
	{
		if ( *itr == DEFVALUES::ERROR )
		{
			string label = IntToString( distance( arr.begin(), itr ), order );

			Exit( "error, default initializaion value was detected in array for index:", message, label, *itr );
		}
		
		if ( *itr < 0 || *itr > 1 )
		{
			string label = IntToString( distance( arr.begin(), itr ), order );

			Exit( "error, probability value out of range 0..1 in array for index:", message, label, *itr );
		}
	}
}
// ----------------------------------------------------
void Model::NormalizeArray( std::vector<double> & arr, std::string const & message )
{
	double sum = 0;

	for( vector<double>::iterator itr = arr.begin(); itr != arr.end(); ++itr )
		sum += *itr;

	if ( sum < 0.99 || sum > 1.01  )
		logger->Print( 1, "warning, sum of probability values is out of the range 0.99-1.01:", message, sum );

	if ( !sum )
		Exit( "error, sum of probability values equals zero", message );

	for( vector<double>::iterator itr = arr.begin(); itr != arr.end(); ++itr )
		*itr /= sum;
}
// ----------------------------------------------------
void Model::CheckForNoZeroValues( std::vector<double> const & arr, unsigned int const order, std::string const & message )
{
	for( vector<double>::const_iterator itr = arr.begin(); itr != arr.end(); ++itr )
	{
		if ( *itr == 0 )
		{
			string label = IntToString( distance( arr.begin(), itr ), order );

			Exit( "error, probability zero found in array at position:", message, label );
		}
	}
}
// ----------------------------------------------------
void Model::AbsoluteToConditionalByLastPosition( std::vector<double> const & source, std::vector<double> & target )
{
	if ( source.size() > target.size() ) Exit("error in 'Model::AbsoluteToConditionalByLastPosition'");

	for( unsigned int i = 0; i < source.size(); i += 4 )
	{
		double sum = source[ i + NT::A ] + source[ i + NT::C ] + source[ i + NT::G ] + source[ i + NT::T ];

		if ( sum != 0 )
		{
			target[ i + NT::A ] = source[ i + NT::A ] / sum;
			target[ i + NT::C ] = source[ i + NT::C ] / sum;
			target[ i + NT::G ] = source[ i + NT::G ] / sum;
			target[ i + NT::T ] = source[ i + NT::T ] / sum;
		}
		else
		{
			target[ i + NT::A ] = 0;
			target[ i + NT::C ] = 0;
			target[ i + NT::G ] = 0;
			target[ i + NT::T ] = 0;
		}
	}
}
// ----------------------------------------------------
void Model::CalculateAbsoluteCodVsNonRatio( std::vector<double> const & c, std::vector<double> const & n, std::vector<double> & target )
{
	if ( c.size() != n.size() )   Exit("error in 'Model::CalculateAbsoluteCodVsNonRatio'");

	for( unsigned int i = 0; i < c.size(); ++i )
	{
		target[i] = c[i] / n[i];
	}
}
// ----------------------------------------------------
void Model::CalculateAbsoluteCodVsNonRatio( std::vector<double> const & c, std::vector<double> const & n, std::vector<double> const & nc, std::vector<double> & target )
{
	if ( c.size() < n.size() )   Exit("error in 'Model::CalculateAbsoluteCodVsNonRatio'");
	if ( n.size() != nc.size() ) Exit("error in 'Model::CalculateAbsoluteCodVsNonRatio'");
	
	int order_c = static_cast<int>( log( 1.0 * c.size())/log(4.0) ) - 1;
	int order_n = static_cast<int>( log( 1.0 * n.size())/log(4.0) ) - 1;

	// order_c = 3  : index build from four letters
	// order_n = 1  : index build from two  letters

	//  c:  P( w, x, y, z, )
	//  n:  Q( w, x ) * Q(y|x) * Q(z|y)

	// 'd' conditional probabilities are needed to make the difference between 'order_c' and 'order_n'

	// 'd' = 3 - 1 = 2  : Q(y|x) and Q(z|y)

	int const d = order_c - order_n;

	// 'mask' - helps to get lower order k-mer from highr order K-mer

	//  index( wx )   & mask == index( wx )
	//  index( wxy )  & mask == index( xy )
	//  index( wxyz ) & mask == index( yz )

	unsigned int const mask = (4<<(2*order_n)) - 1;

	// 'i' is index in higher order
	// 'j' is index in lower order
	
	for( unsigned int i = 0; i < c.size(); ++i )
	{
		unsigned int j = i;

		// index >> 4  ==  w x
		
		j >>= (2*d);

		// 2*'d' = 4 bits
		// i  -> j
		// 0  -> 0
		// 1  -> 0
		// ...
		// 15 -> 0
		// 16 -> 1

		// P( w, x, y, z, ) / Q( w, x )

		target[i] = c[i] / n[j];

		for( int r = d - 1; r >= 0 ; --r )
		{
			// 'k' is index in conditional vector of lower order

			unsigned int k = i;

			// 'r'= 1 <= 'd' - 1 = 2 - 1 :  index >> 4  ==  w x y
			// 'r'= 0 <=  r--            :  index >> 4  ==  w x y z

			k >>= (2*r);

			// 'r' = 1 => 'k' & mask : w x y   => x y
			// 'r' = 0 => 'k' & mask : w x y z => y z

			k &= mask;
		
			// 'r' = 1 =>  ~/Q(y|x)
			// 'r' = 0 =>  ~/Q(z|y)

			target[i] /= nc[k];
		}

		//  P( w, x, y, z, ) / Q( w, x ) / Q(y|x) / Q(z|y)
	}
}
// ----------------------------------------------------
void Model::DivideByConditionalNon( std::vector<double> const & nc, std::vector<double> & target )
{
	assert( target.size() >= nc.size() );

	int order_target = static_cast<int>( log(1.0 * target.size())/log(4.0) ) - 1;
	int order_nc     = static_cast<int>( log(1.0 * nc.size())    /log(4.0) ) - 1;

	unsigned int mask = (4<<(2*order_nc)) - 1;
	unsigned int max  = (4<<(2*order_target));

	for( unsigned int i = 0; i < max; ++i )
	{
		target[i] /= nc[ i & mask ];
	}
}
// ----------------------------------------------------
double Model::LogRatio( double x, double y )
{
	double z = 0;

	if ( x > 0 && y > 0 )
		z = log( x/y );
	else if ( x == 0 && y == 0 )
		z = 0;
	else if ( x == 0 && y != 0 )
		z = LOG_ZERO;
	else
		assert (x >= 0 && y > 0);

	return z;
}
// ----------------------------------------------------
void Model::Log( std::vector<double> & target )
{
	for( vector<double>::iterator itr = target.begin(); itr != target.end(); ++itr )
	{
		if ( *itr > 0 )
			*itr = log( *itr );
		else if ( *itr == 0 )
			*itr = LOG_ZERO;
		else
			assert (*itr >= 0);
	}
}
// ----------------------------------------------------
void Model::SetLogOddForN( std::vector<double> & target, double const p )
{
	if ( !p )
	{
		cerr << "error, 'log' of zero detected" << endl;
		exit(1);
	}

	target.back() = p;
}
// ----------------------------------------------------
std::string  Model::PrintLogOdds(void)
{
	stringstream s;

	s << "# logodd of conditional probabilities: " << model_name << endl;

	if ( logodd_1.size() )
	{
		unsigned int i;

		for( i = 0; i < logodd_1.size() - 1; ++i )
		{	
			s << i << " " << IntToString(i, order_cod) << " " << logodd_1[i] << " " << logodd_2[i] << " " << logodd_3[i] << endl;
		}

		s << i << " " << "N" << " " << logodd_1[i] << " " << logodd_2[i] << " " << logodd_3[i] << endl;
	}
	else
	{
		s << "array is empty" << endl;
	}

	return s.str();
}
// ----------------------------------------------------
std::string Model::PrintLogOdds_abs(void)
{
	stringstream s;

	s << "# logodd of absolute probabilities: " << model_name << endl;

	if ( logodd_1_abs.size() )
	{
		unsigned int i;

		for( i = 0; i < logodd_1_abs.size() - 1; ++i )
		{	
			s << i << " " << IntToString(i, order_cod) << " " << logodd_1_abs[i] << " " << logodd_2_abs[i] << " " << logodd_3_abs[i] << endl;
		}

		s << i << " " << "N" << " " << logodd_1_abs[i] << " " << logodd_2_abs[i] << " " << logodd_3_abs[i] << endl;
	}
	else
	{
		s << "array is empty" << endl;
	}

	return s.str();
}
// ----------------------------------------------------
std::string Model::IntToString( unsigned int index, unsigned int const order )
{
	string key( order + 1, ' ' );

	unsigned int const mask = 3;

	for( unsigned int i = 0; i <= order; ++i )
	{
		switch( index & mask )
		{
			case 0:
				key[order - i] = 'A';
				break;
			case 1:
				key[order - i] = 'C';
				break;
			case 2:
				key[order - i] = 'G';
				break;
			case 3:
				key[order - i] = 'T';
				break;
		}

		index >>= 2;
	}

	return key;
}


//////////////////////// BEGIN NEW DURATION
void Model::InitiateLogoddIncompleDuration(unsigned int const reserved_length)
{
    assert(noncoding_duration_decay);
    assert(coding_duration_decay);
    assert(gene_min_length > 0);

    logodd_incomplete_duration.resize( reserved_length, 0 );

    double q = exp( -1.0 / noncoding_duration_decay );
    double p = exp( -3.0 / coding_duration_decay );

    double A = (1.0 - p)*(1.0 - p)*(1.0 - p) / (1.0 + p);

    double total = 0;

    logodd_incomplete_duration[3] = 1 - A;
    total += logodd_incomplete_duration[3];

    for( unsigned int i = 4; i < reserved_length; ++i )
    {
        if ( i%3 != 0) continue;

        logodd_incomplete_duration[i] = logodd_incomplete_duration[i-3] - A*(i*i/9)*pow(p, i/3-1);
        total += logodd_incomplete_duration[i];
    }

    Log(logodd_incomplete_duration);

    total = log(total);

    for (unsigned int i = 1; i < reserved_length; ++i)
    {
        if (i % 3 != 0) continue;

        logodd_incomplete_duration[i] -= total;
        logodd_incomplete_duration[i] -= log(1-q) + (i-1)*log(q) + log(1/(q*q) + 1/q + 1);
    }
}

double Model::GetIncompleteDurationForORF(unsigned int const length)
{
    double x = 0;
    assert(!(length % 3));

    if (length < logodd_incomplete_duration.size())
        x = logodd_incomplete_duration[length];
    else
    {
        int size = logodd_incomplete_duration.size();
        int index = size - size%3;

        x = logodd_incomplete_duration[index];
    }

    return x;
}
//////////////////////// END NEW DURATION


// ----------------------------------------------------
void Model::InitiateLogoddDuration(unsigned int const reserved_length)
{
	assert(noncoding_duration_decay);
	assert(coding_duration_decay);
	assert(gene_min_length > 0);

	logodd_duration.resize( reserved_length, LOG_ZERO );

	double q = exp( -1.0 / noncoding_duration_decay );
	double p = exp( -3.0 / coding_duration_decay );

	dur_r = log(p)/3 - log(q);
	dur_c = log((1-p)*(1-p)*(1-p)/(p*(1+p))) + log((1-q)/q) - 2*log(3.0);

	for( unsigned int i = gene_min_length; i < reserved_length; ++i )
	{
		if ( i%3 != 0) continue;

		logodd_duration[i] = i*dur_r + 2*log((double)i) + dur_c;
	}
}
// ---------------------------------------------------- OLD
//void Model::InitiateLogoddIncompleDuration(unsigned int const reserved_length)
//{
//	assert(noncoding_duration_decay);
//	assert(coding_duration_decay);
//	assert(gene_min_length > 0);
//
//	logodd_incomplete_duration.resize( reserved_length, LOG_ZERO );
//
//	double q = exp( -1.0 / noncoding_duration_decay );
//	double p = exp( -3.0 / coding_duration_decay );
//
//	dur_r = log(p)/3 - log(q);
//	dur_c = log((1-p)*(1-p)*(1-p)/(p*(1+p))) + log((1-q)/q) - 2*log(3.0);
//
//	for( unsigned int i = 1; i < reserved_length; ++i )
//	{
//		if ( i%3 != 0) continue;
//
//		logodd_incomplete_duration[i] = i*dur_r + 2*log((double)i) + dur_c;
//	}
//}
// ----------------------------------------------------
double Model::GetDurationForORF(unsigned int const length)
{
	double x = 0;
	assert(!(length % 3));

	if ( length < logodd_duration.size() )
		x = logodd_duration[length];
	else
		x = length*dur_r + 2*log((double)length) + dur_c;

	return x;
}
// ---------------------------------------------------- Old
//double Model::GetIncompleteDurationForORF(unsigned int const length)
//{
//	double x = 0;
//	assert(!(length % 3));
//
//	if ( length < logodd_incomplete_duration.size() )
//		x = logodd_incomplete_duration[length];
//	else
//		x = length*dur_r + 2*log((double)length) + dur_c;
//
//	return x;
//}
// ----------------------------------------------------
std::string Model::PrintDuration(void)
{
	stringstream s;

	s << "# logodd of conditional probabilities: " << model_name << endl;

	for( unsigned int i = 0; i < logodd_duration.size() - 1; ++i )
	{	
		s << i << " " << logodd_duration[i] << endl;
	}

	return s.str();
}
// ----------------------------------------------------
std::vector<double> Model::ReduceOrderAbs( std::vector<double> const & arr , unsigned int const order_in, unsigned int const order_out )
{
	assert( order_in >= order_out );

	if ( order_in == order_out )
	{
		return arr;
	}

	unsigned int size = 4<<(2*order_out);

	vector<double> out( size, 0 );

	unsigned int step = order_in - order_out;
	unsigned int L = 0;
	unsigned int R = (4<<(2*(step-1)));

	for( unsigned int i = 0; i < size; ++i )
	{
		for( unsigned int j = L; j < R; ++j )
		{
			out[i] += arr[j];
		}

		L = R;
		R += (4<<(2*(step-1)));
	}

	ReverseCompNoncodingCounts( order_out, out );

	return out;
}
// ----------------------------------------------------
unsigned int Model::ReversCompIndex( unsigned int idx, unsigned int order )
{
	unsigned int new_idx = 0;

	for( unsigned int i = 0; i <= order; ++i )
	{
		new_idx <<= 2;
		new_idx += 3 - ( idx & 3 );
		idx >>= 2;
	}

	return new_idx;
}
// ----------------------------------------------------
void Model::ReverseCompNoncodingCounts( unsigned int order, std::vector<double> & arr )
{
	unsigned int idx;
	unsigned int idx_revcomp;

	double tmp_count;

	unsigned int size = 4<<(2*order);

	for( idx = 0; idx < size; ++idx )
	{
		idx_revcomp = ReversCompIndex( idx, order );

		tmp_count = (arr[idx] + arr[idx_revcomp])/2;

		arr[idx] = tmp_count;
		arr[idx_revcomp] = tmp_count;
	}
}
// ----------------------------------------------------
std::string Model::KL(void)
{
	// find minimum probability of coding in frames 1, 2 and 3

	double min = 1;
	int count = 0;

	for( vector<double>::iterator itr = cod1.begin(); itr != cod1.end(); ++itr )
	{
		if ( *itr < min )
		{
			min = *itr;
			count = 1;
		}
		else if ( *itr == min )
		{
			++count;
		}
	}

	for( vector<double>::iterator itr = cod2.begin(); itr != cod2.end(); ++itr )
	{
		if ( *itr < min )
		{
			min = *itr;
			count = 1;
		}
		else if ( *itr == min )
		{
			++count;
		}
	}

	for( vector<double>::iterator itr = cod3.begin(); itr != cod3.end(); ++itr )
	{
		if ( *itr < min )
		{
			min = *itr;
			count = 1;
		}
		else if ( *itr == min )
		{
			++count;
		}
	}

	// calculate KL

	double kl = 0;

	int count_ignored = 0;

	for( unsigned int i = 0 ; i < cod1.size(); ++i )
	{
		if ( logodd_1[i] == LOG_ZERO || ( cod1[i] == min && count > 1 ) )
			++count_ignored;
		else
			kl += cod1[i] * logodd_1[i];

		if ( logodd_2[i] == LOG_ZERO || ( cod2[i] == min && count > 1 ) )
			++count_ignored;
		else
			kl += cod2[i] * logodd_2[i];

		if ( logodd_3[i] == LOG_ZERO || ( cod3[i] == min && count > 1 ) )
			++count_ignored;
		else
			kl += cod3[i] * logodd_3[i];
	}

	stringstream s;

	s << "KL " << kl/3 << " min_prob " << min << " min_prob_count " << count << " prob_ignored " << count_ignored;

	return s.str();
}
// ----------------------------------------------------
void Model::InitializeSite( Site * ptr, bool with_dur, int norm_order, std::string const & label )
{
	if ( ptr->is_valid == false )
		return;
	
	if ( norm_order == 2 && ptr->order == 2 && order_non >= 2 )
	{
		ptr->Initialize( label, ReduceOrderAbs( non, order_non, 2 ) );
	}
	else // norm_order == 0
	{
		ptr->Initialize( label, ReduceOrderAbs( non, order_non, 0 ) );
	}

	if ( with_dur )
		ptr->InitializeDuration( label, noncoding_duration_decay );

	if (ptr->margin < -3 )
	{
		ORF_start_marging = abs( ptr->margin );
	}
}

void Model::InitializeMultiShiftSite( MultiShiftSite* ptr, bool with_dur, int norm_order, std::string const & label )
{
    if ( ptr->is_valid == false )
        return;
    
    vector<Site* > sites = ptr->get_sites();
    
    for (size_t i = 0; i < sites.size(); i++) {
        InitializeSite(sites[i], with_dur, norm_order, label);
    }
}
// ----------------------------------------------------

