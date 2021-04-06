// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2016
// File: site.h
// Project: GeneMark.hmm-2 (no introns)
// Version:
//   0.1
// Tested with boost 1.48
// ====================================================

#pragma once

#ifndef SITE_2_H
#define SITE_2_H

#include <string>
#include <vector>

#include "logger.h"
#include "common_2.h"

// ----------------------------------------------------
class Site
{
public:
	Site();
	~Site(){};

	bool is_valid;

	void ReserveSpace(void);

	void Initialize( std::string const message, std::vector<double> const & arr );

	void InitializeDuration( std::string const message, double const p );

	unsigned int width;
	unsigned int order;
	int margin;

	std::vector< std::vector<double> > matrix;
	std::vector<double> non;
	std::vector< std::vector<double> > logodd;
	

	double Get( std::vector<unsigned char> const & nt, unsigned int pos, unsigned int status );
	double GetWithDur( std::vector<unsigned char> const & nt, unsigned int pos, unsigned int status );

	double GetDir( std::vector<unsigned char> const & nt, unsigned int pos );
	double GetRev( std::vector<unsigned char> const & nt, unsigned int pos );

	double GetDirWithDur( std::vector<unsigned char> const & nt, unsigned int pos );
	double GetRevWithDur( std::vector<unsigned char> const & nt, unsigned int pos );
    
    double GetMaxDurationScore() const;

	SiteInfo GetDirWithDurFullInfo( std::vector<unsigned char> const & nt, unsigned int pos );
	SiteInfo GetRevWithDurFullInfo( std::vector<unsigned char> const & nt, unsigned int pos );

	unsigned int max_dur;
	std::vector<double> duration;
	double norm_for_duration;
	std::vector<double> duration_logodd;

	bool verbose;
	unsigned int debug;
	Logger* log_here;

private:

	unsigned int height;

	double unknow_letter_logodd;
	
	void NormalizeArr( std::vector< std::vector<double> > & arr, std::string const message );
	void CheckForZeros(std::vector<double> const & arr);
	void CalcLogOdd(void);
	void ClearSitePositions(void);

	void AbsoluteToConditionalByLastPosition( std::vector< std::vector<double> > & arr );
	void AbsoluteToConditionalByLastPosition( std::vector<double> & arr );
};
// ----------------------------------------------------
#endif // SITE_2_H


//        Site model
//
// height is defined by model order
//
//  order 0 - hight  4
//
//   |  0   A     P(A1)  P(A2)  P(A3)  ...
//   |  1   C     P(C1)  P(C2)  P(C3)  ...
//   |  2   G     P(G1)  P(G2)  P(G3)  ...
//   |  3   T     P(T1)  P(T2)  P(T3)  ...
//
//  order 1 - hight 16
//
//  First column values are from 0-order model
//  This allows order independent calculation of probabiliti
//
//        n n+1    n+1    n+1 n     n+1 n
//   |  0  AA     P(A1)  P(A2|A1)  P(A3|A2)  ...
//   |  1  AC     P(C1)  P(C2|A1)  P(C3|A2)  ...
//   |  2  AG     P(G1)  P(G2|A1)  P(G3|A2)  ...
//   |  3  AT     P(T1)  P(T2|A1)  P(T3|A2)  ...
//
//   |  4  CA     P(A1)  P(A2|C1)  P(A3|C2)  ...
//   |  5  CC     P(C1)  P(C2|C1)  P(C3|C2)  ...
//   |  6  CG     P(G1)  P(G2|C1)  P(G3|C2)  ...
//   |  7  CT     P(T1)  P(T2|C1)  P(T3|C2)  ...
//
//   |  8  GA     P(A1)  P(A2|G1)  P(A3|G2)  ...
//   |  9  GC     P(C1)  P(C2|G1)  P(C3|G2)  ...
//   | 10  GG     P(G1)  P(G2|G1)  P(G3|G2)  ...
//   | 11  GT     P(T1)  P(T2|G1)  P(T3|G2)  ...
//
//   | 12  AT     P(A1)  P(A2|T1)  P(A3|T2)  ...
//   | 13  CT     P(C1)  P(C2|T1)  P(C3|T2)  ...
//   | 14  GT     P(G1)  P(G2|T1)  P(G3|T2)  ...
//   | 15  TT     P(T1)  P(T2|T1)  P(T3|T2)  ...
//
//                -- width --
//


