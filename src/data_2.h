// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: data_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#ifndef GMHMMP2_DATA_2_H
#define GMHMMP2_DATA_2_H

#include <string>
#include <vector>

#include "logger.h"
#include "settings_2.h"
#include "pset_2.h"
#include "common_2.h"

// ----------------------------------------------------
class Data
{
public:

	Data( Settings & settings, Pset & pset, Logger * const logger );
	~Data(){};

	std::vector<unsigned char> nt;
	std::vector<unsigned int>  gc;
	std::vector<unsigned int>  flag;

	std::map< int, INTERVAL_EVI > evi_dir_orf;
	std::map< int, INTERVAL_EVI > evi_rev_orf;

	void Set(std::string const & source);

	unsigned int GCasINT( unsigned int const L, unsigned int const R );
	unsigned int CountFlags(std::vector<unsigned int> const & source);

	void ApplyEvidence( std::map< std::string, std::vector<INTERVAL_EVI> > const & evi, std::string const & name );

private:

	void LoadNT(std::string const & source);
	void CountCumulativeGC(std::vector<unsigned char> const & source);
	void SetFlags(std::vector<unsigned char> const & source);
	void SetFlagsAtGaps(unsigned int len);

	std::string PrintNT(void);
	std::string PrintFlags(void);

	void SetNucDataAlphabet( bool const mask, std::vector<char> & target );	
	void SetCodonToFunctionTable( unsigned int const gcode, STRAND_TYPE const st, std::vector<unsigned int> & target );

	std::vector<char> nuc_data_alphabet;
	std::vector<unsigned int> stop_and_start_codon_alphabet;

	void SetNON( int L, int R );
	void SetMASK( int L, int R );
	void SetNoFullOverlap( int L, int R );
	
	void SetCOD( INTERVAL_EVI const & d );
	void SetCDS( INTERVAL_EVI const & d );

	void MergeGaps(std::vector< INTERVAL> & d, int max );

	// settings

	unsigned int  incomplete_at_gaps;
	STRAND_TYPE strand;
	bool softmask;
	unsigned int gcode;
	unsigned int min_gene_length;

	Logger * const logger;
};
// ----------------------------------------------------
#endif // GMHMMP2_DATA_2_H

