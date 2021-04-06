// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// File: output_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS for Visual Studio

#ifndef OUTPUT_2_H
#define OUTPUT_2_H

#include <cstdio>
#include <ctime>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "settings_2.h"
#include "common_2.h"
#include "sequence_file_2.h"
#include "sequence_map_2.h"
#include "logger.h"
#include "pset_2.h"

// ----------------------------------------------------
class Output
{
public:

	Output( Settings const & settings, Logger * const logger, std::string const version );
	~Output(){};

	void Header1(Pset const & set);
	void Header2(FastaVectorItr itr);
	void Footer(void);

	void PrintGenes(std::vector< BestValue > & predictions, FastaVectorItr itr, unsigned int gcode);

	void Stat(FastaVectorItr sitr, double logodd, std::vector< BestValue > & predictions);

	std::map< int, INTERVAL_EVI > * evi_dir_orf;
	std::map< int, INTERVAL_EVI > * evi_rev_orf;

private:

	bool gms2_output;
	bool gms2_training_output;
	bool mgm_output;

	void Header1_LST (void);
	void Header1_GFF (void);
	void Header1_GFF3(void);
	void Header1_shared(Pset const & set);

	void Header2_LST ( FastaVectorItr itr );
	void Header2_GFF3( FastaVectorItr itr );

	void Footer_LST (void);

	unsigned int AssignGeneIDs( std::vector< BestValue > & predictions, unsigned int last_id );

	void PrintGenes_LST ( std::vector< BestValue > & predictions );
	void PrintGenes_GFF ( std::vector< BestValue > & predictions, FastaVectorItr sitr );
	void PrintGenes_GTF ( std::vector< BestValue > & predictions, FastaVectorItr sitr );
	void PrintGenes_GFF3( std::vector< BestValue > & predictions, FastaVectorItr sitr );
	void PrintGenes_EXT (std::vector< BestValue > & predictions, FastaVectorItr sitr);

	void PrintGeneSeqToFile(std::vector< BestValue > & predictions, FastaVectorItr sitr, unsigned int gcode);

	std::string GetNT(std::vector< BestValue >::reverse_iterator itr, std::string & data);
	std::string RevCompNT(std::string & s);
	std::vector<unsigned char> SetComplement(void);
	std::vector<unsigned char> complement;
	std::string Translate(std::string s, unsigned int gcode);
	std::vector<unsigned char> NT2index(void);
	std::vector<unsigned char> index;

	char * ResizeBuffer(char * buf, unsigned int size);
	std::string IntToString(int i);
	std::string DoubleToString(double x, unsigned int i);

	unsigned int gene_count;

	Settings const * st;

	bool gid_per_contig;
	std::string gid_label;
	FILE_FORMAT fp_format;

	std::ofstream out;
	std::ofstream out_nt;
	std::ofstream out_aa;

	std::ofstream * out_nt_ref;
	std::ofstream * out_aa_ref;

	void OpenFile( std::ofstream & ofst, std::string const & name );

	std::string GetTime(void);

	Logger * const logger;
	std::string const version;
};
// ----------------------------------------------------
#endif // OUTPUT_2_H

