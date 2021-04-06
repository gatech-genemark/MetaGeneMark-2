// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Project: GeneMark.hmm-2 (no introns)
// File: sequence_file_2.h
// ====================================================

#ifndef GMHMMP2_SEQUENCE_FILE_2_H
#define GMHMMP2_SEQUENCE_FILE_2_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "common_2.h"
#include "logger.h"
#include "memory_mapped_file.h"
#include "settings_2.h"

// ----------------------------------------------------
class FastaRecord
{
public:
	FastaRecord()
	{
		ID = 0;
		shape = SHAPE_TYPES::NON;
	}
	
	FastaRecord( unsigned int const id, SHAPE_TYPE const s )
	{
		ID = id;
		shape = s;
	}

	~FastaRecord(){}

	unsigned int ID;	// unique consecutive ID assigned to record in the range 1..N
	std::string  data;	// sequence
	std::string  meta;	// defline from FASTA, LOCUS line from GenBank or ID line from EMBL
	std::string  name;	// name of the sequence from FASTA: ^\s*>\s*(\S+)\s+ , from GenBank: ^LOCUS\s+(\S+)\s+ , from EMBL: ^ID\s+(\S+);\s+
	SHAPE_TYPE   shape;	// "incomplete", "linear" or "circular"
};
// ----------------------------------------------------

typedef std::vector< FastaRecord >            FastaVector;
typedef std::vector< FastaRecord >::iterator  FastaVectorItr;

// ----------------------------------------------------
class SequenceParser
{
public:
	SequenceParser( MemoryMappedFile * const f ): fh(f)
	{
		SetNucAlphabet(nuc_alphabet);
		parse_shape = false;
		default_shape = SHAPE_TYPES::INCOMPLETE;
	}
	virtual ~SequenceParser(){};

	void LoadAll(FastaVector & target);
	bool LoadNext( FastaVector & target );
	bool ValidateFormat(void);
	
	bool parse_shape;
	SHAPE_TYPE default_shape;

protected:
	virtual SHAPE_TYPE GetShape( std::string const & str ) { return default_shape; }
	virtual bool FindHeader( char const * & L, char const * & R, bool exit_on_error ) { R = L; return true; }
	virtual std::string NameFromHeader( char const * L, char const * R ) { return std::string(""); }
	virtual std::string MetaFromHeader( char const * L, char const * R ) { return std::string(""); }
	virtual unsigned int FindSeq( char const * & L, char const * & R ) { return 0; }

	void LoadSequenceData( char const * L, char const * R, unsigned int const size, std::string & target );

	void SetNucAlphabet( std::vector<char> & target );
	std::vector<char> nuc_alphabet;
	MemoryMappedFile * const fh;
};

// ----------------------------------------------------
class FastaParser : public SequenceParser
{
public:
	FastaParser( MemoryMappedFile * const fh ) : SequenceParser(fh) {}
	~FastaParser(){}
private:
	SHAPE_TYPE GetShape( std::string const & str );
	bool FindHeader( char const * & L, char const * & R, bool exit_on_error );
	std::string NameFromHeader( char const * L, char const * R );
	std::string MetaFromHeader( char const * L, char const * R );
	unsigned int FindSeq( char const * & L, char const * & R );
};

// ----------------------------------------------------
class GenbankParser : public SequenceParser
{
public:
	GenbankParser( MemoryMappedFile * const fh ) : SequenceParser(fh) {}
	~GenbankParser(){}
private:
	SHAPE_TYPE GetShape( std::string const & str );
	bool FindHeader( char const * & L, char const * & R, bool exit_on_error );
	std::string NameFromHeader( char const * L, char const * R );
	std::string MetaFromHeader( char const * L, char const * R );
	unsigned int FindSeq( char const * & L, char const * & R );
};

// ----------------------------------------------------
class EMBLParser : public SequenceParser
{
public:
	EMBLParser( MemoryMappedFile * const fh ) : SequenceParser(fh) {}
	~EMBLParser(){}
private:
	SHAPE_TYPE GetShape( std::string const & str );
	bool FindHeader( char const * & L, char const * & R, bool exit_on_error );
	std::string NameFromHeader( char const * L, char const * R );
	std::string MetaFromHeader( char const * L, char const * R );
	unsigned int FindSeq( char const * & L, char const * & R );
};

// ----------------------------------------------------
class PlainParser : public SequenceParser
{
public:
	PlainParser( MemoryMappedFile * const fh ) : SequenceParser(fh) {}
	~PlainParser(){}
private:
	unsigned int FindSeq( char const * & L, char const * & R );
};

// ----------------------------------------------------
class SequenceFile
{
public:
	SequenceFile( Settings const & settings, Logger * const logger );
	~SequenceFile(){ delete parser; }

	FastaVector     data;
	FastaVectorItr  Next( FastaVectorItr itr );
	
	bool AllAtOnce(void) { return !one_by_one; }

	int RecordsIn(void);

	bool IgnoreThis(std::string & s);
	bool RunOnThis(std::string & s);

private:

	// one_by_one == false : all sequences are loaded into "data"
	// one_by_one == true  : sequence is parsed one record at a time and current record is stored in "data[0]" location

	bool one_by_one;

	// file parsing

	MemoryMappedFile  mmf;
	
	SequenceParser * parser;
	SequenceParser * SetFileParser( MemoryMappedFile * const mmfile, FILE_FORMAT & file_format );
	
	void TransferShapesFromFileToData( std::string const & name, SHAPE_TYPE const def_shape );

	std::map<std::string, int> ignore_set;
	std::map<std::string, int> run_on_set;
	std::map<std::string, int> GetSeqIds(std::string const & name);

	std::string Summary(unsigned int const debug);
	std::string RecordSummary(unsigned int const id);

	// settings

	std::string   filename;
	FILE_FORMAT   format;
	SHAPE_TYPE    shape;
	std::string   infofile;
	bool          parse_defline;
	unsigned int  file_load_size;

	Logger * const logger;
};
// ----------------------------------------------------
#endif // GMHMMP2_SEQUENCE_FILE_2_H

