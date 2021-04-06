// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Project: GeneMark.hmm-2 (no introns)
// File: sequence_file_2.cpp
// ====================================================

#include <cctype>
#include <cstring>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using std::string;
using std::vector;
using std::map;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::streamoff;
using std::stringstream;

#include "sequence_file_2.h"
#include "exit.h"

// ----------------------------------------------------
namespace LETTERS
{
	char const ERROR  =  0;
	char const IGNORE = -1;

} // namespace LETTERS

// ====================================================
void SequenceParser::SetNucAlphabet( std::vector<char> & target )
{
	target.assign( 256, LETTERS::ERROR );

	string allow( "acgturyswkmbdhvnACGTURYSWKMBDHVN" );

	for( string::iterator itr = allow.begin(); itr != allow.end(); ++itr )
		target[ *itr ] = *itr;

	string ignore( "0123456789" );
	ignore += ' ';
	ignore += '\t';
	ignore += '\n';
	ignore += '\r';
	ignore += '\f';
	ignore += '\v';

	for( string::iterator itr = ignore.begin(); itr != ignore.end(); ++itr )
		target[ *itr ] = LETTERS::IGNORE;
}
// ----------------------------------------------------
void SequenceParser::LoadSequenceData( char const * L, char const * R, unsigned int const size, std::string & target )
{
	// load 'size' nucleotides starting from 'L' until 'R'
	// target.size == size
	// target.size <= R-L

	target.clear();
	target.reserve( size );

	for ( ; L < R; ++L )
	{
		char ch = nuc_alphabet[ *L ];

		if ( ch == LETTERS::ERROR )
			Exit( "error, unexpected symbol found in sequence:", *L );

		if ( ch != LETTERS::IGNORE )
		{
			target.push_back( ch );

			if ( target.size() == size )
				break;
		}
	}

	if ( target.size() != size )
		Exit( "error, mismatch in expected and observed size of sequence:", target.size(), size );
}
// ----------------------------------------------------
void SequenceParser::LoadAll(FastaVector & target)
{
	unsigned int ID = 0;

	fh->SetCurrent( fh->Begin() );
	
	while( ! fh->NothingLeft() )
	{
		++ID;

		target.push_back(FastaRecord(ID, default_shape));

		// find header

		char const * L = fh->Current();
		char const * R = fh->End();

		FindHeader( L, R, true );

		fh->SetCurrent( R );

		// parse header

		target.back().meta = MetaFromHeader( L, R );
		target.back().name = NameFromHeader( L, R );

		if (parse_shape)
			target.back().shape = GetShape(target.back().meta);

		// find sequence
		
		L = fh->Current();
		R = fh->End();

		unsigned int size = FindSeq( L, R );

		fh->SetCurrent( R );

		// parse sequence

		LoadSequenceData(L, R, size, target.back().data);
	}
}
// ----------------------------------------------------
bool SequenceParser::LoadNext( FastaVector & target )
{
	if ( fh->NothingLeft() )
		return false;
	else
	{
		if ( target.empty() )
			target.push_back( FastaRecord( 1, default_shape ) );
		else
			target.at(0).ID += 1;

		// find header

		char const*	L = fh->Current();
		char const*	R = fh->End();

		FindHeader( L, R, true );

		fh->SetCurrent( R );

		// parse header

		target.at(0).meta = MetaFromHeader( L, R );
		target.at(0).name = NameFromHeader( L, R );
		if (parse_shape)
			target.at(0).shape = GetShape(target.at(0).meta );

		// find sequence

		L = fh->Current();
		R = fh->End();

		unsigned int size = FindSeq( L, R );

		fh->SetCurrent( R );

		// parse sequence

		LoadSequenceData( L, R, size, target.at(0).data );
	}

	return true;
}
// ----------------------------------------------------
bool SequenceParser::ValidateFormat(void)
{
	char const * L = fh->Begin();
	char const * R = fh->End();

	return FindHeader( L, R, false );
}

// ====================================================
bool FastaParser::FindHeader( char const * & L, char const * & R, bool exit_on_error )
{
	bool result = false;
	char const * end = R;

	while ( L < end && isspace( *L ) )
		++L;

	// check for BEGIN_TAG

	if ( end - L >= 1 && *L == '>' )
	{
		// get line

		for( R = L; R < end && *R != '\n' && *R != '\r'; ++R )
			;

		result = true;
	}

	if ( !result && exit_on_error )
		Exit( "error, FASTA defline not found" );

	return result;
}
// ----------------------------------------------------
std::string FastaParser::MetaFromHeader( char const * L, char const * R )
{
	return 	string( L, R - L );
}
// ----------------------------------------------------
unsigned int FastaParser::FindSeq( char const * & L, char const * & R )
{
	unsigned int count = 0;
	char const * end = R;

	while ( L < end && isspace( *L ) )
		++L;

	for ( R = L; R < end && *R  != '>' ; ++R )
	{
		char ch = nuc_alphabet[ *R ];

		if ( ch == LETTERS::ERROR )
			Exit( "error, unexpected letter found in FASTA sequence:", *R );

		if ( ch != LETTERS::IGNORE )
			++count;
	}

	return count;
}
// ----------------------------------------------------
string FastaParser::NameFromHeader( char const * L, char const * R )
{
	char const * end = R;

	// skip '<'
	++L;
	while( L < end && isspace( *L ) )
		++L;

	R = L; 
	while( R < end && !isspace( *R ) )
		++R;
	
	return string( L, R - L );
}
// ----------------------------------------------------
SHAPE_TYPE FastaParser::GetShape( std::string const & str )
{
	SHAPE_TYPE s = default_shape;

	char const LABEL[] = "[topology=";

	string::size_type pos = str.find( LABEL );

	if ( pos != string::npos )
	{
		pos += strlen( LABEL );

		if      ( !str.compare( pos, strlen( TAGS::INCOMPLETE ), TAGS::INCOMPLETE ) )  s = SHAPE_TYPES::INCOMPLETE;
		else if ( !str.compare( pos, strlen( TAGS::CIRCULAR )  , TAGS::CIRCULAR ) )    s = SHAPE_TYPES::CIRCULAR;
		else if ( !str.compare( pos, strlen( TAGS::LINEAR )    , TAGS::LINEAR ) )      s = SHAPE_TYPES::LINEAR;
		else
			Exit( "error, unexpected shape label found:", str );
	}

	return s;
}

// ====================================================
bool GenbankParser::FindHeader( char const * & L, char const * & R, bool exit_on_error )
{
	bool result = false;
	char const * end = R;

	while ( L < end && isspace( *L ) )
		++L;

	// check for BEGIN_TAG

	if ( end - L >= 6 && !strncmp( "LOCUS ", L, 6 ) )
	{
		R = L;
		while( R < end )
		{
			// get line and check for END_TAG

			while( R < end && isspace( *R ) ) 
				++R;

			char const * tag_begin = R;

			// find end-of-line

			while( R < end && *R != '\n' && *R != '\r' )
				++R;

			if ( R - tag_begin >= 6 && !strncmp( "ORIGIN", tag_begin, 6 ) )
			{
				result = true;
				break;
			}
		}
	}

	if ( !result && exit_on_error )
		Exit( "error, GenBank header not found" );

	return result;
}
// ----------------------------------------------------
unsigned int GenbankParser::FindSeq( char const * & L, char const * & R )
{
	unsigned int count = 0;
	char const * end = R;

	while ( L < end && isspace( *L ) )
		++L;

	for ( R = L; R + 1 < end && *R != '/' && *(R+1) != '/'; ++R )
	{
		char ch = nuc_alphabet[ *R ];

		if ( ch == LETTERS::ERROR )
			Exit( "error, unexpected letter found in GenBank sequence:", *R );

		if ( ch != LETTERS::IGNORE )
			++count;
	}

	if ( R + 1 >= end || *R != '/' || *(R+1) != '/' )
		Exit( "error in GenBank format" );

	R += 2;

	return count;
}
// ----------------------------------------------------
string GenbankParser::MetaFromHeader( char const * L, char const * R )
{
	char const * end = R;

	R = L;
	while ( R < end &&  *R != '\n' && *R != '\r' )
		++R;

	return string( L, R - L );
}
// ----------------------------------------------------
string GenbankParser::NameFromHeader( char const * L, char const * R )
{
	char const * end = R;

	// skip "LOCUS "
	L += 6;
	while( L < end && isspace( *L ) )
		++L;

	R = L; 
	while( R < end && !isspace( *R ) )
		++R;

	return 	string( L, R - L );
}
// ----------------------------------------------------
SHAPE_TYPE GenbankParser::GetShape( std::string const & str )
{
	SHAPE_TYPE s = default_shape;

	char const circular[] = " circular ";
	char const linear[]   = " linear ";
	
	string::size_type pos = str.find( circular );
	
	if ( pos != string::npos )
		s = SHAPE_TYPES::CIRCULAR;
	else
	{
		pos = str.find( linear );

		if ( pos != string::npos )
			s = SHAPE_TYPES::LINEAR;
	}
	
	return s;
}

// ====================================================
bool EMBLParser::FindHeader( char const * & L, char const * & R, bool exit_on_error )
{
	bool result = false;
	char const * end = R;

	while ( L < end && isspace( *L ) )
		++L;

	// check for BEGIN_TAG
	
	if ( R - L >= 3 && !strncmp( "ID ", L, 3 ) )
	{
		R = L;
		while( R < end )
		{
			// get line and check for END_TAG

			while( R < end && isspace( *R ) ) 
				++R;

			char const * tag_begin = R;

			// find end-of-line

			while( R < end && *R != '\n' && *R != '\r' )
				++R;

			if ( R - tag_begin >= 3 && !strncmp( "SQ ", tag_begin, 3 ) )
			{
				result = true;
				break;
			}
		}
	}

	if ( !result && exit_on_error )
		Exit( "error, EMBL header not found" );

	return result;
}
// ----------------------------------------------------
unsigned int EMBLParser::FindSeq( char const * & L, char const * & R )
{
	unsigned int count = 0;
	char const * end = R;

	while ( L < end && isspace( *L ) )
		++L;

	for ( R = L; R + 1 < end && *R != '/' && *(R+1) != '/'; ++R )
	{
		char ch = nuc_alphabet[ *R ];

		if ( ch == LETTERS::ERROR )
			Exit( "error, unexpected letter found in EMBL sequence:", *R );

		if ( ch != LETTERS::IGNORE )
			++count;
	}

	if ( R + 1 >= end || *R != '/' || *(R+1) != '/' )
		Exit( "error in EMBL format" );

	R += 2;

	return count;
}
// ----------------------------------------------------
string EMBLParser::MetaFromHeader( char const * L, char const * R )
{
	char const * end = R;

	R = L;
	while ( R < end &&  *R != '\n' && *R != '\r' )
		++R;

	return string( L, R - L );
}
// ----------------------------------------------------
string EMBLParser::NameFromHeader( char const * L, char const * R )
{
	char const * end = R;

	// skip "ID "
	L += 3;
	while( L < end && isspace( *L ) )
		++L;

	R = L; 
	while( R < end && !isspace( *R ) )
		++R;

	return 	string( L, R - L );
}
// ----------------------------------------------------
SHAPE_TYPE EMBLParser::GetShape( std::string const & str )
{
	SHAPE_TYPE s = default_shape;

	char const circular[] = " circular;";
	char const linear[]   = " linear;";
	
	string::size_type pos = str.find( circular );
	
	if ( pos != string::npos )
		s = SHAPE_TYPES::CIRCULAR;
	else
	{
		pos = str.find( linear );

		if ( pos != string::npos )
			s = SHAPE_TYPES::LINEAR;
	}
	
	return s;
}

// ====================================================
unsigned int PlainParser::FindSeq( char const * & L, char const * & R )
{
	unsigned int count = 0;

	char const * end = R;

	while ( L < end && isspace( *L ) )
		++L;

	for ( R = L; R < end; ++R )
	{
		char ch = nuc_alphabet[ *R ];

		if ( ch == LETTERS::ERROR )
			Exit( "error, unexpected letter found in sequence:", *R );

		if ( ch != LETTERS::IGNORE )
			++count;
	}

	return count;
}

// ====================================================
SequenceParser * SequenceFile::SetFileParser( MemoryMappedFile * const mmfile, FILE_FORMAT & file_format )
{
	SequenceParser * parser_ptr = 0;

	switch( file_format )
	{
		case FILE_FORMATS::FASTA :
			parser_ptr = new FastaParser(mmfile);
			break;

		case FILE_FORMATS::GENBANK :
			parser_ptr = new GenbankParser(mmfile);
			break;

		case FILE_FORMATS::EMBL :
			parser_ptr = new EMBLParser(mmfile);
			break;

		case FILE_FORMATS::PLAIN :
			parser_ptr = new PlainParser(mmfile);
			break;

		case FILE_FORMATS::AUTO :

			parser_ptr = new FastaParser(mmfile);

			if ( parser_ptr->ValidateFormat() )
			{
				file_format = FILE_FORMATS::FASTA;
				break;
			}
			else
			{
				delete parser_ptr;
				parser_ptr = 0;
			}

			parser_ptr = new GenbankParser(mmfile);

			if ( parser_ptr->ValidateFormat() )
			{
				file_format = FILE_FORMATS::GENBANK;
				break;
			}
			else
			{
				delete parser_ptr;
				parser_ptr = 0;
			}

			parser_ptr = new EMBLParser(mmfile);

			if ( parser_ptr->ValidateFormat() )
			{
				file_format = FILE_FORMATS::EMBL;
				break;
			}
			else
			{
				delete parser_ptr;
				parser_ptr = 0;
			}
			
			parser_ptr = new PlainParser(mmfile);
			file_format = FILE_FORMATS::PLAIN;
			break;

		default:
			Exit( "error, unsupported file format in 'SequenceFile::SetFileParser'" );
	}

	return parser_ptr;
}
// ----------------------------------------------------
FastaVectorItr SequenceFile::Next( FastaVectorItr itr )
{
	if ( one_by_one )
	{
		// sequnce is loaded one-by-one into the same memory location data[0], so iterator stays at the same location

		if ( parser->LoadNext( data ) )
		{
			// next sequence is now in default position "data[0]"
			// do nothig with iterator

			;
		}
		else
		{
			// end-of-sequence, so move iterator
			// iterator now is equal to ".end()"

			++itr;
		}
	}
	else
	{
		// all the sequences are in the vector
		// standard iterator over vector

		++itr;
	}

	return itr;
}
// ----------------------------------------------------
SequenceFile::SequenceFile( Settings const & settings, Logger * const logger ) : logger(logger)
{
	filename       = settings.in.filename;
	format         = settings.in.format_t;
	shape          = settings.in.shape_t;
	parse_defline  = settings.in.parse_defline;
	infofile       = settings.in.infofile;
	file_load_size = settings.hmm.file_load_size;

	// --

	if( ! settings.in.ignore_records.empty() )
		ignore_set = GetSeqIds(settings.in.ignore_records);

	if (!settings.in.run_on_records.empty())
		run_on_set = GetSeqIds(settings.in.run_on_records);

	mmf.MapFile(filename);

	parser = SetFileParser( &mmf, format );

	parser->default_shape = shape;
	parser->parse_shape = parse_defline;

	if ( mmf.Size() <= file_load_size )
	{
		parser->LoadAll(data);

		one_by_one = false;
	}
	else
	{
		parser->LoadNext( data );

		if ( mmf.NothingLeft() )
			one_by_one = false;
		else
			one_by_one = true;
	}
	
	if ( ! one_by_one )
	{
		if ( !infofile.empty() )
			TransferShapesFromFileToData( infofile, shape );

		mmf.FreeFile();
		delete parser;
		parser = 0;

		if ( logger->verbose ) logger->Print( 0, "# Loading sequence ... done" );
	}
	else
	{
		if ( logger->verbose ) logger->Print( 0, "# Loading sequence one-by-one ..." );
	}

	if ( logger->debug ) logger->Print( Summary(logger->debug) );
}
// ----------------------------------------------------
int SequenceFile::RecordsIn(void)
{
	if (!one_by_one)
		return data.size();
	else
	{
		int counts = 0;
		for (char const *i = mmf.Begin(); i < mmf.End(); ++i)
		{
			if (*i == '>')
				counts += 1;
		}

		return counts;
	}
}
// ----------------------------------------------------
bool SequenceFile::RunOnThis(std::string & s)
{
	bool found = false;

	if (run_on_set.empty())
		return true;

	if (run_on_set.find(s) != run_on_set.end())
		found = true;

	return found;
}
// ----------------------------------------------------
bool SequenceFile::IgnoreThis(std::string & s)
{
	bool found = false;

	if (ignore_set.empty())
		return false;

	if (ignore_set.find(s) != ignore_set.end())
		found = true;

	return found;
}
// ----------------------------------------------------
std::map<std::string, int> SequenceFile::GetSeqIds(std::string const & name)
{
	ifstream file;
	file.open(name.c_str());
	if (!file.is_open())
		Exit("error on open file for reading: ", name);

	string line;
	map<string, int> ignore_set;

	while (getline(file, line))
	{
		std::size_t found = line.find_first_of(" \t");
		if (found != std::string::npos)
		{
			line.erase(found, line.size() - found);
			ignore_set[line] = 1;
		}
	}
	file.close();

	return ignore_set;
}
// ----------------------------------------------------
string SequenceFile::Summary(unsigned int const debug)
{
	stringstream s;

	s << "# START SequenceFile :" << endl;
	s << "file_format : " << format << endl;

	if ( !one_by_one )
	{
		s << "# Load all" << endl;
		s << "data.size() : " << data.size() << endl;

		if (debug > 1)
		{
			for( unsigned int i = 0; i < data.size(); ++i )
				s << RecordSummary(i);
		}

		s << "# END SequenceFile " << endl;
	}
	else
	{
		s << "# Load one-by-one mode" << endl;
		s << "data.size() : " << data.size() << endl;

		if ( debug > 1 && data.size() )
			s << RecordSummary(0);
	}

	return s.str();
}
// ----------------------------------------------------
string SequenceFile::RecordSummary(unsigned int const id)
{
	stringstream s;

	s << "# START SequenceRecord :"               << endl;
	s << "ID : "          << data[id].ID          << endl;
	s << "data.size() : " << data[id].data.size() << endl;
	s << "meta : "        << data[id].meta        << endl;
	s << "name : "        << data[id].name        << endl;
	s << "shape : "       << data[id].shape       << endl;
	s << "# END SequenceRecord :"                 << endl;

	return s.str();
}
// ----------------------------------------------------
void SequenceFile::TransferShapesFromFileToData( std::string const & name, SHAPE_TYPE const def_shape )
{
	ifstream in( name.c_str() );
	if ( !in.is_open() )
		Exit( "error on open file:", name );

	map<string, SHAPE_TYPE> shape_map;
	
	// load seq_name -to- seq_shape information to the map (if any)

	string line;

	while( getline( in, line ) )
	{
		string::size_type pos = line.find_first_not_of( " \t" );
	
		if ( pos == string::npos ) continue;
		if ( line[pos] == '#' )    continue;

		string seq_name = line.substr( pos, line.find_first_of( " \t", pos ) - pos );

		char LABEL[] = "[topology=";

		pos = line.find( LABEL, pos + seq_name.size() );

		SHAPE_TYPE shape = def_shape;

		if ( pos != string::npos )
		{
			pos += strlen( LABEL );

			if      ( !line.compare( pos, strlen( TAGS::INCOMPLETE ), TAGS::INCOMPLETE ) )  shape = SHAPE_TYPES::INCOMPLETE;
			else if ( !line.compare( pos, strlen( TAGS::CIRCULAR )  , TAGS::CIRCULAR ) )    shape = SHAPE_TYPES::CIRCULAR;
			else if ( !line.compare( pos, strlen( TAGS::LINEAR )    , TAGS::LINEAR ) )      shape = SHAPE_TYPES::LINEAR;
			else
				Exit( "error, unsupported shape found in file:", name, line );

			shape_map[ seq_name ] = shape;
		}
	}

	// transfer to vector "data"

	if ( !shape_map.empty() )
	{
		for( FastaVectorItr itr = data.begin(); itr != data.end(); ++itr )
		{
			map<string, SHAPE_TYPE>::iterator i = shape_map.find(itr->name);
			
			if ( i != shape_map.end() )
				itr->shape = i->second;
			else
			{
				if ( logger->verbose )
					logger->Print( 0, "warning, sequence ID from infile not found in sequence file", name, itr->name );
			}
		}
	}

	if ( logger->debug ) logger->Print( 10, "shape_map.size():", shape_map.size() );
}
// ----------------------------------------------------

