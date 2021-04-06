// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: parameter_parser_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================
// to do: move to 'strtod' and 'strtol' ? or not
// to do: in load arrays no checking for the index duplication
// to do: in AsString no cheking for "\0" inside the string
// ====================================================

#pragma once

#ifndef PARAMETER_PARSER_2_H
#define PARAMETER_PARSER_2_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <utility>

typedef std::map< std::string, std::pair< std::size_t, std::size_t > > parameter_map;

//--------------------------------------------------
class ParameterParser
{
public:

	ParameterParser(){};
	~ParameterParser(){};

	std::string buffer;

	void PutFileToBuffer( std::string const file_name ) { PutFileToBuffer( buffer, file_name ); }
	void PutFileToBuffer( std::string & buffer, std::string const file_name );

	void LoadFromFile     ( parameter_map & target, std::string & buffer, std::string const separator, std::string const file_name );
	void LoadFromString   ( parameter_map & target, std::string & buffer, std::string const separator );
	void LoadFromSubString( parameter_map & target, std::string & buffer, std::string const separator, std::size_t const L, std::size_t const R );

	bool IsKey( parameter_map const & m, std::string const key ) { return ( m.find(key) != m.end() ) ? true : false; };

	std::string  asString  ( parameter_map const & m, std::string const key, std::string const & source );
	int          asInt     ( parameter_map const & m, std::string const key, std::string const & source );
	unsigned int asPInt    ( parameter_map const & m, std::string const key, std::string const & source );
	double       asDouble  ( parameter_map const & m, std::string const key, std::string const & source );
	double       asPDouble ( parameter_map const & m, std::string const key, std::string const & source );
	bool         asBool    ( parameter_map const & m, std::string const key, std::string const & source );

	void asVectorOfDoublesWithLabel ( parameter_map const & m, std::string const key, std::string const & source, std::vector<double> & arr );
	void as3VectorOfDoublesWithLabel( parameter_map const & m, std::string const key, std::string const & source, std::vector<double> & arr1, std::vector<double> & arr2, std::vector<double> & arr3 );
	void asMatrixOfDoublesWithLabel ( parameter_map const & m, std::string const key, std::string const & source, std::vector< std::vector<double> > & arr );
	void asVectorOfDoublesWithPos   ( parameter_map const & m, std::string const key, std::string const & source, std::vector<double> & arr );

	std::string Print( parameter_map const & m, std::string const & source );

private:

	void WhiteoutComments( std::string & str, char const ch );
	void WhiteSpaceToChar( std::string & str, char const ch );
	void Mapenize( parameter_map & target, std::string const & str, std::string const separator, char const space, std::size_t const L, std::size_t const R );
	void StopIfNotSeparator( std::size_t const pos, std::string const & str, std::string const separator, const char space );

	unsigned int StringToInt( std::string const str, unsigned int const size );

	void ErrOnValue( parameter_map const & m, std::string const key, std::string const & source );
	void ErrOnKey(std::string const key);
};
//--------------------------------------------------
#endif // PARAMETER_PARSER_2_H

