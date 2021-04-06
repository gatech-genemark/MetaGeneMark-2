// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: parameter_parser_2.cpp
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#include <cctype>
#include <cstring>

#include <string>
#include <vector>
#include <iostream>
#include <utility>

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::cout;
using std::size_t;
using std::pair;

#include "parameter_parser_2.h"
#include "memory_mapped_file.h"
#include "exit.h"

// ----------------------------------------------------
void ParameterParser::PutFileToBuffer( std::string & buffer, std::string const file_name )
{
	MemoryMappedFile mmf;
	mmf.MapFile(file_name);
	buffer.assign( mmf.Begin(), mmf.End() );
	mmf.FreeFile();

	WhiteoutComments( buffer, 0 );
	WhiteSpaceToChar( buffer, 0 );
}
// ----------------------------------------------------
void ParameterParser::LoadFromFile( parameter_map & target, std::string & buffer, std::string const separator, std::string const file_name )
{
	PutFileToBuffer( buffer, file_name );
	Mapenize( target, buffer, separator, 0, 0, buffer.size() );
}
// ----------------------------------------------------
void ParameterParser::LoadFromString( parameter_map & target, std::string & buffer, std::string const separator )
{
	Mapenize( target, buffer, separator, 0, 0, buffer.size() );
}
// ----------------------------------------------------
void ParameterParser::LoadFromSubString( parameter_map & target, std::string & buffer, std::string const separator, std::size_t const L, std::size_t const R )
{
	Mapenize( target, buffer, separator, 0, L, R );
}
// ----------------------------------------------------
void ParameterParser::WhiteoutComments( std::string & str, char const ch )
{
	size_t commentStart = string::npos;
	size_t commentEnd   = string::npos;

	for( commentStart = str.find( '#', 0 ); commentStart != string::npos; commentStart = str.find( '#', commentEnd ) )
	{
		commentEnd = str.find_first_of( "\n\r", commentStart );
		fill( str.begin() + commentStart, str.begin() + commentEnd, ch );
	}
}
// ----------------------------------------------------
void ParameterParser::WhiteSpaceToChar( std::string & str, char const ch )
{
	for( unsigned int i = 0; i < str.length(); ++i )
	{
		if ( str[i] == '"' )
		{
			unsigned int j = i;

			for( ++j; j < str.length(); ++j )
			{
				if ( str[j] == '"' )
					break;
			}

			i = j;
		}
		else
		{
			if (isspace( str[i] ))
				str[i] = ch;
		}
	}
}
// ----------------------------------------------------
void ParameterParser::StopIfNotSeparator( std::size_t const pos, std::string const & str, std::string const separator, const char space )
{
	bool found = true;

	if ( pos == string::npos )
		found = false;
	else 
	{
		if ( str.compare( pos, separator.size(), separator ) ) 
			found = false;
		else
		{
			// if not the first in string, then space should be before

			if (( pos != 0 )&&( str[pos - 1] != space )) 
				found = false;
	
			// key should follow the separator, not a space

			if (( pos + separator.size() >= str.length() )||( str[ pos + separator.size() ] == space )) 
				found = false;
		}
	}

	if (!found)
		Exit( "error in file format: separator was expected at position: ", separator, pos );
};
// ----------------------------------------------------
void ParameterParser::Mapenize( parameter_map & target, std::string const & str, std::string const separator, char const space, std::size_t const L, std::size_t const R )
{
	assert( L <= str.size() );
	assert( R <= str.size() );
	assert( L <= R );

	size_t keyStart;
	size_t keyEnd;
	size_t valueStart;
	size_t valueEnd;

	size_t current;
	string key;

	// skip space

	current = str.find_first_not_of( space, L );

	while( current != string::npos && current < R )
	{
		StopIfNotSeparator( current, str, separator, space );

		keyStart = current;
		keyEnd = str.find_first_of( space, keyStart );

		if ( keyEnd > R )
			keyEnd = R;

		key = str.substr( keyStart, keyEnd - keyStart );

		// search for next separator

		current = str.find( separator, keyEnd );

		if ( current == string::npos )
			valueEnd = str.size();
		else if ( current > R )
			valueEnd = R;
		else
			valueEnd = current - 1;

		// skip space between key and value
		
		valueStart = str.find_first_not_of( space, keyEnd );

		if ( valueStart > valueEnd )
			valueStart = valueEnd;

		// skip space from end

		valueEnd = str.find_last_not_of( space, valueEnd ) + 1;

		if ( valueEnd < valueStart)
			valueEnd = valueStart;

		if (( valueEnd - valueStart > 1 )&&( str[valueStart] == '"' )&&( str[valueEnd-1] == '"' ))
		{
			++valueStart;
			--valueEnd;
		}

		target[ key ] = std::make_pair( valueStart, valueEnd );
	}
}
// ----------------------------------------------------
void ParameterParser::ErrOnValue( parameter_map const & m, std::string const key, std::string const & source )
{
	Exit( "error in value of the parameter with key:", key, asString( m, key, source ) );
}
// ----------------------------------------------------
void ParameterParser::ErrOnKey(std::string const key)
{
	Exit( "error, parameter not found:", key );
}
// ----------------------------------------------------
string ParameterParser::Print( parameter_map const & m, std::string const & source )
{
	string str;
	for( parameter_map::const_iterator itr = m.begin(); itr != m.end(); ++itr )
	{
		str += itr->first;
		str += ' ';
		str += asString( m, itr->first, source );
		str += '\n';
	}
	return str;
}
// ----------------------------------------------------
std::string ParameterParser::asString( parameter_map const & m, std::string const key, std::string const & source )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);
	
	return source.substr( itr->second.first, itr->second.second - itr->second.first );
}
// ----------------------------------------------------
int ParameterParser::asInt( parameter_map const & m, std::string const key, std::string const & source )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	return atoi(source.c_str() + itr->second.first);
}
// ----------------------------------------------------
unsigned int ParameterParser::asPInt( parameter_map const & m, std::string const key, std::string const & source )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	int i = atoi(source.c_str() + itr->second.first);

	if ( i < 0 ) ErrOnValue( m, key, source );

	return i;
}
// ----------------------------------------------------
double ParameterParser::asDouble( parameter_map const & m, std::string const key, std::string const & source )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	return atof(source.c_str() + itr->second.first);
}
// ----------------------------------------------------
double ParameterParser::asPDouble( parameter_map const & m, std::string const key, std::string const & source )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	double x = atof(source.c_str() + itr->second.first);

	if ( x < 0 ) ErrOnValue( m, key, source );

	return x;
}
// ----------------------------------------------------
bool ParameterParser::asBool( parameter_map const & m, std::string const key, std::string const & source )
{
	return ( asDouble(m,key,source) != 0 );
}
// ----------------------------------------------------
void ParameterParser::asVectorOfDoublesWithLabel( parameter_map const & m, std::string const key, std::string const & source, std::vector<double> & arr )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	// search in this range (L..R)

	size_t L = itr->second.first;
	size_t R = itr->second.second;

	unsigned int size = arr.size();
	
	// count the number of pairs "index-values" found

	unsigned int count = 0;

	char sep = 0;

	unsigned int index;

	// just in case

	size_t current = source.find_first_not_of( sep, L );

	do
	{
		// find index

		if ( current == string::npos || current > R )
			Exit( "error, index was expected for key:", key );

		size_t idxStart = current;
		size_t idxEnd   = source.find( sep, idxStart );
		if ( idxEnd > R ) idxEnd = R;

		index = StringToInt( source.substr( idxStart, idxEnd - idxStart ), size );

		// find value

		current = source.find_first_not_of( sep, idxEnd );

		if ( current == string::npos || current > R )
			Exit( "error, no value for index:", index, key );

		size_t valueStart = current;
		size_t valueEnd   = source.find( sep, valueStart );
		if ( valueStart > R ) valueStart = R;

		arr[index] = atof(source.c_str() + valueStart);

		// all found

		current = source.find_first_not_of( sep, valueEnd );

		++count;
	}
	while( current != string::npos && current < R );

	if ( count != size )
		Exit( "error, diff between expected and observed number of values for key " + key, size, count );
}
// ----------------------------------------------------
void ParameterParser::as3VectorOfDoublesWithLabel( parameter_map const & m, std::string const key, std::string const & source, std::vector<double> & arr1, std::vector<double> & arr2, std::vector<double> & arr3 )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	size_t L = itr->second.first;
	size_t R = itr->second.second;

	unsigned int size = arr1.size();
	unsigned int count = 0;

	char sep = 0;

	unsigned int index;

	size_t current = source.find_first_not_of( sep, L );

	do
	{
		// find index

		if ( current == string::npos || current > R )
			Exit( "error, index was expected for key:", key );

		size_t idxStart = current;
		size_t idxEnd = source.find( sep, idxStart );
		if ( idxEnd > R ) idxEnd = R;
		
		index = StringToInt( source.substr( idxStart, idxEnd - idxStart ), size );

		// value1

		current = source.find_first_not_of( sep, idxEnd );

		if ( current == string::npos || current > R )
			Exit( "error, no value for index:", index, key );

		size_t valueStart = current;
		size_t valueEnd   = source.find( sep, valueStart );
		if ( valueStart > R ) valueStart = R;
		
		arr1[index] = atof(source.c_str() + valueStart);

		// value2
		
		current = source.find_first_not_of( sep, valueEnd );
		
		if ( current == string::npos || current > R )
			Exit( "error, no value for index:", index, key );

		valueStart = current;
		valueEnd = source.find( sep, valueStart );
		if ( valueStart > R ) valueStart = R;
		
		arr2[index] = atof(source.c_str() + valueStart);

		// value3
		
		current = source.find_first_not_of( sep, valueEnd );
		
		if ( current == string::npos || current > R )
			Exit( "error, no value for index:", index, key );

		valueStart = current;
		valueEnd = source.find( sep, valueStart );
		if ( valueStart > R ) valueStart = R;
		
		arr3[index] = atof(source.c_str() + valueStart);

		// all found

		++count;

		current = source.find_first_not_of( sep, valueEnd );
	}
	while( current != string::npos && current < R );

	if ( count != size )
		Exit( "error, diff between expected and observed number of values for key " + key, size, count );
}
// ----------------------------------------------------
unsigned int ParameterParser::StringToInt( string const str, unsigned int const size )
{
	unsigned int size_from_str = 4<<(2*(str.size() - 1));

	if ( size_from_str != size )
		Exit( "error, label has unexpected length:", size, str);

	unsigned int i = 0;

	for( string::const_iterator itr = str.begin(); itr != str.end(); ++itr )
	{
		i <<= 2;
		switch( *itr )
		{
			case 'A':
				i += 0;
				break;
			case 'C':
				i += 1;
				break;
			case 'G':
				i += 2;
				break;
			case 'T':
				i += 3;
				break;
			default:
				Exit( "error, unexpected letter found in label: ", *itr, str );
		}
	}

	return i;
}
// ----------------------------------------------------
void ParameterParser::asMatrixOfDoublesWithLabel( parameter_map const & m, std::string const key, std::string const & source, std::vector< std::vector<double> > & arr )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	// search in this range (L..R)

	size_t L = itr->second.first;
	size_t R = itr->second.second;

	unsigned int hight = arr.size();
	unsigned int width = arr[0].size();

	unsigned int count = 0;

	char sep = 0;
	unsigned int index;

	size_t current = source.find_first_not_of( sep, L );

	do
	{
		// find index

		if ( current == string::npos || current > R )
			Exit( "error, index was expected for key:", key );

		size_t idxStart = current;
		size_t idxEnd = source.find( sep, idxStart );
		if ( idxEnd > R ) idxEnd = R;
		
		index = StringToInt( source.substr( idxStart, idxEnd - idxStart ), hight );

		// values
		current = source.find_first_not_of( sep, idxEnd );

		for( unsigned int i = 0; i < width; ++i )
		{	
			if ( current == string::npos || current > R )
				Exit( "error, no value for index:", index, key );

			size_t valueStart = current;
			size_t valueEnd = source.find( sep, valueStart );
			if ( valueStart > R ) valueStart = R;

			arr[index][i] = atof(source.c_str() + valueStart);

			current = source.find_first_not_of( sep, valueEnd );
		}

		// all found

		++count;

	}
	while( current != string::npos && current < R );

	if ( count != hight )
		Exit( "error, diff between expected and observed number of values for key " + key, hight, count );
}
// ----------------------------------------------------
void ParameterParser::asVectorOfDoublesWithPos( parameter_map const & m, std::string const key, std::string const & source, std::vector<double> & arr )
{
	parameter_map::const_iterator itr = m.find(key);
	if( itr == m.end()) ErrOnKey(key);

	// search in this range (L..R)

	size_t L = itr->second.first;
	size_t R = itr->second.second;

	unsigned int size = arr.size();
	unsigned int count = 0;
	char sep = 0;
	unsigned int index;

	size_t current = source.find_first_not_of( sep, L );

	do
	{
		// find index

		if ( current == string::npos || current > R )
			Exit( "error, index was expected for key:", key );

		size_t idxStart = current;
		size_t idxEnd = source.find( sep, idxStart );
		if ( idxEnd > R ) idxEnd = R;
		
		index = atoi(source.c_str() + idxStart);

		if ( index >= size )
			Exit ( "error, index is out of the range for key:", key, index );

		// values

		current = source.find_first_not_of( sep, idxEnd );

		if ( current == string::npos || current > R )
			Exit( "error, no value for index:", index, key );

		size_t valueStart = current;
		size_t valueEnd = source.find( sep, valueStart );
		if ( valueStart > R ) valueStart = R;

		arr[index] = atof(source.c_str() + valueStart);

		current = source.find_first_not_of( sep, valueEnd );

		// all found

		++count;
	}
	while( current != string::npos && current < R );

	if ( count != size )
		Exit( "error, diff between expected and observed number of values for key " + key, size, count );
}
// ----------------------------------------------------

