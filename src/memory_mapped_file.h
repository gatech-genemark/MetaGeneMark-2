// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Copyright: GeneProbe Inc.
// Release: 2018
// Tested with boost 1.48 <boost/iostreams/device/mapped_file.hpp>
// File: memory_mapped_file.h
// Project: AL library
// ====================================================

#ifndef AL_MEMORY_MAPPED_FILE_H
#define AL_MEMORY_MAPPED_FILE_H

#include <string>
#include <iostream>
#include <fstream>

#include <boost/iostreams/device/mapped_file.hpp>

#include "exit.h"

// ----------------------------------------------------
class MemoryMappedFile
{
public:
	MemoryMappedFile()
	{
		begin = 0;
		end = 0;
		current = 0;
	}
	~MemoryMappedFile(){}

	void MapFile(std::string const & name)
	{
		if ( name.empty() )
			Exit( "error, file name is empty in 'MemoryMappedFile::MapFile'" );

		try
		{
			if ( mfile.is_open() ) mfile.close();

			// Memory mapping of empty files is treated as error by boost
			// Empty files are allowed, so there should be no error message
			// For empty files set begin=end=0 ; do not map such files

			std::streamoff  file_size = 0;

			// open with pointer at the end of file

			std::ifstream in( name.c_str(), std::ifstream::binary | std::ifstream::ate );

			if ( ! in.is_open() )
				Exit( "error on open file:", name );

			file_size = in.tellg();

			if ( file_size == -1 )
				Exit( "error on read file:", name );

			in.close();

			params.path = name;

			if ( file_size )
			{
				mfile.open( params );

				begin = mfile.data();
				end = begin + mfile.size();
			}
			else
			{
				begin = 0;
				end = 0;
			}

			current = begin;
		}
		catch( std::exception const & e )
		{
			Exit( "error on open file:", params.path, e.what() );
		}
	}

	char const * Begin(void)   { return begin; }
	char const * End(void)     { return end; }
	char const * Current(void) { return current; }

	std::streamoff Size(void) { return mfile.is_open() ? mfile.size() : 0; }
	
	bool EndOfFile(void) { return current < end ? false : true; }

	bool NothingLeft(void)
	{
		while( current < end && isspace(*current) ) 
			++current;

		return ( current < end ) ? false : true; 
	}
	
	void FreeFile(void)
	{
		if ( mfile.is_open() )
			mfile.close();
		params.path.clear();
		begin = 0;
		end = 0;
		current = 0;
	}
	
	void SetCurrent( char const * pos )
	{
		if( pos < begin || pos > end )
			Exit( "error, out of the memory mapped region 'MemoryMappedFile::SetCurrent'" );

		current = pos;
	}

	char const * GetLineEND( char const * L )
	{
		while( L <= end )
		{
			if ( *L == '\n' || *L == '\r' )
				break;
			else
				++L;
		}

		return L;
	}

private:
	boost::iostreams::mapped_file_source  mfile;
	boost::iostreams::mapped_file_params  params;

	char const *  begin;
	char const *  end;
	char const *  current;
};
// ----------------------------------------------------
#endif // AL_MEMORY_MAPPED_FILE_H

