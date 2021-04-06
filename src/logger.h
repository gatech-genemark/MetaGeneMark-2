// ====================================================
// Author: Alex Lomsadze
// Project: AL library
// File: logger.h
// ====================================================

// ----------------------------------------------------
// verbose messages are logged into STDOUT
// debug messages are logged into STDERR or log file
// debug mode turns all verbose messages into the debug ones
// ----------------------------------------------------

#ifndef AL_LOGGER_H
#define AL_LOGGER_H

#include <string>
#include <iostream>
#include <fstream>

#include "exit.h"

// ----------------------------------------------------
class Logger
{
public:
	Logger()
	{
		verbose = false;
		debug = 0;
		verbose_logstream = &std::cout;
		debug_logstream = &std::cerr;
	}
	~Logger(){}

	bool verbose;
	unsigned int debug;

	template < class T >
	void Print(T const & t)
	{
		*debug_logstream << t << std::endl;
	}

	template < class T >
	void Print( unsigned int const level, T const & t )
	{
		if ( verbose && level == 0 )
			*verbose_logstream << t << std::endl;

		if ( debug && level <= debug )
			*debug_logstream   << t << std::endl;
	}

	template < class T1, class T2 >
	void Print( unsigned int const level, T1 const & t1, T2 const & t2 )
	{
		if ( verbose && level == 0 )
			*verbose_logstream << t1 << " " << t2 << std::endl;

		if ( debug && level <= debug )
			*debug_logstream   << t1 << " " << t2 << std::endl;
	}

	template < class T1, class T2, class T3 >
	void Print( unsigned int const level, T1 const & t1, T2 const & t2, T3 const & t3 )
	{
		if ( verbose && level == 0 )
			*verbose_logstream << t1 << " " << t2 << " " << t3 << std::endl;

		if ( debug && level <= debug )
			*debug_logstream   << t1 << " " << t2 << " " << t3 << std::endl;
	}

	template < class T1, class T2, class T3, class T4 >
	void Print(unsigned int const level, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4)
	{
		if (verbose && level == 0)
			*verbose_logstream << t1 << " " << t2 << " " << t3 << " " << t4 << std::endl;

		if (debug && level <= debug)
			*debug_logstream   << t1 << " " << t2 << " " << t3 << " " << t4 << std::endl;
	}

	void Open( std::string const & name )
	{
		filestream.open( name.c_str() );

		if ( !filestream.is_open() )
			Exit( "error on open log file:", name );

		debug_logstream = &filestream;
	}

private:
	std::ofstream filestream;
	std::ostream * verbose_logstream;
	std::ostream * debug_logstream;
};
// ----------------------------------------------------
#endif // AL_LOGGER_H

