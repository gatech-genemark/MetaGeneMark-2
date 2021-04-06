// ====================================================
// Author: Alex Lomsadze
// Project: AL library
// File: exit.h
// ====================================================

// ----------------------------------------------------
// exit() vs throw()
// Using exit() as there is no recovery on error in this project
//
// Each Exit() statement must have clear error message on
// "why" and "where in the code" it was called.
// ----------------------------------------------------

#ifndef AL_EXIT_H
#define AL_EXIT_H

#include <cstdlib>

template < class T >
void Exit( T const & t )
{
	std::cerr << t << std::endl;
	exit(1);
}

template < class T1, class T2 >
void Exit( T1 const & t1, T2 const & t2 )
{
	std::cerr << t1 << " " << t2 << std::endl;
	exit(1);
}

template < class T1, class T2, class T3 >
void Exit( T1 const & t1, T2 const & t2, T3 const & t3 )
{
	std::cerr << t1 << " " << t2 << " " << t3 << std::endl;
	exit(1);
}

template < class T1, class T2, class T3, class T4 >
void Exit( T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4 )
{
	std::cerr << t1 << " " << t2 << " " << t3 << " " << t4 << std::endl;
	exit(1);
}
// ----------------------------------------------------
#endif // AL_EXIT_H

