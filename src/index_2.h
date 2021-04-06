// ====================================================
// Author: Alex Lomsadze
// Copyright: Georgia Institute of Technology
// Release: 2017
// File: index_2.h
// Project: GeneMark.hmm-2 (no introns)
// ====================================================

#ifndef GMHMMP2_INDEX_2_H
#define GMHMMP2_INDEX_2_H

//  Index
//
//  order = 0
//  width = 11 binary
//  has_N = 00 or 11 binary
//
//  order = 1
//  width = 1111 binary
//  has_N = 0000 , 0011 , 1100 or 1111 binary

// ----------------------------------------------------
class Index
{
public:

	Index (unsigned int order) : order(order), width((4 << (2*order)) - 1)
	{
		idx   = 0;
		has_N = width;
	}

	~Index(){}

	void Reset(void)
	{
		idx   = 0;
		has_N = width;
	}

	inline unsigned int GetRC(unsigned char i)
	{
		if (i < 4)
			i = 3 - i;

		idx   <<= 2;
		has_N <<= 2;

		idx   &=  width;
		has_N &=  width;

		if (i < 4 )  idx   += i;
		else         has_N += 3;

		return (has_N) ? npos : idx;
	}

	inline unsigned int GetRCwithN(unsigned char i)
	{
		if (i < 4)
			i = 3 - i;

		idx   <<= 2;
		has_N <<= 2;

		idx   &= width;
		has_N &= width;

		if (i < 4 )  idx   += i;
		else         has_N += 3;

		return (has_N) ? (width + 1) : idx;
	}

	inline unsigned int Get(unsigned char i)
	{
		idx   <<= 2;
		has_N <<= 2;

		idx   &= width;
		has_N &= width;

		if (i < 4)  idx   += i;
		else        has_N += 3;

		return (has_N) ? npos : idx;
	}

	inline unsigned int GetWithN(unsigned char const i)
	{
		idx   <<= 2;
		has_N <<= 2;

		idx   &= width;
		has_N &= width;

		if (i < 4)  idx   += i;
		else        has_N += 3;

		return (has_N) ? (width + 1) : idx;
	}

	static const unsigned int npos = -1;

private:

	unsigned int idx;
	unsigned int has_N;

	unsigned int const order;
	unsigned int const width;
};
// ----------------------------------------------------
#endif // GMHMMP2_INDEX_2_H


