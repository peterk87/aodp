#ifndef __Fold_h__
#define __Fold_h__

#include <array>
#include <deque>

#include <iomanip>
#include <typeinfo>

#include "Types.h"
#include "Thermo.h"

#include "util.h"

#include "ThermoStructure.h"

using namespace std;

/**
 * Support for single-strand self-folding melting temperature calculations
 * 
 * Calculations are based on maximizing the pseudo-melting temperature (\see Thermo::Tp) through dynamic programming
 */
class Fold : public ThermoStructure {
private:
	const static int min_hairpin = 3;

	static_assert( min_hairpin < max_length, "min_hairpin must be smaller than max_length" );

	sq<int, max_length> O;   // total number of outermost unstacked nucleotides
	sq<int, max_length> H;   // number of outermost helices
	/**
	 * L[i,j] = l (length of subsequence) if the subsequence is a "strand"
	 *         = length left dangling on the left side of the fold otherwise
	 * 
	 *    -----------      strand: T[i,j] = l ;  H[i,j] = 0
	 *    0          l
	 *                     dangle: L[i,j] = length on left side of the outermost match (< l; can be 0)
	 *        x..x         dangle: R[i,j] = length on right side of the outermost match (< l; can be 0)
	 *     L /    \ R
	 */
	sq<int, max_length> L;
	sq<int, max_length> R;

public:
	Fold( const string& _s, const Position l, const Position h, deque<Length>& _o, const int _m, const int _M, const Thermo& _th, ostream& fo );

	virtual void fold();

	/**
	 * Retrieve thermodynamic paramaeters of the optimum folded structure at position i/j
	 * 
	 * \returns false if there is no folded structure (DG is default *zero*)
	 */
	bool foldParameters( int i, int j, double& dg, double& dh, double& ds, double& tm );
private:
	void show( const int i );
	inline void showFold( const int i, const int j );

	void fold( int i, int j );
	void join( int i, int k, int j );
	void match( int i, int j );
	void mismatch( int i, int j );
	void loopEnergy( int i, int j, int& t0, int& ds );

	void _showFold( const int i, const int j );
};

#endif

// This file is part of aodp (the Automated Oligonucleotide Design Pipeline)
// 
// (C)	HER MAJESTY THE QUEEN IN RIGHT OF CANADA (2014-2018)
// (C)	Manuel Zahariev mz@alumni.sfu.ca (2000-2008,2014-2018)
// 
// aodp is free software: you can redistribute it and/or
// modify it under the terms of version 3 of the GNU General Public
// License as published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License (version 3) for more details.
// 
// You should have received a copy of the GNU General Public License
// (version 3) along with this program. If not, see
// http://www.gnu.org/licenses/.
