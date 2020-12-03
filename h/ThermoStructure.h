#ifndef __ThermoStructure_h__
#define __ThermoStructure_h__

#include <deque>

#include "Types.h"
#include "array2.h"

#include "Thermo.h"

/**
 * Recursive local thermodynamic structure for a nucleotide sequence, used for
 * MFE dynamic programming algorithms
 * 
 * NOTE: abstract class; specific classes define specific branch comparison, match and mismatch mechanism
 * 
 * \see Fold, HomoDimer, CoDimer
 */
class ThermoStructure {
public:
	ThermoStructure( const string& _s, const Position _lo, const Position _le, deque<Length>& _o, const int _mi, const int _ma, const Thermo& _th, ostream& _out );

	/**
	 * Virtual destructor is necessary for deleting derived objects with pointer type ThermoStructure*
	 */
	virtual ~ThermoStructure() = default;
	virtual void fold() = 0;

protected:
	const static int max_length = 256; // maximum length in nucleotides of a single strand in the (folded or annealing) structure

	const Thermo& th;       // thermodynamic constants (SantaLucia and Hicks 2004)

	const string& s;        // underlying nucleotide-encoded string
	const Position lo;      // start of portion of string to analyze
	const Position le;      // length of portion of string to analyze

	deque<Length>& o;       // lengths of oligo candidates at given position in underlying nucleotide-encoded string

	const Length mi;           // minimum strand length
	const Length ma;           // maximum strand length

	ostream& out;           // output stream for displaying the secondary structure and melting temperature


	sqm<int,max_length> X;   // optimization criterion: DG (free energy) or TP (pseudo-melting temperature)
	                         // Objective: MINIMIZE this using dynamic programming
	sqm<int,max_length> AX;  // accumulated stack value for the optimization criterion (if the stack would be closed by a match)

	sq<int, max_length> S;   // optimum fold entropy (necessary for calculating the actual melting temperature)
	sq<int, max_length> AS;  // accumulated stack entropy (for a possible stack of depth D)

	sq<int, max_length> D;   // depth of a possible stack ending at (i,j) containing only matches and internal single mismatches
	sq<int, max_length> K;   // dynamic programming traceback
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
