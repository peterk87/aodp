#ifndef __NNParameters_h__
#define __NNParameters_h__

#include <iostream>

#include <array>

#include <cassert>

#include "Types.h"
#include "array2.h"

using namespace std;

/**
 * Collection of nearest-neighbor thermodynamic parameters
 * 
 * Used for storing sequence thermodynamic parameters (initiation, terminal A/T penalty, symmetry)
 * and nearest-neighbor matching values for:
 *  - DS -- entropy  (10x entorpy units - e.u. = cal/mol/K)
 *  - Tp -- pseudo-melting temperature (1000x K)
 * 
 * WARNING: floating point values times a decimal multiplier are stored as signed integers (for faster calculations):
 *  - DS : 10x
 *  - Tp : 1000x
 * 
 * \see references:
 * 
 * SantaLucia J Jr, Hicks D. 2004. The Thermodynamics of DNA Structural Motifs. Annu. Rev. Biophys. 33:415–40
 * Allawi HT, SantaLucia J Jr. 1997. Thermodynamics and NMR of Internal GT Mismatches in DNA. Biochemistry 36:10581-10594
 * Allawi HT, SantaLucia J Jr. 1998a. Nearest Neighbor Thermodynamic Parameters for Internal GA Mismatches in DNA. Biochemistry 37:2170-2179
 * Allawi HT, SantaLucia J Jr. 1998b. Thermodynamics of internal CT mismatches in DNA. Nucleic Acid Res. 26:2694-2701
 * Allawi HT, SantaLucia J Jr. 1998c. Nearest-Neighbor Thermodynamics of Internal A‚C Mismatches in DNA: Sequence Dependence and pH Effects. Biochemistry 37:9435-9444
 * Peyret N, Seneviratne AP, Allawi HT, SantaLucia J Jr. 1999. Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal AA, CC, GG, and TT Mismatches. Biochemistry, 38(12):3468–3477
 * SantaLucia J Jr. 1998 A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc. Natl. Acad. Sci. USA 95:1460–1465
 * Bommarito S, Peyret N, SantaLucia J Jr. 2000. Thermodynamic parameters for DNA sequences with dangling ends. Nucleic Acids Res. 28/9:1929–1934
 * 
 */
class NNParameters {
public:
	NNParameters() : _initiation( 0 ), _terminal_at_penalty( 0 ), _symmetry_correction( 0 ) {
		_nn.fill( 0 );
	};

	/**
	 * Calculate a theromodynamic value for two equal length, human-readable DNA strings,
	 * incrementally, using Thermo values and dimer parameters
	 * 
	 * Use for test driver 
	 */
	int calc( const string& a, const string& b ) const {
		return _calc( convertAsc2Nu( a ), convertAsc2Nu( b ));
	};

	/**
	 * Calculate a theromodynamic value for two equal length, nucleotide encoded DNA strings,
	 * incrementally, using Thermo values and dimer parameters
	 */
	int _calc( const string& a, const string& b ) const {
		int la = a.size();
		int lb = b.size();

		assert( la > 1 );   // minimum length 2
		assert( la == lb );

		int r = 0;

		r += _initiation;

		for( int i = 0 ; i < ( la-1 ) ; i++ ) {
			assert( !nu2ambig.at( a.at( i )));
			assert( !nu2ambig.at( a.at( i+1 )));
			assert( !nu2ambig.at( b.at( i )));
			assert( !nu2ambig.at( b.at( i+1 )));

			r += nn( a.at( i ), a.at( i+1 ), b.at( i ), b.at( i+1 ));
		}

		r += terminal_at_penalty( a, b );
		r += symmetry_correction( a, b );

		return r;
	};

	array<int,256> _nn;
	array<int,256> _terminal_mismatch;
	array<int,512> _loop;  // size of the loop is the sum of the lengths
	array<int,256> _bulge;
	array<int,256> _hairpin;
	array2<int,256,5> _multiloop; // length internal x (helices-3) external

	array<int,64>  _dangx;
	array<int,64>  _dangy;

	unordered_map<string,int> _hairpin_increments;

	int _initiation;
	int _terminal_at_penalty;
	int _symmetry_correction;
	int _hairpin_at_penalty;

	void clear() {
		_initiation = 0;
		_terminal_at_penalty = 0;
		_symmetry_correction = 0;
		_hairpin_at_penalty = 0;

		_nn.fill( 0 );
		_loop.fill( 0 );

		_dangx.fill( 0 );
		_dangy.fill( 0 );
	};
	/**
	 * \param s1 nucleotide encoded string
	 * \param s2 nucleotide encoded string
	 * 
	 * \returns value for dimer pair s1[i]s1[i+1]/s2[j]s2[j+1]
	 * 
	 * WARNING: use in self-duplexes
	 * 
	 *        i
	 *   s1 --**----
	 *        ||
	 *   s2 --**----
	 *        j
	 */
	inline const int nn( const string& s1, int i, const string& s2, int j ) const {
		return nn( s1.at( i ), s1.at( i+1 ), s2.at( j ), s2.at( j+1 ));
	};
	/**
	 * \param s nucleotide encoded string
	 * 
	 * \returns value for dimer pair s[i]s[i+1]/s[j]s[j-1]
	 * 
	 * WARNING: use in self-folding configurations (hairpins)
	 * 
	 *       i
	 *  s  --**----+
	 *       ||    |
	 *     --**----+
	 *       j
	 */
	inline const int nn( const string& s, int i, int j ) const {
		return nn( s.at( i ), s.at( i+1 ), s.at( j ), s.at( j-1 ));
	};
	/**
	 * \returns value for dimer XY/UV in nearest neighbor array
	 */
	inline const int nn( Symbol x, Symbol y, Symbol u, Symbol v ) const {
		assert( !nu2ambig.at( x ));
		assert( !nu2ambig.at( y ));
		assert( !nu2ambig.at( u ));
		assert( !nu2ambig.at( v ));

		int r = _nn.at( nu2pre.at( x ) * 4 + nu2pre.at( y ) + nu2pre.at( u ) * 64 + nu2pre.at( v ) * 16 );

		return r;
	};
	/**
	 * \param s1 nucleotide encoded string
	 * \param s2 nucleotide encoded string
	 * 
	 * \returns value for dimer pair s1[i]s1[i+1]/s2[j]s2[j+1]
	 * 
	 * WARNING: use in self-duplexes
	 * 
	 *        i
	 *   s1 --**----
	 *        ||
	 *   s2 --**----
	 *        j
	 */
	inline const int terminalMismatch( const string& s1, int i, const string& s2, int j ) const {
		return terminalMismatch( s1.at( i ), s1.at( i+1 ), s2.at( j ), s2.at( j+1 ));
	};
	/**
	 * \param s nucleotide encoded string
	 * 
	 * \returns value for dimer pair s[i]s[i+1]/s[j]s[j-1]
	 * 
	 * WARNING: use in self-folding configurations (hairpins)
	 * 
	 *       i
	 *  s  --**----+
	 *       ||    |
	 *     --**----+
	 *       j
	 */
	inline const int terminalMismatch( const string& s, int i, int j ) const {
		return terminalMismatch( s.at( i ), s.at( i+1 ), s.at( j ), s.at( j-1 ));
	};
	/**
	 * \returns value for terminal mismatch XY/UV in nearest neighbor array
	 */
	inline const int terminalMismatch( Symbol x, Symbol y, Symbol u, Symbol v ) const {
		assert( !nu2ambig.at( x ));
		assert( !nu2ambig.at( y ));
		assert( !nu2ambig.at( u ));
		assert( !nu2ambig.at( v ));

		int r = _terminal_mismatch.at( nu2pre.at( x ) * 4 + nu2pre.at( y ) + nu2pre.at( u ) * 64 + nu2pre.at( v ) * 16 );

		return r;
	};

	inline const int loop( int l1, int l2 ) const {
		assert( l1 > 0 ); // not a bulge
		assert( l2 > 0 );

		return _loop.at( l1+l2 ) + abs( l1-l2 ) * 3;
	};
	inline const int bulge( int l ) const {
		return _bulge.at( l );
	};
	/**
	 * 
	 * \param s is the underlying string
	 * \param i index of first character of hairpin (excluding matching characters)
	 * \param j index of last character of hairpin (excluding matching characters)
	 */
	inline const int hairpin( const string& s, int i, int j ) const {
		const int l = j - i + 1;

		assert( l > 2 );

		switch( l ) {
			case 3:
				return hairpin3( s, i, j );
			case 4:
				return hairpin4( s, i, j );
		}

// 		const string h = convertNu2Asc( s.substr( i-1, j-i+2 ));
// 		cout << h << '\t' << _hairpin.at( l ) << endl;

		assert( l > 4 );

		const int h = _hairpin.at( l );
		const int t = terminalMismatch( s, i-1, j+1 );

		return h+t;
	};
	/**
	 * Terminal A/T penalty
	 */
	inline const int terminalATPenalty( const string& s, int i, int j ) const {
		if( !bpAT( s.at( i ), s.at( j ))) return 0;
		return _terminal_at_penalty;
	};
	/**
	 * For hairpins of length 3 (SantaLucia and Hicks 2004; eq. 8):
	 * 
	 * DG37 (total) = DG37 (Hairpin of 3) + DG37 (triloop bonus) + closing AT penalty
	 * 
	 * where DG37 (Hairpin of 3) is +3.5 kcal mol −1 (Table 4) and the closing AT
	 * penalty is +0.5 kcal mol −1 and is applied only to hairpin sequences that are closed by AT
	 * 
	 * WARNING: not using the 0.50 "hairpin closing" AT penalty (as in paper), but the
	 * A/T terminal penalty 0.05 (incl. DH and DS values), as in DINAMelt
	 */
	inline const int hairpin3( const string& s, const int i, const int j ) const {
		const auto hpi = _hairpin_increments.find( s.substr( i-1, j-i+3 ) );
		const int bonus = (( hpi == _hairpin_increments.end()) ? 0 : hpi->second );

		return _hairpin.at( 3 ) + bonus + terminalATPenalty( s, i-1, j+1 );
	};
	/**
	 * For hairpins of length 4 (SantaLucia and Hicks 2004; eq. 9):
	 * DG37 (total) = DG37 (Hairpin of 4) + DG37 (tetraloop bonus) + DG37 (terminal mismatch),
	 * 
	 * where DG37 (Hairpin of 4) is +3.5 kcal mol −1 (Table 4) and DG37 (terminal mismatch)
	 * is the increment for terminal mismatches
	 * 
	 * NOTE: Terminal mismatches are in unpublished results (Varma and SantaLucia).
	 * We ignore them here.
	 * 
	 */
	inline const int hairpin4( const string& s, const int i, const int j ) const {
		const auto hpi = _hairpin_increments.find( s.substr( i-1, j-i+3 ) );

		const int bonus = (( hpi == _hairpin_increments.end()) ? 0 : hpi->second );

		return _hairpin.at( 4 ) + bonus + terminalMismatch( s, i-1, j+1 );
	};
	/**
	 * Energy of a multiloop with \param l free bases and \param x helices
	 * 
	 * NOTE: helices are stored as index+3 in the _multiloop array
	 * NOTE: for helices more than 8, default to 8
	 */
	inline const int multiloop( const int l, const unsigned int x ) const {
		assert( l < 256 );

		if( x < ( _multiloop.size()-3 )) { // helices in table
			return _multiloop.at( l, x-3 );
		}

		return _multiloop.back().at( l );
	};
	/**
	 * \returns value for dangling end XY/V
	 */
	inline const int dangX( const Symbol x, const Symbol y, const Symbol v ) const {
		assert( !nu2ambig.at( x ));
		assert( !nu2ambig.at( y ));
		assert( !nu2ambig.at( v ));

		return _dangx.at( nu2pre.at( x ) * 4 + nu2pre.at( y ) + nu2pre.at( v ) * 16 );
	};
	/**
	 * \returns value for dangling end Y/UV
	 */
	inline const int dangY( Symbol y, Symbol u, Symbol v ) const {
		assert( !nu2ambig.at( y ));
		assert( !nu2ambig.at( u ));
		assert( !nu2ambig.at( v ));

		return _dangy.at( nu2pre.at( y ) + nu2pre.at( u ) * 16 + nu2pre.at( v ) * 4 ) ;
	};
	
	/**
	 * \returns terminal A/T penalty (where applicable) for two nucleotide-encoded strings
	 */
	const int terminal_at_penalty( const string& a, const string& b ) const {
		int r = 0;

		if( bpAT( a.front(), b.front())) r += _terminal_at_penalty;
		if( bpAT( a.back(), b.back()))  r += _terminal_at_penalty;

		return r;
	};
	/**
	 * \returns value for symmetry correction (where applicable) for two nucleotide-encoded strings
	 */
	int symmetry_correction( const string& a, const string& b ) const {
		if( a.compare( -b ) == 0 ) {
			return _symmetry_correction;
		}
		return 0;
	};

	friend ostream& operator<< ( ostream& o, const NNParameters& nn ) { // TEST
		for( int i = 1 ; i < 256 ; i++ ) {
			o << i
				<< '\t' << nn._loop.at( i )
				<< '\t' << nn._bulge.at( i )
				<< '\t' << nn._hairpin.at( i )
			<< endl;
		}

		for( int i = 256 ; i < 512 ; i++ ) {
			o << i
				<< '\t' << nn._loop.at( i )
			<< endl;
		}

		return o;
	};
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
