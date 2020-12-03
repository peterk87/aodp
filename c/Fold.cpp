#define __Fold_cpp__

#include "Fold.h"

/**
 * 
 * i: index in string
 * j: ~relative~ index of last element
 * k: ~relative~ index of last element of the first section of a split
 * 
 * 
 *     i              i+k  i+k+1       i+j
 * 
 *   (10)  11  12  13 [14][15] 16  17 [18]
 *     |   |   |   |   |   |   |   |   |
 *     0   1   2   3  (4)  5   6   7  (8)
 *                     k               j
 *           [5]               [4]
 *    |<----l1=k+1 ---->| |<--l2=j-k -->|
 *    (  i, k           ) ( i+k+1,j-k-1 )
 *                  [9]
 *    |<----------- l=j+1 ------------->|
 *    (        i, j                     )
 */
Fold::Fold( const string& _s, const Position _lo, const Position _le, deque<Length>& _o, const int _mi, const int _ma, const Thermo& _th, ostream& _out )
	: ThermoStructure( _s, _lo, _le, _o, _mi, _ma, _th, _out )
{
	for( int i = 0 ; i < max_length ; i++ ) { // initialize hairpin values
		for( int j = 0 ; j < ( min_hairpin+1 ) ; j++ ) {
			O.at( i, j ) = j+1;
			L.at( i, j ) = j+1;
			R.at( i, j ) = j+1;

// 			X.at( i, j ) = 0; // loose single-strand pseudo-melting temperature 0; from free energy = 0 (Rose 2011)
		}
	}
}

bool Fold::foldParameters( int i, int j, double& dg, double& dh, double& ds, double& tm ) {
	if( X.zeroAt( i, j )) return false;

	dg = X.at( i, j ) / 100.0;
	ds = S.at( i, j ) / 10.0;

	dh = th.dgds2dh( X.at( i, j ), S.at( i, j ));
	tm = th.dg2tm( X.at( i, j ), S.at( i, j ), j ) - Thermo::K;

	return true;
}

void Fold::fold(){
	assert( lo+le-ma >= 0 );

	const ios_base::fmtflags floatflag = out.setf( ios::fixed, ios::floatfield ); // remember + set float display to "fixed"
	const streamsize prec = out.precision( 1 ); // remember set precision to 1

	for( Position i = lo+le-min_hairpin-1 ; i > lo+le-ma ; i-- ) { // initialize triangular matrices
		for( Length j = ( min_hairpin+1 ) ; j <= ( ma-1 ) ; j++ ) {
			if( j >= lo+le-i ) break;
			fold( i, j );
		}

		if( !i ) break; // WARNING: special case: avoid skipping over 0U--
	}

	for( Position i = lo+le-ma ; i >= lo ; i-- ) { // main loop: calculate one column at a time, right to left
		for( Length j = ( min_hairpin+1 ) ; j <= ( ma-1 ) ; j++ ) {
			fold( i, j );
		}

		if( !i ) break; // WARNING: special case: for loop test would will never succeed for lo=0 and unsigned Position
	}

	out.precision( prec );
	out.setf( floatflag , ios::floatfield );
}

// for each split k:
//   if( (i,k) is left dangle
//      join( 
//   if( (i+k+1,j-k-1) is right dangle
//
// if( bp( i, j )) // match
// 
void Fold::fold( int i, int j ) {
	for( int k = 0 ; k < j ; k++ ) {
		const int i1 = i;
		const int j1 = k;
		const int i2 = i+k+1;
		const int j2 = j-k-1;

		if( !k ) {
			join( i, k, j );
			continue;
		}

		if( X.copyLess( i, j, X.at( i1, j1 ), X.at( i2, j2 ))) join( i, k, j );
	}

	if( bp( s.at( i ), s.at( i+j ))) {
		match( i, j );
	} else {
		mismatch( i, j );
	}

	if( j >= ( mi-1 )) {   // top element has been computed; show result
		if( o.at( i ) >= j ) { // a smaller length has not already been chosen
			if( th.dgSalt( X.at( i, j ), j ) >= 0 ) {
				o.at( i ) = j; // the maximum length is ~less~ than the length at which the melting temperature is higher than the target
			}

			showFold( i, j );
		}
	}
}

/**
 * \param i position
 * \param j relative index of last element
 * \param k relative index of the last element of the first section of the split
 * 
 * left edge of first strand: i
 * relative index of right edge of first strand : k
 * length of the first strand: k+1
 * 
 * position of second left edge: i+k+1
 * relative index of the right edge of the second strand: j-k-1
 * length of second strand: j-k
 */
void Fold::join( int i, int k, int j ) {
	assert( k >= 0 ); assert( k < j );

	const int i1 = i;
	const int j1 = k;
	const int i2 = i+k+1;
	const int j2 = j-k-1;

	X.set( i, j, X.at( i1, j1 ), X.at( i2, j2 ));

	AX.at( i, j ) = 0;

	S.at( i, j ) = S.at( i1, j1 ) + S.at( i2, j2 );
	AS.at( i, j ) = 0;

	D.at( i, j ) = 0;
	K.at( i, j ) = k;

	O.at( i, j ) = O.at( i1, j1 ) + O.at( i2, j2 );
	H.at( i, j ) = H.at( i1, j1 ) + H.at( i2, j2 );

	if( !H.at( i1, j1 )) { // first section is a single strand
		L.at( i, j ) = L.at( i1, j1 ) + L.at( i2, j2 );

		if( !H.at( i2, j2 )) { // ... and second section is a single strand
			R.at( i, j ) = R.at( i1, j1 ) + R.at( i2, j2 );
		} else {
			R.at( i, j ) = R.at( i2, j2 );
		}

		return;
	}

	if( !H.at( i2, j2 )) { // second section is a single strand
		L.at( i, j ) = L.at( i1, j1 );
		R.at( i, j ) = R.at( i1, j1 ) + j2;

		return;
	}

// 	none of the sections are single strands
	L.at( i, j ) = L.at( i1, j1 );
	R.at( i, j ) = R.at( i2, j2 );
}

/**
 * K[i,j] = 0..j-1 : split
 *        = j : match
 */
void Fold::match( int i, int j ) {
	assert( j >= 2 );

	const int ii = i+1; // coordinates of diagonal
	const int jj = j-2;

	D.at( i, j ) = D.at( ii, jj ) + 1;

	if( !D.at( ii, jj )) { // this is the first match
		loopEnergy( ii, jj, AX.at( i, j ), AS.at( i, j )); // initialize the potential stack energies
// 		all other variables remain unchanged
		return;
	}
		
// here, this is not the first match
	AX.at( i, j ) = AX.at( ii, jj ) + th.DG.nn( s, i, i+j );
	AS.at( i, j ) = AS.at( ii, jj ) + th.DS.nn( s, i, i+j );

	if( X.copyLess( i, j, AX.at( i, j ), th.DG.terminalATPenalty( s, i, j ))) { // this match is the optimum; apply the AT terminal penalty
		S.at( i, j ) = AS.at( i, j ) + th.DS.terminalATPenalty( s, i, j );

		K.at( i, j ) = j;

		L.at( i, j ) = 0;
		R.at( i, j ) = 0;
		O.at( i, j ) = 0;

		H.at( i, j ) = 1; // closing a hairpin creates a helix
	}
}

void Fold::mismatch( int i, int j ) {
	assert( j >= 2 );

	const int ii = i+1; // coordinates of diagonal
	const int jj = j-2;

	if( !D.at( ii, jj ) // cannot start a stack with a mismatch
	  || !bp( s.at( ii ), s.at( ii+jj ))) { // bail after two consecutive mismatches --> internal loop
		D.at( i, j ) = 0;

		AX.at( i, j ) = 0;
		AS.at( i, j ) = 0;

		return; // break the stack
	}

	D.at( i, j ) = D.at( ii, jj ) + 1;

	AX.at( i, j ) = AX.at( ii, jj ) + th.DG.nn( s, i, i+j ); // add single mismatch nn
	AS.at( i, j ) = AS.at( ii, jj ) + th.DS.nn( s, i, i+j ); // add single mismatch nn

	if( D.at( ii, jj ) < 2 ) return; // cannot close a stack with just one match

	if( X.copyLess( i, j, AX.at( ii, jj ), th.DG.terminalMismatch( s, i, i+j ))) { // terminal mismatch is the optimum
		S.at( i, j ) = AS.at( ii, jj ) + th.DS.terminalMismatch( s, i, i+j );

		K.at( i, j ) = j;

		L.at( i, j ) = 1;
		R.at( i, j ) = 1;
		O.at( i, j ) = 2;

		H.at( i, j ) = 1; // closing a hairpin creates a helix

		return;
	}

// 	cout << "#3" << endl;
}

/**
 * Energy if closing the loop (if any) starting at i of length j
 */
void Fold::loopEnergy( int i, int j, int& dg, int& ds ) {
	assert( j>=0 );

	if( j < 2 ) { // cannot close a hairpin of length less than 3
		dg = 0;
		ds = 0;

		return;
	}

	if( !H.at( i, j )) { // no helices so far: HAIRPIN
		dg = th.DG.hairpin( s, i, i+j );
		ds = th.DS.hairpin( s, i, i+j );

		return;
	}

	if( H.at( i, j ) > 1 ) { // two or more helices so far; 3 or more counting the one being closed: MULTILOOP
		dg = th.DG.multiloop( O.at( i, j ), H.at( i, j ) + 1 ); // "+1" comes from counting the current helix
		ds = th.DS.multiloop( O.at( i, j ), H.at( i, j ) + 1 ); // "+1" comes from counting the current helix

		return;
	}

// 	exactly one helix so far
	if( !L.at( i, j ) || !R.at( i, j )) { // at least one of the strands is zero: BULGE
		dg = th.DG.bulge( L.at( i, j )+R.at( i, j ));
		ds = th.DS.bulge( L.at( i, j )+R.at( i, j ));

		return;
	}

// 	exactly one helix so far and non-zero length strands on both sides: INTERNAL LOOP
	dg = th.DG.loop( L.at( i, j ), R.at( i, j ));
	ds = th.DS.loop( L.at( i, j ), R.at( i, j ));
}


// TEST
void Fold::showFold( const int i, const int j ) {
	if( &out == &onull ) return;

	static mutex lock; // serialize printing the folds between multiple threads

	lock.lock();
	double tm = X.zeroAt( i, j ) ? 0.0 : th.dg2tm( X.at( i, j ), S.at( i, j ), j );

	const double zeroC = Thermo::K;
	const double hundredC = 100.0 + Thermo::K;

	out << convertNu2Asc( s.substr( i, j+1 )) << '\t';
	_showFold( i, j );

	out << '\t';

	if( tm > 1 ) out << tm - Thermo::K; // do not show absolute zero

	if(( tm >= hundredC ) || ( tm <= zeroC )) out << "*"; // mark invalid melting temperatures
	out << '\t';

// 	DEBUG
	out
		<< "DG = " << X.at( i, j ) / 100.0 << "  "
		<< "DH = "   << th.dgds2dh( X.at( i, j ), S.at( i, j )) << "  "
		<< "DS = "   << S.at( i, j ) / 10.0  << "  "
		<< "Tm = "   << tm - Thermo::K
	<< endl;

	lock.unlock();
}

// TEST
void Fold::_showFold( const int i, const int j ) {
	const int k = K.at( i, j );

	if( !j ) {
		out << '.';
		return;
	}

	if( k == j ) {
		for( int m = 0 ; m < D.at( i, j ) ; m++ ) {
			if( bp( s.at( i+m ), s.at( i+j-m ))) out << '(';
			else out << '*';
		}
		_showFold( i+D.at( i, j ), j-D.at( i, j )-D.at( i, j ));
		for( int m = D.at( i, j )-1 ; m >= 0  ; m-- ) {
			if( bp( s.at( i+m ), s.at( i+j-m ))) out << ')';
			else out << '*';
		}
		return;
	}

	_showFold( i, k );
	_showFold( i+k+1, j-k-1 );
}

// TEST
void Fold::show( const int i ) {
	out << "=============== X =================" << endl;
	X.put(  out, ma, i, convertNu2Asc( s ), 5 );
	out << "=============== AX=================" << endl;
	AX.put( out, ma, i, convertNu2Asc( s ), 5 );
	out << "=============== S =================" << endl;
	S.put( out, ma, i, convertNu2Asc( s ), 5 );
	out << "=============== AS=================" << endl;
	AS.put( out, ma, i, convertNu2Asc( s ), 5 );
	out << "=============== D =================" << endl;
	D.put(  out, ma, i, convertNu2Asc( s ));
	out << "=============== K =================" << endl;
	K.put(  out, ma, i, convertNu2Asc( s ));

	out << "=============== O =================" << endl;
	O.put(  out, ma, i, convertNu2Asc( s ));
	out << "=============== L =================" << endl;
	L.put(  out, ma, i, convertNu2Asc( s ));
	out << "=============== R =================" << endl;
	R.put(  out, ma, i, convertNu2Asc( s ));
	out << "=============== H =================" << endl;
	H.put(  out, ma, i, convertNu2Asc( s ));
}

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
