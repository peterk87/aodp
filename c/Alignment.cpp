#define __Alignment_cpp__

#include "Alignment.h"

void Alignment::init( const LLength co, const LLength ro ) {
	for( LLength i = 0 ; i < ro ; i++ ) e.at( i ) = { int( i+1 ) * ga, 0 , 0 };
}

/**
 * Align the strings s1 and s2 shifting parameters so that longest string comes first (columns)
 * 
 * \see Align::_align
 */
Alignment::Element Alignment::align( const string& s1, const Position p1, const LLength l1, const string& s2, const Position p2, const LLength l2 ) {
	if( l1 < l2 ) return _align( s2, p2, l2, s1, p1, l1 );
	return _align( s1, p1, l1, s2, p2, l2 );
}

/**
 * Modified Needleman-Wunsch global alignment algorithm
 * 
 * \param c string in columns of the dynamic programming algorithm (longer string)
 * \param r string in rows of the dynamic programming algorithm
 * 
 * Uses standard match/mismatch/gap metrics (+1/-1/-1); maximize "edit" score (S)
 * 
 * Modifications:
 *  - gap penalties are ZERO on the first and last row (skipping the edges of the longest string) WARNING
 *  - maintain number of matches (M) and overlap length (L)
 * 
 * \returns the final alignment score triplet:
 *  - S: actual dynamic programming score (maximize this)
 *  - M: number of matches in the optimal configuration
 *  - L: length of the overlap region between the two strings including gaps
 * 
 *  --AGAC-TAGTTAC
 *    |||| |      <- M=5     S=2 = -1*2    +   1*5   +  -1*1
 *  CGAGACGT------               edge gaps   matches    gaps
 *    ------
 *     L=6         M/L = 5 / 6 = 83.3%
 * 
 * \see Alignment::Element
 */
Alignment::Element Alignment::_align( const string& c, Position pc, LLength lc, const string& r, Position pr, LLength lr ) {
	assert( lc >= lr ); // more columns than rows!

	if( lr >= max_length_smallest_sequence ) error( "both strings to align are longer than maximum length (", max_length_smallest_sequence, ")" );

	init( lc, lr );

	for( LLength j = 0 ; j < lc ; j++ ) { // foreach column
		Element NW = {0,0,0};
		Element N  = {0,0,0};

		for( LLength i = 0 ; i < lr ; i++ ) { // foreach row
			Element W = e.at( i );
			match( i, j, c, pc, lc, r, pr, lr, N, NW, W );
		}
	}

	assert( e.at( lr-1 ).L > 0 ); // overlap length will be used to calculate the overlap percentage

	return e.at( lr-1 );
}

char Alignment::match( const LLength i, const LLength j,
		const string& c, const Position pc, const LLength lc,
		const string& r, const Position pr, const LLength lr,
		Element& N, Element& NW, Element& W
) {
	char d = '-';

	if( i < lr - 1 ) e.at( i ) += { ga, 0, 1 }; // add gap penalty, except on last row

// NW - match or mismatch
	if( r.at( pr+i ) & c.at( pc+j )) NW += { ma, 1, 1 }; // add match bonus
	else NW += { mi, 0, 1 }; // add mismatch penalty

	if( e.at( i ) < NW ) {
		e.at( i ) = NW;
		d = '\\';
	}

// N - vertical gap
	if( j < lc - 1 ) N += { ga, 0, 1 }; // do not count edge gaps
	else N += { ga, 0, 0 };

	if( e.at( i ) < N ) {
		e.at( i ) = N;
		d = '|';
	}

	N  = e.at( i );
	NW = W;

	return d;
}
