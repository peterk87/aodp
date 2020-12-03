#ifndef __Alignment_h__
#define __Alignment_h__

#include <iostream>
#include <iomanip>

#include <cmath>

#include <cassert>

#include "util.h"
#include "Types.h"

using namespace std;

/**
 * Modified Needleman-Wunsch global alignment algorithm
 */
class Alignment {
public:
	/**
	 * Element in dynamic programming single-column table
	 */
	struct Element {
		int S; // actual matching score
		LLength M; // matches
		LLength L; // length of the overlap

		/**
		 * Compares two elements based on their matching scores 
		 */
		bool operator< ( const Element& e ) const {
			return this->S < e.S;
		};
		/**
		 * Element addition
		 */
		Element& operator+= ( const Element& e ) {
			this->S += e.S;
			this->M += e.M;
			this->L += e.L;

			return *this;
		}

		friend ostream& operator<< ( ostream& out, const Element& e ) {
			assert( e.L > 0 );

			out << 100.0 * e.M / e.L << '\t' << e.L << endl;
			return out;
		}
	};

	/**
	 * Maximum length of the smallest sequence
	 */
	static const size_t max_length_smallest_sequence = 4096;

private:
	static const int ma = +1;
	static const int mi = -1;
	static const int ga = -1;

	array<Element,max_length_smallest_sequence> e;
public:
	Element align( const string& s1, const Position p1, const LLength l1, const string& s2, const Position p2, const LLength l2 );

protected:
	Element _align( const string& c, const Position pc, const LLength lc, const string& r, const Position pr, const LLength lr );

	virtual void init( const LLength co, const LLength ro );
	virtual char match( const LLength i, const LLength j,
		const string& c, const Position pc, const LLength lc,
		const string& r, const Position pr, const LLength lr,
		Element& N, Element& NW, Element& W );
};


#endif
