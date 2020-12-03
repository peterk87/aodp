#ifndef __AlignmentPrint_h__
#define __AlignmentPrint_h__

#include "Alignment.h"

/**
 * Modified Needleman-Wunsch global alignment algorithm; prints alignment and traceback table
 */
class AlignmentPrint: public Alignment {
private:
	array<array<char,max_length_smallest_sequence>,max_length_smallest_sequence> T; // traceback path (* - | \)
public:
	void print( const string& s1, Position p1, LLength l1, const string& s2, Position p2, LLength l2 );

protected:
	virtual void init( const LLength co, const LLength ro );
	void _print( const string& c, Position pc, LLength lc, const string& r, Position pr, LLength lr );
	virtual char match( const LLength i, const LLength j,
		const string& c, const Position pc, const LLength lc,
		const string& r, const Position pr, const LLength lr,
		Element& N, Element& NW, Element& W );

template<typename T, size_t l>
	static void printTable( ostream& out, const array<array<T,max_length_smallest_sequence>,max_length_smallest_sequence>& a, const string& c, const string& r ) {
		const LLength co = c.size() + 1;
		const LLength ro = r.size() + 1;

		out << '\t' << '\t';
		for( LLength j = 1 ; j < co ; j++ ) {
			cout << nu2asc.at( c.at( j-1 ));
			out << '\t';
		}
		out << endl;

		for( LLength i = 0 ; i < ro ; i++ ) {
			if( i ) out << nu2asc.at( r.at( i-1 ));

			for( LLength j = 0 ; j < co ; j++ ) {
				out << '\t' << a.at( i ).at( j );
			}

			out << endl;
		}
	}
};

#endif
