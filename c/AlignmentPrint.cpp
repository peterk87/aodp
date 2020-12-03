#define __AlignmentPrint_cpp__

#include "AlignmentPrint.h"

void AlignmentPrint::init( const LLength co, const LLength ro ) {
	Alignment::init( co, ro );

	T.at( 0 ).at( 0 ) = '*';
	for( LLength j = 0 ; j < co ; j++ )  T.at( 0 ).at( j+1 ) = '-';
	for( LLength i = 0 ; i < ro ; i++ )  T.at( i+1 ).at( 0 ) = '|';
}

char AlignmentPrint::match( const LLength i, const LLength j,
		const string& c, const Position pc, const LLength lc,
		const string& r, const Position pr, const LLength lr,
		Element& N, Element& NW, Element& W
){
	return T.at( i+1 ).at( j+1 ) = Alignment::match( i, j, c, pc, lc, r, pr, lr, N, NW, W );
}


void AlignmentPrint::print( const string& s1, Position p1, LLength l1, const string& s2, Position p2, LLength l2 ) {
	assert( s1.size() < max_length_smallest_sequence );
	assert( s2.size() < max_length_smallest_sequence );

	if( s1.size() < s2.size()) _print( s2, p2, l2, s1, p1, l1 );
	else _print( s1, p1, l1, s2, p2, l2 );
}

void AlignmentPrint::_print( const string& c, const Position pc, const LLength lc, const string& r, const Position pr, const LLength lr ) {
	assert( T.at( 0 ).at( 0 ) == '*' );

// 	PRINT
	string ca( "" );
	string ra( "" );
	string ee( "" );

	LLength i = lr;
	LLength j = lc;

	while( true ) {
		switch( T.at( i ).at( j )) {
			case '\\':
				if( c.at( pc+j-1 ) & r.at( pr+i-1 )) ee += '|';
				else ee += ' ';
				
				ca += nu2asc.at( c.at( pc + ( --j )));
				ra += nu2asc.at( r.at( pr + ( --i )));
				break;
			case '-':
				ca += nu2asc.at( c.at( pc + ( --j )));
				ra += '-';
				ee += ' ';
				break;
			case '|':
				ca += '-';
				ra += nu2asc.at( r.at( pr + ( --i )));
				ee += ' ';
				break;
			case '*': goto DONE;
			default:
				assert( false );
		}
	}

DONE:
	reverse( ca.begin(), ca.end());
	reverse( ra.begin(), ra.end());
	reverse( ee.begin(), ee.end());

	cout << ca << endl;
	cout << ee << endl;
	cout << ra << endl;

// 	DEBUG
// 	printTable( cout, T, c, r ); // prints the traceback table (WARNING: this can be quite large)
}
