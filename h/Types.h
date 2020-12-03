#ifndef __Types_h__
#define __Types_h__

#include <iostream>
#include <iomanip>
#include <bitset>
#include <map>
#include <algorithm>
#include <limits>

#include <cassert>

#include "Error.h"

using namespace std;

// TODO: change types to classes (LOW)

// NOTE: 4-bit LSB portion used only (ACGT) and compound symbols
// 15 symbols total
// WARNING: 0000 is an invalid symbol
typedef unsigned char Symbol;

typedef unsigned int Sequence;
typedef unsigned int TypeFragment;
/**
 * Position in the source array (string)
 * 
 * Supports up to 2^32 = 4,294,967,296 total length of the source
 */
typedef unsigned int Position;
typedef unsigned char Depth;
typedef unsigned char Length;
typedef unsigned int  LLength;

typedef unsigned int Species;
typedef unsigned int Cluster;

const Cluster Cluster_invalid = numeric_limits<Cluster>::max();

/**
 * Index in the array (deque) of Trie Slices
 */
typedef unsigned short Slice;
/**
 * Encoded 4-base subsequence prefix; used to calculate the slice corresponding to a subsequence
 * 
 * WARNING: fixed length
 */
typedef unsigned short Prefix4;

// TODO: generate exception for nodes over limit (LOW)
// WARNING:	Supports up to 2^(32-4) = 2^28 = 268,435,456 nodes
// NOTE: 	There are less nodes than nucleotides in the database
typedef unsigned int Node;
// symbol | ( node << 4 )
typedef unsigned int SymbolNode;

inline SymbolNode symbolNode( Symbol sy, Node no ) {
// 	cout << ":" << co2ba( sy ) << endl;
	return
// 	WARNING: overwrite MSB of nodes over 2^28
		( no << 4 ) |
// 	NOTE: protect against incorrect symbols (over 16)
		( sy & 0xF )
		;
}

inline Symbol symbol( SymbolNode syno ) {
	return Symbol( syno & 0xF );
}

inline Node node( SymbolNode syno ) {
	return syno >> 4;
}

typedef unsigned long long PositionLength;

inline PositionLength positionLength( Position p, Length l ) {
	return
		( PositionLength( p ) << 8 ) |
		l;
}

inline Position position( PositionLength pl ) {
	return Position( pl >> 8 );
}

inline Length length( PositionLength pl ) {
	return Length( pl & 0xFF );
}

// TEST
// inline ostream& operator<<( ostream& out, unsigned char c ) {
// 	out << int( c );
// 	return out;
// }

typedef unsigned long long PositionDepthLength;

inline PositionDepthLength positionDepthLength( Position p, Depth d, Length l ) {
	return
		( PositionDepthLength( p ) << 16 ) |
		( PositionDepthLength( d ) << 8 )  |
		PositionDepthLength( l )
		;
}

inline Position pdlPosition( PositionDepthLength pdl ) {
	return Position( pdl >> 16 );
}

inline Depth pdlDepth( PositionDepthLength pdl ) {
	return Length(( pdl >> 8 ) & 0xFF );
}

inline Length pdlLength( PositionDepthLength pdl ) {
	return Depth( pdl & 0xFF );
}

inline bool pdlCompare( PositionDepthLength pdl1, PositionDepthLength pdl2 ) {
	Position start1 = pdlPosition( pdl1 )-pdlDepth( pdl1 );
	Position start2 = pdlPosition( pdl2 )-pdlDepth( pdl2 );

	if( start1 < start2 ) {
		return true;
	}

	if( start1 > start2 ) {
		return false;
	}

	return (( pdlDepth( pdl1 ) + pdlLength( pdl1 )) < ( pdlDepth( pdl2 ) + pdlLength( pdl2 )));
}

class Distribution : public map<int,int> {};

inline ostream& operator<< ( ostream& out, Distribution& di ) {
	for( auto& el: di ) {
		out << el.first << '\t' << el.second << endl;
	}

	return out;
}

/*
 Conversion table

      asc                            nu              pre
   ----------    ---------      ------------     ----------
   Nucleotide    Expansion       4-bit code      2-bit code
(human-readable)                (fast ambig)    (Trie prefix
                                                  encoding) 
   ----------    ---------      ------------     ----------
      A a            -            0001   0x1        00  0x0
      C c            -            0010   0x2        01  0x1
      G g            -            0100   0x4        10  0x2
      T t            -            1000   0x8        11  0x3

      R r          A or G         0101   0x5           -
      Y y          C or T         1010   0xA           -
      W w          A or T         1001   0x9           -
      S s          G or C         0110   0x6           -
      M m          A or C         0011   0x3           -
      K k          G or T         1100   0xC           -

      H h        A or C or T      1011   0xB           -
      B b        C or G or T      1110   0xE           -
      V v        A or C or G      0111   0X7           -
      D d        A or G or T      1101   0xD           -

      N n      A or C or G or T   1111   0xF           -

  -----------------------------------------------------------+
      \0                         0000   0x0           0xF    |   invalid

 WARNING: array[index] has undefined behaviour on overflow; use array.at(index)
*/

/**
 * Conversion table from human-readable nucleotides to 4-bit representation
 * 
 * WARNING: representation must be ASCII
 */
const array<Symbol,128> asc2nu = {
   0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 
   0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 
   0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 
   0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 

//      *A   *B   *C   *D    E    F   *G   *H    I    J   *K    L   *M   *N    O
   0x0, 0x1, 0xE, 0x2, 0xD, 0x0, 0x0, 0x4, 0xB, 0x0, 0x0, 0xC, 0x0, 0x3, 0xF, 0x0, 

//  P    Q   *R   *S   *T    U   *V    W    X   *Y    Z
   0x0, 0x0, 0x5, 0x6, 0x8, 0x0, 0x7, 0x9, 0x0, 0xA, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 

//      *a   *b   *c   *d    e    f   *g   *h    i    j   *k    l   *m   *n    o
   0x0, 0x1, 0xE, 0x2, 0xD, 0x0, 0x0, 0x4, 0xB, 0x0, 0x0, 0xC, 0x0, 0x3, 0xF, 0x0, 

//  p    q   *r   *s   *t    u   *v    w    x   *y    z
   0x0, 0x0, 0x5, 0x6, 0x8, 0x0, 0x7, 0x9, 0x0, 0xA, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0
};

/**
 * Conversion table: 4-bit nucleotides to ASCII human-readable codes
 */
const array<char,16> nu2asc = {
//  0000  0001  0010  0011  0100  0101  0110  0111  1000  1001  1010  1011  1100  1101  1110  1111
//   0x0   0x1   0x2   0x3   0x4   0x5   0x6   0x7   0x8   0x9   0xA   0xB   0xC   0xD   0xE   0xF
	 '*',  'A',  'C',  'M',  'G',  'R',  'S',  'V',  'T',  'W',  'Y',  'H',  'K',  'D',  'B',  'N'
};

/**
 * Conversion table: 4-bit nucleotides to ambiguous (true/1) or unambiguous (false/0)
 */
const array<bool,16> nu2ambig = {
//           A     C     M     G     R     S     V     T     W     Y     H     K     D     B     N 
//  0000  0001  0010  0011  0100  0101  0110  0111  1000  1001  1010  1011  1100  1101  1110  1111
//   0x0   0x1   0x2   0x3   0x4   0x5   0x6   0x7   0x8   0x9   0xA   0xB   0xC   0xD   0xE   0xF
   false,false,false, true,false, true, true, true,false, true, true, true, true, true, true, true
};

/**
 * Conversion table: 4-bit nucleotides to 2-bit unambiguous nucleotides for Trie prefix representation
 */
const array<Symbol,16> nu2pre = {
//    -    'A'   'C'    -    'G'    -     -     -    'T'    -     -     -     -     -     -     -
//  0000  0001  0010  0011  0100  0101  0110  0111  1000  1001  1010  1011  1100  1101  1110  1111
//   0x0   0x1   0x2   0x3   0x4   0x5   0x6   0x7   0x8   0x9   0xA   0xB   0xC   0xD   0xE   0xF
     0xF,  0x0,  0x1,  0xF,  0x2,  0xF,  0xF,  0xF,  0x3,  0xF,  0xF,  0xF,  0xF,  0xF,  0xF,  0xF
//     -    00    01    -     10    -     -     -     11    -     -     -     -     -     -     -
};

/**
 * Conversion table: 4-bit nucleotide to its 4-bit complement
 * 
 * WARNING: generated table
 */
const array<Symbol,16> nu2compl = {
//    -    'A'   'C'   'M'   'G'   'R'   'S'   'V'   'T'   'W'   'Y'   'H'   'K'   'D'   'B'   'N'
//  0000  0001  0010  0011  0100  0101  0110  0111  1000  1001  1010  1011  1100  1101  1110  1111
//   0x0   0x1   0x2   0x3   0x4   0x5   0x6   0x7   0x8   0x9   0xA   0xB   0xC   0xD   0xE   0xF
     0x0,  0x8,  0x4,  0xC,  0x2,  0xA,  0x6,  0xE,  0x1,  0x9,  0x5,  0xD,  0x3,  0xB,  0x7,  0xF,
//    -    'T'   'G'   'K'   'C'   'Y'   'S'   'B'   'A'   'W'   'R'   'D'   'M'   'H'   'V'   'N'
};

/**
 * Lookup table: 4-bit nucleotide x 4-bit nucleotide -> boolean: whether the nucleotides can be complementary
 * 
 * WARNING: generated table
 */
const array<bool,256> nunu2compl = {
// 	*     A     C     M     G     R     S     V     T     W     Y     H     K     D     B     N    
// 	0000  0001  0010  0011  0100  0101  0110  0111  1000  1001  1010  1011  1100  1101  1110  1111    
	0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    //  *  0000
	0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    //  A  0001
	0,    0,    0,    0,    1,    1,    1,    1,    0,    0,    0,    0,    1,    1,    1,    1,    //  C  0010
	0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  M  0011
	0,    0,    1,    1,    0,    0,    1,    1,    0,    0,    1,    1,    0,    0,    1,    1,    //  G  0100
	0,    0,    1,    1,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  R  0101
	0,    0,    1,    1,    1,    1,    1,    1,    0,    0,    1,    1,    1,    1,    1,    1,    //  S  0110
	0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  V  0111
	0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    //  T  1000
	0,    1,    0,    1,    0,    1,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  W  1001
	0,    1,    0,    1,    1,    1,    1,    1,    0,    1,    0,    1,    1,    1,    1,    1,    //  Y  1010
	0,    1,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  H  1011
	0,    1,    1,    1,    0,    1,    1,    1,    0,    1,    1,    1,    0,    1,    1,    1,    //  K  1100
	0,    1,    1,    1,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  D  1101
	0,    1,    1,    1,    1,    1,    1,    1,    0,    1,    1,    1,    1,    1,    1,    1,    //  B  1110
	0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  N  1111
};

/**
 * \returns true if \param sy1 and \param \sy2 can form a Watson-Crick base pair
 * 
 * NOTE: most optimistic case (maximizes stability)
 */
inline bool bp( const Symbol sy1, const Symbol sy2 ) {
	return nunu2compl.at( sy1 * 16 + sy2 );
};

/**
 * Lookup table: 4-bit nucleotide x 4-bit nucleotide -> boolean: whether the nucleotides can be an A/T match
 * 
 * WARNING: generated table
 */
const array<bool,256> nunu2at = {
// 	*     A     C     M     G     R     S     V     T     W     Y     H     K     D     B     N    
// 	0000  0001  0010  0011  0100  0101  0110  0111  1000  1001  1010  1011  1100  1101  1110  1111    
	0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    //  *  0000
	0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    //  A  0001
	0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    //  C  0010
	0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    //  M  0011
	0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    //  G  0100
	0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    //  R  0101
	0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    //  S  0110
	0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,    //  V  0111
	0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    //  T  1000
	0,    1,    0,    1,    0,    1,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  W  1001
	0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    //  Y  1010
	0,    1,    0,    1,    0,    1,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  H  1011
	0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    //  K  1100
	0,    1,    0,    1,    0,    1,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  D  1101
	0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    0,    1,    //  B  1110
	0,    1,    0,    1,    0,    1,    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,    //  N  1111

};

/**
 * \returns true if \param sy1 and \param \sy2 can form a A/T base pair
 * 
 * NOTE: most optimistic case (maximizes stability)
 */
inline bool bpAT( const Symbol sy1, const Symbol sy2 ) {
	return nunu2at.at( sy1 * 16 + sy2 );
};

/**
 * Conversion table: 2-bit unambiguous nucleotides for Trie prefix representation to 4-bit nucleotides
 */
const array<Symbol,4> pre2nu = {
//   00   01   10   11
//  0x0  0x1  0x2  0x3
    0x1, 0x2, 0x4, 0x8
};

inline string convertAsc2Nu( const string& s ) {
	string r = s;
	int i = 0;
	for( char& c: r ) {
		c = asc2nu.at( c );
		if( !c ) { // BUG
			error( "incorrect nucleotide symbol(", s.at( i ), ")");
		}
		i++;
	}
	return r;
}

/**
 * Convert a 4-bit nucleotide string to its reverse complement
 */
inline string convertNu2ReverseComplement( const string& s ) {
	string rc = s;

	reverse( rc.begin(), rc.end()); // first, reverse in place

	for( char& c: rc ) {
		assert( c > 0 );

		c = nu2compl.at( c );
	}

	return rc;
}

inline string convertNu2Asc( const string& s ) {
	string r = s;
	for( char& c: r ) {
		c = nu2asc.at( c );
	}
	return r;
}

/**
 * Return the maximum self-homologous contiguous suffix of string \param s
 * 
 * \example suffixHomolo( "ACGTTT" ) = pair<Symbol,Length>{ 'T', 3 )
 *                            ***
 */
inline pair<Symbol,Length> suffixHomolo( const string& s ) {
	Length l = 0;
	Symbol c = 0;

	bool first = true;

	for( auto it = s.crbegin() ; it != s.crend() ; it++ ) {
		if( first ) {
			first = false;
			c = *it;
			l = 1;
			continue;
		}

		if( c != *it ) {
			break;
		}

		l++;
	}

	return pair<Symbol,Length>{c,l};
}

/**
 * Return the length of the maximum self-homologous contiguous section of string \param s
 * 
 * \example maxHomolo( "ACCCGTA" ) = 3
 *                       ***
 */
inline Length maxHomolo( const string& s ) {
	Symbol sy0 = 0;
	Length l=1;
	Length maxl=0;
	bool first = true;

	for( Symbol sy: s ) {
		if( first ) {
			first = false;
			sy0 = sy;
			continue;
		}

		if( sy==sy0 ) {
			if( ++l > maxl ) {
				maxl = l;
			}
			continue;
		}

		sy0 = sy;
		l = 1;
	}

	return maxl;
}

/**
 * Create a 4-base prefix from four unambiguous nucleotides
 */
inline Prefix4 nu2p4( const Symbol n1, const Symbol n2, const Symbol n3, const Symbol n4 ) {
	assert( !nu2ambig.at( n1 ));
	assert( !nu2ambig.at( n2 ));
	assert( !nu2ambig.at( n3 ));
	assert( !nu2ambig.at( n4 ));

	Prefix4 r = n1 & 0xF;
	r <<= 4;
	r |= ( n2 & 0xF );
	r <<= 4;
	r |= ( n3 & 0xF );
	r <<= 4;
	r |= ( n4 & 0xF );

	return r;
};

/**
 * Convert a 4-base prefix starting at position \param p in string \param s to an encoded prefix
 */
inline Prefix4 nu2p4( const string& s, Position p ) {
	assert( s.size() >= p+4 );

	Prefix4 r = s.at( p ) & 0xF;
	r <<= 4;
	r |= ( s.at( p+1 ) & 0xF );
	r <<= 4;
	r |= ( s.at( p+2 ) & 0xF );
	r <<= 4;
	r |= ( s.at( p+3 ) & 0xF );

	return r;
};

/**
 * Convert an encoded prefix to a 4-base prefix string
 * 
 * \return reference to shared storage; will be overwritten by subsequent calls (WARNING)
 */
inline const string& p42nu( const Prefix4& pr ) {
 	thread_local string s = "    ";

	s.at( 3 ) = pr     & 0xF;
	s.at( 2 ) = pr>>4  & 0xF;
	s.at( 1 ) = pr>>8  & 0xF;
	s.at( 0 ) = pr>>12 & 0xF;

	return s;
};

/**
 * Number of matches between two encoded prefixes
 */
inline Length p4ma( const Prefix4& pr1, const Prefix4& pr2 ) {
	Prefix4 m = pr1 & pr2;
	return !!( m & 0x000F ) + !!( m & 0x00F0 ) + !!( m & 0x0F00 ) + !!( m & 0xF000 );
};

inline pair<Length,Symbol> p4ho( const Prefix4& pr ) {
	Symbol sy = pr & 0xF;
	if( nu2ambig.at( sy ))     { return pair<Length,Symbol>{ 0, 0 };  } // ambiguous last symbol; no start of a homologous section

	if( sy != ( pr>>4  & 0xF )) { return pair<Length,Symbol>{ 1, sy }; } // one 
	if( sy != ( pr>>8  & 0xF )) { return pair<Length,Symbol>{ 2, sy }; } // two
	if( sy != ( pr>>12 & 0xF )) { return pair<Length,Symbol>{ 3, sy }; } // three
	                              return pair<Length,Symbol>{ 4, sy };   // four

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
