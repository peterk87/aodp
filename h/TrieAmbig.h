#ifndef __TrieAmbig_h__
#define __TrieAmbig_h__

#include "Trie.h"

using namespace std;

class TrieAmbig: public Trie {
private:
public:
	TrieAmbig( Source& so, Length m, Length M ): Trie( so, m, M ) {};
protected:
	virtual void _cover( bool& first ){
		loop( first, range, lengthRange, elementAdd );
	};

	virtual void _touch( bool& first ){
		loop( first, range, lengthRange, elementMark );
	}

	virtual void buildSlices() {
		for( const auto& e2fr: source.fragments.to ) { // capture all prefixes
 			assert( e2fr.first.getRange().size() > minim );

			const Position first = e2fr.first.getRange().lo();
			const Position last = e2fr.first.getRange().hi() - minim;

			for( Position p = first ; p <= last ; ++p ) {
				newSlice( p );
			}
		}

		cake.resize( prefixes.size());

// 	populate the ambig prefix match
		for( const auto& e2pr1: prefixes ) { // cross-match all prefixes
			prefix_match.emplace( e2pr1.first, e2pr1.second ); // match a prefix with its own slice

			for( const auto& e2pr2: prefixes ) { // unambiguous prefixes match only themselves
				if( p4ma( e2pr1.first, e2pr2.first ) == 4 ) {
					prefix_match.emplace( e2pr1.first, e2pr2.second );
				}
			}
		}

// 	populate the diff1 prefix match
		for( const auto& e2pr1: prefixes ) {
			for( const auto& e2pr2: prefixes ) { // add all diff1 unambiguous prefixes
				if( p4ma( e2pr1.first, e2pr2.first ) == 3 ) {
					prefix_diff1.emplace( e2pr1.first, e2pr2.second );
				}
			}
		}
	}
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
