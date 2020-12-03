#ifndef __Reference_h__
#define __Reference_h__

#include <cmath>

#include <cassert>

#include "ParserFasta.h"
#include "Trie.h"

class Reference : public ParserFasta {
public:
	Reference( Trie& t, unsigned int thr ) :
		trie( t ),
		threads( thr ),
		sequence_buffer_size( max( thr * 4, 100U ))
	{};

	/**
	 * Event listener for the successful processing of a fragment
	 * \param rc whether this is a "reverse complement" is ignored
	 */
	virtual void onFragment( const string& na, const string& fn, const Cover<Position>& amb, const bool rc = 0 ){
		if( !trie.source.reference.from.count( na )) {
			error( "cannot find sequence with id:", na, "(read in database file", fn, ") in associated taxonomy file" );
		}

		if(( trie.source.max_ambiguities > 0 ) && ( amb.length() > trie.source.max_ambiguities )) {
			trie.source.excluded.emplace_back( na + '\t' + fn );
			return;
		}

		if(( trie.source.max_crowded_ambiguities > 0 ) && ( amb.window( trie.source.maxim ) > trie.source.max_crowded_ambiguities )) {
			trie.source.excluded.emplace_back( na + '\t' + fn );
			return;
		}

		reference_sequences.emplace_back( make_pair( na, amb.range()));

		if( reference_sequences.size() > sequence_buffer_size ) {
			trie.confirm( threads, content, reference_sequences );
			reference_sequences.clear();
			content.clear(); // forget/reuse the storage for the current reference sequence
		}
	};
	inline virtual void onFinish() {
		this->ParserFasta::onFinish();

		trie.confirm( threads, content, reference_sequences );
		reference_sequences.clear();
		content.clear(); // forget/reuse the storage for the current reference sequence
	}
private:
	Trie& trie;
	const unsigned int threads;
	const unsigned int sequence_buffer_size;

	vector<pair<string,Range<Position>>> reference_sequences;
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
