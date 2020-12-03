#ifndef __Match_h__
#define __Match_h__

#include <iostream>

#include <unordered_map>
#include <deque>
#include <array>
#include <string>

#include <mutex>

#include "TrieSlice.h"
#include "Trie.h"
#include "util.h"
#include "Alignment.h"

/**
 * Finds cluster matches in a Trie to entries in a FASTA file, then compares the source with the target sequence
 */
class Match : public ParserFasta {
public:
	Match( ostream& o, unsigned int thr, Trie& t, const Length l, bool rc = false ) :
		ParserFasta ( rc ), out( o ), threads( thr ), 
		match_buffer_size( max( thr * 4U, 256U )),
		trie( t ), le( l ) {};
protected:
	virtual void onFragment( const string& na, const string& fn, const Cover<Position>& amb, const bool rc = false );
	virtual void onFinish();
private:
	inline void processFragments();
	inline ostream& print(
		ostream& out,

		const string& target_sequence_name,
		const string& source_sequence_name,

		const double match_percentage,
		const LLength overlap,

		const LLength target_sequence_length,
		const LLength min_set_size,
		const LLength max_set_size = 0
	);
	void sequenceLoop( bool& first );
	void onSequence( const string& na, const Cover<Position>& amb, Alignment& al );

	ostream& out; // output stream
	const unsigned int threads;
	const unsigned int match_buffer_size;

	Trie& trie;
	const Length le;

	deque<pair<string,Cover<Position>>> match_sequences;
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
