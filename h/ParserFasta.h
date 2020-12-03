#ifndef __fasta_h__
#define __fasta_h__

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <cassert>

#include "Cover.h"
#include "Types.h"

#include "Error.h"

using namespace std;

class ParserFasta;

extern int  fasta_parse( ParserFasta&, const string& );
extern void fasta_restart( FILE* );

class ParserFasta {
private:
	const bool b_reverse_complement = false; // whether to generate the reverse complement of sequences in the FASTA file being read
	bool b_fragment = false;

	/**
	 * Indicate that a fragment has been read (event listener)
	 */
	inline void _onFragment() {
		Position le = content.size()-lo;
		Range<Position> ra{ lo, le };

		ambig &= ra;

		onFragment( name, file_name, ambig ); // record the fragment

		if( b_reverse_complement ) { // append the reverse complement
// 	WARNING: The reverse complement fragment is appended immediately after the direct fragment
// 	         Necessary for displaying the location of the reverse complement fragment in the range of the direct fragment
			content += convertNu2ReverseComplement( string( content, lo, le )); // append the reverse complement of the last fragment to the content
			ambig >> le;  // shift the cover of the fragment to the right (to overlap the reverse complement fragment)
			ambig.flip(); // flip the cover (to match the ambiguities in the last fragment)

			onFragment( name, file_name, ambig, true ); // record the reverse complement fragment

// 	WARNING: from here on, ambig represents the cover of the reverse complement, the original cover is long gone
		}

		lo = content.size();

		name.clear();
		file_name.clear();
		ambig.clear();
	};
protected:
	string name;
	string file_name;

	/**
	 * Content storage
	 * 
	 * WARNING: derived classes are responsible with the management of the content storage
	 */
	string content;

	/**
	 * Low boundary for the current fragment
	 */
	Position lo;
	/**
	 * Cover with all ambiguous ranges for the fragment
	 */
	Cover<Position> ambig;

	/**
	 * Event listener for the successful processing of a fragment
	 * 
	 * Derived classes will know what to do here
	 * 
	 * \param na  name of the sequence originating the fragment; multiple fragments with the same name will correspond to the same sequence
	 * \param fn  file name originating the fragment
	 * \param amb cover representing the ambiguous bases of the fragment; NOTE: the range of the cover represents the position of the whole fragment
	 * \param rc  whether the fragment comes from the reverse complement of the sequence
	 */
	virtual void onFragment( const string& na, const string& fn, const Cover<Position>& amb, const bool rc = false ) = 0;
public:
	ParserFasta( const bool rc = false ) : b_reverse_complement( rc ), lo( 0 ) {};

	/**
	 * Encounter a fragment name (event listener)
	 */
	inline void addFragmentName( const string& na, const string& fn ) {
		if( b_fragment ) {
			_onFragment();
		} else {
// 	encountering the first FASTA sequence name means the start of a the first fragment
			b_fragment = true;
		}
		name = na.substr( 0, na.find_first_of( " \t\n\r" ));
		file_name = fn;
		lo = content.size();
	};
	/**
	 * Add a contiguous section of non-ambiguous nucleotides (event listener)
	 */
	inline void addNucleotides( const string& co ) {
		assert( b_fragment );
		content += co;
	};
	/**
	 * Add a contiguous section of ambiguous nucleotides (event listener)
	 */
	inline void addAmbig( const string& co ) {
		assert( b_fragment );

		Position start = content.size();
		Position size = co.size();

		content += co;
		ambig += Range<Position>{ start, size };
	};
	/**
	 * Event listener for the end of processing
	 */
	inline virtual void onFinish() {
		if( b_fragment ) {
// 	if a fragment is being read, trigger its processing
			_onFragment();

			b_fragment = false;
		}
	};

	/**
	 * Interface with the FASTA flex/bison parser
	 * 
	 * The parser will call the event listeners
	 */
	virtual void parse( const string& path ) {
		FILE* fasta_in = fopen( path.c_str(), "rt" );

		if( !fasta_in ) {
			error( "cannot open sequence file ( ", path, " )" );
		}

// 	connect flex to open input file
		fasta_restart( fasta_in );

// 	bison parser call
		if( fasta_parse( *this, path )) {
// 	NOTE: yyparse will return 0 when everything is ok; 1 on YYABORT; 2 on memory exhaustion
// 	reference: http://www.gnu.org/software/bison/manual/html_node/Parser-Function.html
			error( "cannot parse sequence file ( ", path, " )" );
		}

		fclose( fasta_in );
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
