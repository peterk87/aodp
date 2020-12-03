#ifndef __Fragment_h__
#define __Fragment_h__

#include "Cover.h"

/**
 * A fragment from a Sequence. This can be a gene from a genome, a reverse complement of a Sequence, etc
 */
class Fragment {
public:
	/**
	 * Whether this is a reverse complement fragment
	 * Necessary for encoding the name
	 */
	const bool b_reverse_complement;
	Position start;
	string file_name;

	Cover<Position> ambig;
	Cover<Position> ambig_plus;
	Cover<Position> ambig_compl;
	Cover<Position> cover_range;

	Fragment( const string& fn, const Cover<Position>& a, Length M, const bool rc )
		: b_reverse_complement( rc ),
		start( ambig.range().lo()), file_name( fn ), ambig( a ), ambig_plus( a ), ambig_compl( -a ), cover_range() {
		ambig_plus += ( M-1 );

		cover_range += ambig.range();
		cover_range &= ambig.range();
	};

	const Cover<Position>& getAmbig() const        { return ambig; }
	const Cover<Position>& getAmbigPlus() const    { return ambig_plus; }
	const Cover<Position>& getAmbigCompl() const   { return ambig_compl; }
	const Cover<Position>& getRangeAsCover() const { return cover_range; }

	const Range<Position>& getRange() const 		{ return ambig.range(); }

	/**
	 * Print the reverse complement part of a Signature Oligo id:
	 *   ""    if the SO comes from the direct strand
	 *   "-rc" otherwise
	 */
	inline const string& rcId() const {
		const static string rc_string( "-rc" );
		const static string direct_string( "" );

		return b_reverse_complement ? rc_string : direct_string;
	};

// 	necessary for inclusion in map
	friend inline bool operator< ( const Fragment& c1, const Fragment& c2 ) {
		const Range<Position>& r1 = c1.ambig.range();
		const Range<Position>& r2 = c2.ambig.range();

		if( r1.lo() < r2.lo()) {
			return true;
		}

		if( r1.lo() > r2.lo()) {
			return false;
		}

		return r1.hi() < r2.hi();
	};

// 	TEST
	friend ostream& operator<<( ostream& o, const Fragment& co ) {
		o << "{" << co.file_name << ":" << co.ambig << "}";
		return o;
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
