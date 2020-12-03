#ifndef __Trie_h__
#define __Trie_h__

#include <iostream>

#include <string>
#include <set>
#include <map>
#include <deque>
#include <unordered_map>
#include <vector>
#include <algorithm>

// TEST
#include <bitset>

#include <thread>
#include <mutex>
#include <atomic>

#include <cassert>

#include "util.h"

#include "Cover.h"
#include "Source.h"

#include "Types.h"

using namespace std;

class TrieSlice;

template<typename T> ostream& operator<< ( ostream& o, const vector<T>& v ) {
	for( const T& e : v ) {
		o << e << ' ';
	}

	return o;
};

class Trie
{
//=============================================================================================
// 	CONSTANTS
//=============================================================================================
public:
	/**
	 * Depth to which the Trie is a lookup table of TrieSlices
	 * 
	 * WARNING: must be 4, since Prefix4 is hardcoded
	 */
	static const Depth fixed_depth  = 4;

//=============================================================================================
// 	DATA members
//=============================================================================================
public:
	Source& source;      // sequences the Trie is read from
	const Length minim;  // minimum length of a signature oligo
	const Length maxim;  // maximum length of a signature oligo

	/**
	 * All collected "group" oligonucleotides
	 * 
	 * Group = set of keys to the sequence in the Source (database)
	 * 
	 * NOTE: an olgonucleotide for a sequence is a oligonucleotide for a group
	 *       with only one element: the sequence itself
	 * NOTE: oligonucleotides for a group are stored in the order they are collected,
	 *       sorting will put them in the order of occurrence in the database
	 * 
	 */
	map<Cluster,vector<PositionDepthLength>> matches;

protected:
	/**
	 * Convert from an encoded prefix to a "cake" sequential index
	 */
	unordered_map<Prefix4,Slice> prefixes;
	/**
	 * List of "cake" indexes that match a prefix (including ambiguous)
	 */
	unordered_multimap<Prefix4,Slice> prefix_match;
	/**
	 * List of "cake" indexes that match a prefix with exactly one difference
	 */
	unordered_multimap<Prefix4,Slice> prefix_diff1;

	/**
	 * Lookup table of TrieSlice-s
	 * 
	 * Allows for multithreading of operations
	 */
	deque<TrieSlice> cake;

//=============================================================================================
// 	FUNCTIONALS; needed for using the same Trie::loop function for iterating over the Trie
//=============================================================================================
protected:

	/**
	 * Functional describing the cover (set of ranges) to iterate over in a generic loop of the Trie
	 * 
	 * \see Trie::loop
	 */
	struct CoverFunction {
		virtual const Cover<Position>& operator() ( const Trie& t, const Fragment& fr ) const = 0;
	};

	const struct : CoverFunction {
		virtual const Cover<Position>& operator() ( const Trie& t, const Fragment& fr ) const {
			return fr.getAmbigPlus();
		};
	} ambigCover;

	const struct : CoverFunction {
		virtual const Cover<Position>& operator() ( const Trie& t, const Fragment& fr ) const {
			return fr.getAmbigCompl();
		};
	} ambigCoverComplement;

	const struct : CoverFunction {
		virtual const Cover<Position>& operator() ( const Trie& t, const Fragment& fr ) const {
			return fr.getRangeAsCover();
		};
	} range;

	/**
	 * Functional describing the function to call for each element in the generic loop of the Trie
	 */
	struct ElementFunction {
		virtual void operator() ( Trie& t, Position po, Length le, Sequence se ) const = 0;
	};

// 	WARNING: some loops use Cover::cover, some Range::cover
	const struct : ElementFunction {
		virtual void operator() ( Trie& t, Position po, Length le, Sequence se ) const {
			t.add( po, le, se );
		};
	} elementAdd;
	
	const struct : ElementFunction {
		virtual void operator() ( Trie& t, Position po, Length le, Sequence se ) const {
			t.mark( po, le, se );
		};
	} elementMark;
	
	const struct : ElementFunction {
		virtual void operator() ( Trie& t, Position po, Length le, Sequence se ) const {
			t.diff( po, le, se );
		};
	} elementDiff;

	/**
	 * Functional returning the length of a subsequence to pass to the generic loop of the Trie
	 */
	struct LengthFunction {
		virtual Length operator() ( const Source& so, const Cover<Position>& co, Cover<Position>::iterator& r, Position p, Length minim, Length maxim ) const = 0;
	};

	const struct :public LengthFunction {
		virtual Length operator() ( const Source& so, const Cover<Position>& co, Cover<Position>::iterator& r, Position p, Length minim, Length maxim ) const {
// 	calling Cover<Position>::cover
			return co.cover( p, minim, maxim, r );
		};
	} lengthCover;

	const struct :public LengthFunction {
		virtual Length operator() ( const Source& so, const Cover<Position>& co, Cover<Position>::iterator& r, Position p, Length minim, Length maxim ) const {
// 	calling Range<Position>::cover
			return r->cover( p, minim, maxim );
		};
	} lengthRange;

	const struct :public LengthFunction {
		virtual Length operator() ( const Source& so, const Cover<Position>& co, Cover<Position>::iterator& r, Position p, Length minim, Length maxim ) const {
// 	using pre-calculated max_length_at (possibly filtered for melting temperatures)
			return so.max_length_at.at( p );
		};
	} lengthMax;

//=============================================================================================
// 	CODE
//=============================================================================================
public:
	Trie( Source& so, Length m, Length M );

	void cover( unsigned int threads );
	void touch( unsigned int threads );
	void filterHomolo( unsigned int threads, const Length max_homolo );
	void smallDiff( unsigned int threads );
	void confirm( unsigned int threads, const string& s, const vector<pair<string,Range<Position>>>& );

	void encodeClusters();
	void collectClusters( unsigned int threads );
	void collectMatches( unsigned int threads );
	void sortMatches( unsigned int threads );

	/**
	 * Build the array of TrieSlice based on prefixes encountered in the source
	 * Create the multimaps:
	 *  - prefix_match: ambiguous to unambiguous prefixes
	 *  - prefix_diff1: all prefixes with exactly one difference
	 */
	virtual void buildSlices();

	/**
	 * Return the names of all targets that have at least one match
	 */
	const set<string> getNodesWithMatches() const;

	/**
	 * \returns the cluster id associated with a site in the Source (EXPERIMENTAL)
	 * 
	 * WARNING: must be called AFTER collectClusters
	 */
	const Cluster getClusterId( Position p, Length l );

	/**
	 * Get the cluster id associated with a substring
	 * 
	 * \param s source string
	 * \param p start position
	 * \param l length of substring
	 * 
	 * \param clu cluster id (output)
	 * 
	 * \returns false if no cluster was found
	 */
	const bool getCluster( const string& s, Position p, Length l, Cluster& clu );

// 	/**
// 	 * Finds subsequence without ambiguities in Trie
// 	 * 
// 	 * WARNING: \param s must not have any ambiguous bases
// 	 */
// 	vector<Sequence> find( string s );
	virtual void _find( const string& s ) const;
	void _find1( const string& s ) const;
	void find( const vector<string>& s ) const;
	void find1( const vector<string>& s ) const;

	void measure(
		int& nodes, int& leaves, int& length, int& occurrences, int& clusters,
		Distribution& nucleo_distribution, Distribution& prefix_distribution, Distribution& depth_distribution, Distribution& length_distribution, Distribution& occurrence_distribution, Distribution& cluster_distribution
	) const;

protected:
	virtual void _cover( bool& );
	virtual void _touch( bool& );
	virtual void _filterHomolo( bool&, const Length& max_homolo );
	virtual void _smallDiff( bool& );
	virtual void _confirm( bool&, const string&, const vector<pair<string,Range<Position>>>& );
	virtual void __confirm( const string&, Sequence re, Position p, Length l );
	virtual void _collectClusters( bool& );
	virtual void _collectMatches( bool& );
	virtual void _sortMatches( bool& );

	void add( Position p, Length l, Sequence s );
	virtual void mark( Position p, Length l, Sequence s );
	virtual void diff( Position p, Length l, Sequence s );

	void _mark( Position p, Length l, Sequence s, Depth d, Slice k );

	/**
	 * Get the trie slice associated with the prefix found at \param p
	 */
	inline TrieSlice& getSlice( Position p ) {
		return getSlice( source.getSource(), p );
	}

	/**
	 * Get the trie slice associated with the prefix found in string \param s at \param p
	 */
	inline TrieSlice& getSlice( const string& s, Position p ) {
		return cake.at( prefixes.at( nu2p4( s, p )));
	}

	/**
	 * Pair of iterators describing all slices that match the prefix at position \param p in the source
	 */
	pair<unordered_multimap<Prefix4,Slice>::const_iterator,unordered_multimap<Prefix4,Slice>::const_iterator> getSlicesMatching( Position p ) const {
		return prefix_match.equal_range( nu2p4( source.getSource(), p ));
	};

	/**
	 * Pair of iterators describing all slices that have exactly one difference compared to the prefix at position \param p in the source
	 */
	pair<unordered_multimap<Prefix4,Slice>::const_iterator,unordered_multimap<Prefix4,Slice>::const_iterator> getSlicesDiff1( Position p ) const {
		return prefix_diff1.equal_range( nu2p4( source.getSource(), p ));
	};

	/**
	 * Create a new TrieSlice
	 * WARNING: not multi-threaded !
	 */
	void newSlice( Position p );

	/**
	 * Generic multithread-safe loop thay will execute an ElementFunction (functional)
	 * for every subsequence of a Cover returned by a CoverFunction (functional)
	 * of length returned by a LengthFunction (functional)
	 * 
	 * NOTE: the same loop is used by add, mark and diff for Trie and TrieAmbig
	 */
	void loop( bool& first, const CoverFunction& cov, const LengthFunction& len, const ElementFunction& ele );

//=============================================================================================
// 	TEST
//=============================================================================================
	void show() const ;
	void print() const;
	void verify() const;
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
