#ifndef __TrieSlice_h__
#define __TrieSlice_h__

#include <iostream>

#include <mutex>
#include <atomic>

#include <unordered_map>
#include <vector>
#include <map>

using namespace std;

#include "util.h"

#include "Types.h"
#include "Error.h"
#include "Source.h"

#include "Trie.h"

class TrieSlice {
//=======================================
// 	DATA
//=======================================
private:
	const Depth fixed_depth;

	unordered_map<Node, PositionLength>  source;
	unordered_map<SymbolNode, Node>      children_hash;
	unordered_multimap<Node, SymbolNode> children_list;

	unordered_multimap<Node,Sequence>    occurrences; // temporary; removed after collect

	unordered_map<Node, Cluster>  cluster; // id of the cluster of sequences that correspond to the node

	Node next_node;
	mutex lock;

//=======================================
// 	CODE interface
//=======================================
public:
	TrieSlice( Depth depth = Trie::fixed_depth );

	inline Depth getDepth() const { return fixed_depth; };

	/**
	 * Add subsequence to TrieSlice
	 * 
	 * Simple thread-safe wrapper around thread-unsafe recursive worker
	 * 
	 * \see TrieSlice::_add
	 */
	inline void add( Source& src, Sequence s, Position p, Length l, Length minim ) {
		lock.lock();
		_add( src, s, 0, p+fixed_depth, fixed_depth, l-fixed_depth, minim );
		lock.unlock();
	};

	inline void mark( Source& src, Sequence s, Position p, Length l, Length minim ) {
		lock.lock();
		_mark( src, s, 0, p+fixed_depth, fixed_depth, l-fixed_depth, minim );
		lock.unlock();
	};

	inline void confirm( Trie& trie, const Source& src, const string& s, Sequence re, Position p, Length l, Length minim ) {
		lock.lock();
		_confirm( trie, src, s, 0, re, p+fixed_depth, fixed_depth, l-fixed_depth, minim );
		lock.unlock();
	};

	/**
	 * Calls recursive _smallDiff for the 
	 * 
	 * \param diff 0 no difference detected so far
	 */
	inline void smallDiff( Source& src, Position p, Length l, Sequence s, Length minim, Length diff ) {
		lock.lock();
		if( diff ) {
// 	void _smallDiff( Source& src, Node n0, Position p   , Depth d    ,      Length l, Sequence s, Length minim  );
			_smallDiff1(     src,       0, p+fixed_depth, fixed_depth, l-fixed_depth,          s,        minim );
		} else {
// 	void _smallDiff( Source& src, Node n0, Position p   , Depth d    ,      Length l, Sequence s, Length minim );
			_smallDiff0(     src,       0, p+fixed_depth, fixed_depth, l-fixed_depth,          s,        minim );
		}
		lock.unlock();
	};

	/**
	 * Remove nodes in the Trie that have more than \param max_homolo
	 * homologous [identical sequential] bases (thread safe)
	 * 
	 * NOTE: calls recursive version \see filterHomolo
	 */
	void filterHomolo( Source& src, Depth max_homolo, Symbol sy0, Length h0, Length minim );

	/**
	 * Collect all clusters (sets of matching sequences) from all nodes in the trie slice, recursively
	 */
	void collectClusters( Trie& trie, Depth d );

	/**
	 * Generate a unique identifier for each Cluster
	 */
	void encodeClusters( set<set<Sequence>>& cluster_set, Depth d, Node n ) const;

	/**
	 * Collect all matches (
	 */
	void collectMatches( Trie& trie, Depth d );

//=======================================
// 	CODE workers
//=======================================
private:
	void  _add( Source& src, Sequence s, Node n0, Position p, Depth d, Length l, Length minim );
	void _mark( Source& src, Sequence s, Node n0, Position p, Depth d, Length l, Length minim );

	void _smallDiff0( Source& src, Node n0, Position p, Depth d, Length l, Sequence s, Length minim );
	void _smallDiff1( Source& src, Node n0, Position p, Depth d, Length l, Sequence s, Length minim );

	void _confirm( Trie& trie, const Source& src, const string& s, Node n0, Sequence re, Position p, Depth d, Length l, Length minim );

	/**
	 * Remove, recursively, nodes in the Trie that have more than \param max_homolo
	 * homologous [identical sequential] bases
	 * 
	 * \param max_homolo maximum length of a contiguous homologous section
	 * \param n0 starting node for the search
	 * \param d0 depth of node n0
	 * \param sy0 repeated symbol up to hitting node n0
	 * \param h0 accumulator (count) of contiguous section of symbol sy0 exactly before (above) n0
	 * \param minim minimum length of an oligo
	 * 
	 * NOTE: filterHomolo will erase (make inaccessible) whole chunks of the TrieSlice
	 */
	void _filterHomolo( Source& src, Depth max_homolo, Node n0, Depth d0, Symbol sy0, Length h0, Length minim );

	/**
	 * Collect all clusters (sets of matching sequences) from all nodes in the trie slice, recursively
	 */
	void _collectClusters( Trie& trie, Depth d, Node n );

	void _collectMatches( Trie& trie, deque<pair<Cluster,PositionDepthLength>>& m, Depth d, Node n );

//=======================================
// 	CODE inline workers
//=======================================
	inline bool childAt( Node n0, Symbol sy, Node& n ) const {
		auto pn = children_hash.find( symbolNode( sy, n0 ));
		if( pn == children_hash.end()) {
			return false;
		}
		n = pn->second;
		return true;
	};

	inline vector<SymbolNode> children( Node n0 ) const {
		vector<SymbolNode> r;
		auto it = children_list.equal_range( n0 );
		for( auto& c = it.first ; c != it.second ; c++ ) {
			r.emplace_back( c->second );
		}
		return r;
	};

	inline vector<Sequence> occ( Node n0 ) const {
		vector<Sequence> r;
		auto it = occurrences.equal_range( n0 );
		for( auto& c = it.first ; c != it.second ; c++ ) {
			r.emplace_back( c->second );
		}
		return r;
	};

	inline Node newChild( Node n0, Symbol sy, Position p, Length l ) {
		Node n = next_node++;

		source.emplace( n, positionLength( p, l ));
		children_hash.emplace( symbolNode( sy, n0 ), n );
		children_list.emplace( n0, symbolNode( sy, n ));

		return n;
	};

	/**
	 * Split an existing node into two nodes, duplicating occurrences
	 * 
	 *           l     l0
	 * n0 |------======|E[...]    =>   n0 |------|[...]
	 *                                        n1  +|======|E[...]
	 * 
	 * \return new node
	 */
	inline Node splitNode( Node n0, Symbol sy, Length l ) {
		PositionLength pl0 = source.at( n0 );
		Position p0 = position( pl0 );
		Length l0 = length( pl0 );

		assert( l  > 0 );
		assert( l0 > l );

		vector<SymbolNode> ch = children( n0 );
		vector<Sequence>   oc = occ( n0 );

		source.at( n0 ) = positionLength( p0, l );

		eraseChildren( n0 );
		Node n1 = newChild( n0, sy, p0+l, l0-l );

		for( SymbolNode sn: ch ) {
			children_list.emplace( n1, sn );
			children_hash.emplace( symbolNode( symbol( sn ), n1 ), node( sn ));
		}

		for( Sequence s: oc ) {
			occurrences.emplace( n1, s );
		}

		return n1;
	}

	/**
	 * Erase a child with symbol \param sy of node \param n0
	 * 
	 * Only the navigation elements in the children_* maps are erased:
	 * the payload remains, but is inaccessible
	 * 
	 * TODO: erase descendents
	 */
	inline void eraseChild( Node n0, Symbol sy ) {
		children_hash.erase( symbolNode( sy, n0 ));

// 	ch is a pair<iterator_first, iterator_past_last>
		auto ch = children_list.equal_range( n0 );
// 	c is a pair<Node,SymbolNode>
		for( auto c=ch.first ; c!=ch.second ; c++ ) {
			if( symbol( c->second ) == sy ) {
				auto n=c;
// 	n is the next element after the child to be erased
				children_list.erase( c, ++n );
// 	At this point, the job is done; Leave now!
// 	WARNING: the iterators are invalidated after "erase"
				return;
			}
		}
	};
	/**
	 * Erase all children of node \param n0
	 * 
	 * Only the navigation elements in the children_* maps are erased:
	 * the payload remains, but is inaccessible
	 * 
	 * TODO: erase descendents
	 */
	inline void eraseChildren( Node n0 ) {
// 	ch is a pair<iterator_first, iterator_past_last>
		auto ch = children_list.equal_range( n0 );
// 	c is a pair<Node,SymbolNode>
		for( auto c=ch.first ; c!=ch.second ; c++ ) {
			children_hash.erase( symbolNode( symbol( c->second ), n0 ));
		}
		children_list.erase( n0 );
	};

	/**
	 * Erase all occurrences of a given node, unless all are of Sequence s
	 * 
	 * Useful for with --diff=no to avoid erasing a oligonucleotide if a subsequence
	 * with a small difference occurs in the same sequence.
	 */
	inline void eraseOccurrencesUnlessOwn( Node n0, Sequence s ) {
		auto r = occurrences.equal_range( n0 );

		for( auto se = r.first ; se != r.second ; se++ ) {
			if( se->second != s ) {
				occurrences.erase( n0 );
// 	WARNING: invalidates the iterator; don't do anything more with it!
				return;
			}
		}
	};

//=======================================
// 	TEST
//=======================================
public:
	void show( const Source& src, Node n, Depth d ) const;
	void measure(
		Node n0, Depth d,
		int& nodes, int& leaves, int& l, int& o,
		Distribution& depth_distribution, Distribution& length_distribution, Distribution& occurrence_distribution
	) const;

	void print( const Source& src, Node n, Depth d, string root="" ) const;
	string getNodeSource( Source& src, Node n0, Node n, Depth d=0 ) const;
	void verify( const Source& src, Node n ) const;

	inline Node size() const { return children_hash.size(); };

	unordered_multimap<string,Sequence>& find( const Source& so, const string& s, unordered_multimap<string,Sequence>& a ) const;
	unordered_multimap<string,Sequence>& find( const Source& so, const string& s, Node n0, Depth d, unordered_multimap<string,Sequence>& a ) const;

	/**
	 * \returns cluster id associated with oligo candidate at position \param p of length \param l
	 */
	Cluster getClusterId( Source& src, Position p, Length l ) const;
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
	const bool getCluster( const Trie& trie, const string& s, Position p, Length l, Cluster& clu ) const;

private:
	/**
	 * Recursive worker function
	 * 
	 * \returns cluster id associated with oligo candidate at position \param p of length \param l
	 */
	Cluster _getClusterId( Source& src, Node n0, Position p, Depth d, Length l ) const;

	const bool _getCluster( const Trie& trie, Node n0, Depth d, const string& s, Position p, Length l, Cluster& clu ) const;

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
