#define __Trie_cpp__

#include "Trie.h"
#include "TrieSlice.h"

Trie::Trie( Source& so, Length m, Length M ) : source( so ), minim( m ), maxim( M ), cake() {
};

void Trie::buildSlices() {
	Slice cake_index( 0 );

	for( Symbol n1: { 1, 2, 4, 8 }) {
		for( Symbol n2: { 1, 2, 4, 8 }) {
			for( Symbol n3: { 1, 2, 4, 8 }) {
				for( Symbol n4: { 1, 2, 4, 8 }) {
					prefixes.emplace( nu2p4( n1, n2, n3, n4 ), cake_index++ );
				}
			}
		}
	}

	for( const auto& e2fr: source.fragments.to ) {
		for( const Range<Position>& r: e2fr.first.getAmbigCompl()) { // capture unambiguous prefixes
			if( r.size() < minim ) continue; // range will not contribute to the creation of slices

			const Position first = r.lo();
			const Position last = r.hi() - minim;

			for( Position p = first ; p <= last ; ++p ) {
				newSlice( p );
			}
		}
	}

	cake.resize( prefixes.size());

	set<Prefix4> ap4; // temporary storage for ambiguous prefixes

	for( const auto& e2fr: source.fragments.to ) {
		Cover<Position> c = e2fr.first.getAmbig();
		c += ( fixed_depth-1 );

		for( const Range<Position>& r: c ) {
			Position p = r.lo();
			while( true ) {
				Length l = r.cover( p, fixed_depth, fixed_depth );
				if( !l ) {
					break;
				}

				ap4.emplace( nu2p4( source.getSource(), p++ ));
			}
		}
	}

// 	populate the ambiguous prefix match
	for( const auto& e2pr: prefixes ) { // unambiguous prefixes match only themselves
		prefix_match.emplace( e2pr.first, e2pr.second );
	}

	for( const Prefix4& pr1: ap4 ) {  // ambiguous prefixes match any unambiguous prefixes
		for( const auto& e2pr2: prefixes ) {
			if( p4ma( pr1, e2pr2.first ) == 4 ) {
				prefix_match.emplace( pr1, e2pr2.second );
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
		for( const Prefix4& pr3: ap4 ) { // add all matching diff1 unambiguous prefixes
			if( p4ma( e2pr1.first, pr3 ) == 3 ) {
				prefix_diff1.emplace( pr3, e2pr1.second );
			}
		}
	}
};

void Trie::newSlice( Position p ) {
	static Slice cake_index = 0; // incremental index

	Prefix4 pr = nu2p4( source.getSource(), p );
	if( prefixes.emplace( pr, cake_index ).second ) {
// 	here, the insertion has actually taken place
		++cake_index;
	}
};

const set<string> Trie::getNodesWithMatches() const {
	set<string> result;

	for( const auto& e2ma: matches ) {
		if( !source.targets.has( source.clusters.at( e2ma.first ))) continue;
		result.emplace( source.targets.at( source.clusters.at( e2ma.first )));
	}
	return result;
};

/**
 * Load all subsequences from the associated "source" database into the Trie; multi-threaded
 */
void Trie::cover( unsigned int threads ){
	spin( threads, *this, &Trie::_cover );
}

void Trie::_cover( bool& first ){
	loop( first, ambigCoverComplement, lengthMax, elementAdd );
	return;
}

/**
 * Add a subsequence starting at to the corresponding TrieSlice
 * 
 * \param p is the position in the single string source
 */
void Trie::add( Position p, Length l, Sequence s ) {
	getSlice( p ).add( source, s, p, l, minim );
};

/**
 * Mark all ambiguous subsequences from the associated "source" database into the Trie
 */
void Trie::touch( unsigned int threads ){
	spin( threads, *this, &Trie::_touch );
}

void Trie::_touch( bool& first ){
	loop( first, ambigCover, lengthCover, elementMark );
}

/**
 * Marks occurrences of an ambiguous subsequence in the Trie
 * 
 */
void Trie::mark( Position p, Length l, Sequence s ) {
	for( auto i2sl = this->getSlicesMatching( p ) ; i2sl.first != i2sl.second ; ++i2sl.first ) {
		cake.at( i2sl.first->second ).mark( source, s, p, l, minim );
	}
};

void Trie::smallDiff( unsigned int threads )
{
	spin( threads, *this, &Trie::_smallDiff );
}

void Trie::_smallDiff( bool& first )
{
	loop( first, range, lengthRange, elementDiff );
}

/**
 * Wrapper around recursive function _diff
 */
void Trie::diff( Position p, Length l, Sequence s )
{
	for( auto i2pr = prefix_match.equal_range( nu2p4( source.getSource(), p ))
			; i2pr.first != i2pr.second ; ++i2pr.first ) { // look at matching prefixes
		cake.at( i2pr.first->second ).smallDiff( source, p, l, s, minim, 0 );
	}

	for( auto i2pr = prefix_diff1.equal_range( nu2p4( source.getSource(), p ))
			; i2pr.first != i2pr.second ; ++i2pr.first ) { // look at prefixes with small difference
		cake.at( i2pr.first->second ).smallDiff( source, p, l, s, minim, 1 );
	}
}

/**
 * Confirm the trie against a vector of reference sequences
 * 
 * \param s : storage for the content of all reference sequences
 * \param a : vector of name ; range within the storage string
 */
void Trie::confirm( unsigned int threads, const string& s, const vector<pair<string,Range<Position>>>& a ) {
	if( !s.size() || !a.size()) return; // if there is no content

	spin( threads, *this, &Trie::_confirm, s, a );  // simply spins the _confirm worker function
}
/**
 * Iterate through all subsequences of a vector of reference strings \param a
 */
void Trie::_confirm( bool& first, const string& s, const vector<pair<string,Range<Position>>>& a ) {
	static vector<pair<string,Range<Position>>>::const_iterator i, e;

	static Position p;
	static mutex lock;

	lock.lock();

	if( first ) { // initialization
		first = false;

		i = a.cbegin(); // iterator for all reference sequences
		e = a.cend(); // end of the vector of range sequences

		if( i == e ) goto last_element; // the vector of references may have been empty (defence)
		p = i->second.lo(); // position at the start of the sequence
	} else {
		if( i == e ) goto last_element; // some other thread may have finished processing
	}

	while( true ) { // sequence loop
		while( true ) { // position loop
			Length le = i->second.cover( p, minim, maxim ); // length '0' means the 

			if( !le ) break; // to sequence loop

			Position po = p;
			Sequence re = source.reference.at( i->first );

			lock.unlock();

			__confirm( s, re, po, le ); // call the worker for each subsequence

			lock.lock();
			if( i == e ) goto last_element; // the sequence counter may have been changed by other threads outside the lock
			++p;
		}

		assert( i != e );

		if( ++i == e ) goto last_element;

		p = i->second.lo();
	}


last_element:
	lock.unlock();
}

/**
 * Confirm all matching slices of the trie against a subsequence of a reference sequence
 */
void Trie::__confirm( const string& s, Sequence re, Position p, Length l ) {
	for( auto i2pr = prefixes.equal_range( nu2p4( s, p )) ; i2pr.first != i2pr.second ; ++i2pr.first ) { // for each prefix matching the subsequence
		cake.at( i2pr.first->second ).confirm( *this, source, s, re, p, l, minim ); // call the worker of that trie slice
	}
}

/**
 * Collect all occurrences (clusters) from the Trie into a map of matches
 */
void Trie::collectMatches( unsigned int threads )
{
	spin( threads, *this, &Trie::_collectMatches );
}

/**
 * Worker function: collect all clusters from each TrieSlice into the map of matches
 * 
 * WARNING: filling Trie::matches (shared resource) must be done under lock
 */
void Trie::_collectMatches( bool& first ) {
	static decltype( cake )::iterator cr( cake.begin()), en( cake.end());
	static mutex lock;

	lock.lock();

	while( cr != en ) {
		TrieSlice& sl = *cr++;

		lock.unlock();
		sl.collectMatches( *this, sl.getDepth());
		lock.lock();
	}

	lock.unlock();
}

/**
 * Generates a unique, persistent identifier for each unique cluster (set of Sequences) contained
 * in Trie nodes
 * 
 * Since the unique identifier is the order in the [sorted] set of Sequence sets,
 * the cluster identifier is persistent
 * 
 * WARNING: single threaded
 */
void Trie::encodeClusters() {
	set<set<Sequence>> cluster_set;
	for( const TrieSlice& slice: cake ) { // collect all 
		slice.encodeClusters( cluster_set, slice.getDepth(), 0 );
	}

	assert( source.clusters.from.size() == 0 );

	int cluster_id = 0; // cluster counter
	for( const auto& cluster: cluster_set ) {
		source.clusters.emplace( cluster, cluster_id++ );
	}
}

/**
 * Collects and writes a unique cluster identifier in each node associated with it
 */
void Trie::collectClusters( unsigned int threads )
{
	spin( threads, *this, &Trie::_collectClusters );
}

/**
 * Worker function: collects and writes a unique cluster identifier in each node associated with it
 * by calling each slice's collectClusters method
 * 
 * NOTE: operations on each slice will be in isolated threads or sequentially
 * WARNING: incrementing the cluster id needs to be done under lock in TrieSlice::collectClusters
 */
void Trie::_collectClusters( bool& first ) {
	static mutex lock;

	static decltype( cake )::iterator cr( cake.begin()), en( cake.end());

	lock.lock();

	while( cr != en ) {
		TrieSlice& slice = *cr++;

		lock.unlock();
		slice.collectClusters( *this, slice.getDepth());
		lock.lock();
	}
	lock.unlock();
}

/**
 * Sort all matches by the sequence and the position where they occur
 */
void Trie::sortMatches( unsigned int threads )
{
	spin( threads, *this, &Trie::_sortMatches );
}

/**
 * Worker function: sort all matches by the sequence and the position where they occur
 */
void Trie::_sortMatches( bool& first ) {
	static decltype( matches )::iterator cr( matches.begin()), en( matches.end());
	static mutex lock;

	lock.lock();
	while( cr !=en ) {
		decltype( cr->second )& li = cr->second;
		cr++;

		lock.unlock();
		::sort( li.begin(), li.end(), pdlCompare );
		lock.lock();
	}
	lock.unlock();
}

/**
 * Remove occurrences of subsequences that have self-homologous sections
 * longer than \param max_homolo
 */
void Trie::filterHomolo( unsigned int threads, const Length max_homolo )
{
	spin( threads, *this, &Trie::_filterHomolo, max_homolo );
}

/**
 * Worker thread: remove occurrences of subsequences that have self-homologous sections
 * longer than \param max_homolo
 */
void Trie::_filterHomolo( bool& first, const Length& max_homolo ) {
	static mutex lock;
	static decltype( prefixes )::iterator cr, en;

	lock.lock();

	if( first ) { // initialization
		cr = prefixes.begin();
		en = prefixes.end();

		first = false;
	}

	while( cr != en ) { // for each prefix/slice
// 	Local variables; use outside lock
		const auto e2ho = p4ho( cr->first ); // the homologous section of the prefix
		const Slice sl  = cr->second; // the slice associated with the prefix

		++cr;

		lock.unlock();

		if( e2ho.first > max_homolo ) { // the whole slice is busted
			cake.at( sl ).filterHomolo( source, max_homolo, 0, max_homolo+1, minim );
		} else { // call recursive filterHomolo for the slice
			cake.at( sl ).filterHomolo( source, max_homolo, e2ho.second, e2ho.first, minim );
		}

		lock.lock();
	}

	lock.unlock();
}

const Cluster Trie::getClusterId( Position p, Length l ) {
	return getSlice( p ).getClusterId( source, p, l );
}

const bool Trie::getCluster( const string& s, Position p, Length l, Cluster& clu ) {
	return getSlice( s, p ).getCluster( *this, s, p, l, clu );
}

void Trie::loop( bool& first, const CoverFunction& cov, const LengthFunction& len, const ElementFunction& ele ) {
// // 	Thread-shared counters
// // 	WARNING: any changes must be done under lock
	static decltype( source.instance_fragments.from )::const_iterator fr, en;
	static Cover<Position> c;
	static Cover<Position>::iterator r;
	static Position p ;

	static mutex lock;

	lock.lock();
// 	initialization
	if( first ) {
// 	Only the first thread will execute the initialization
		fr = source.instance_fragments.from.begin();
		en = source.instance_fragments.from.end();

		while( true ) {
			if( fr == en ) {
// 	exit condition
				goto last_sequence;
			}

			c = cov( *this, source.fragments.at( fr->second ));
			if( !!c ) {
				break;
			}

			++fr;
			continue;
		}

		r = c.begin();
// 	HERE: c is non-empty
		assert( !!c );

		p = r->lo();
		first = false;
	} else {
// 	For subsequent threads, make sure that there still is processing to do

// 	NOTE: the first threads may have finished processing the whole input before
// 	      the last threads get to initialization
		if( fr == en ) {
// 	repeat exit condition
			goto last_sequence; // to last_sequence
		}
	}

	assert( !!c );

// 	Loop
	while( true ) { // fragment loop
		while( true ) { // range loop
			while( true ) { // position loop
				Length le = len( source, c, r, p, minim, maxim );

				if( !le ) {
					break;
				}

				Position po = p++; // increment the position iterator
				Sequence se = fr->first;

				lock.unlock();

// 	Communicate with the cover function using local, thread-specific variables po le se
				ele( *this, po, le, se );

// 	NOTE: position(p), range (r) and sequence (s) may have been changed magically by other threads outside the lock
				lock.lock();

				if( fr == en ) {
// 	repeat exit condition
					goto last_sequence; // to last_sequence
				}

				if( r == c.end()) {
					break;
				}
			} // position loop


			if(( r == c.end()) || ( ++r == c.end())) {
				break; // to fragment loop
			}

			p = r->lo();
		} // range loop

		while( true ) { // increment sequence
			if( ++fr == en ) {
// 	exit condition
				goto last_sequence; // to last_sequence
			}

			c = cov( *this, source.fragments.at( fr->second ));
			if( !c ) {
				continue;
			}

			r = c.begin();
			assert( r != c.end());
			p = r->lo();

			assert( *r & p );
			break;
		}
	}  // fragment loop

last_sequence:
	lock.unlock();
};

// TEST

/**
 * Get the human readable prefix from a trie slice
 */
void Trie::show() const {
	for( auto& slice: cake ) {
		slice.show( source, 0, Depth( fixed_depth ));
	}
	cout << endl;
}

void Trie::measure(
	int& nodes, int& leaves, int& length, int& occurrences, int& clusters,
	Distribution& nucleo_distribution, Distribution& prefix_distribution,
	Distribution& depth_distribution, Distribution& length_distribution, Distribution& occurrence_distribution, Distribution& cluster_distribution
) const {
	for( Slice slice=0 ; slice<cake.size() ; slice++ ) {
		cake.at( slice ).measure(
			0, cake.at( slice ).getDepth(),
			nodes, leaves, length, occurrences,
			depth_distribution, length_distribution, occurrence_distribution
		);
	}

	for( const auto& ipr: prefixes ) {
		prefix_distribution[ ipr.first ] = cake.at( ipr.second ).size();
	}

	Position l = source.length();
	for( Position p=0 ; p<l ; p++ ) {
		nucleo_distribution[ source.symbol( p )]++;
	}
}

void Trie::print() const {
	for( const auto& sl: cake ) {
		sl.print( source, 0, Depth( fixed_depth ));
	}
}

void Trie::verify() const {
	for( const auto& sl: cake ) {
		sl.verify( source, 0 );
	}
}

void Trie::find( const vector<string>& v ) const {
	for( const string& s: v ){
		_find( s );
	}
}

void Trie::_find( const string& s ) const {
	string nu = convertAsc2Nu( s );
	unordered_multimap<string,Sequence> r;

	cout << "--------------------------" << endl;
	cout << s << endl;
	assert( nu.size() >= fixed_depth );

	Prefix4 pr = nu2p4( nu, 0 );

	cake.at( prefixes.at( pr )).find( source, nu, r );

	for( auto i2pr = prefix_match.equal_range( pr ) ; i2pr.first != i2pr.second ; ++i2pr.first ) {
		cake.at( i2pr.first->second ).find( source, nu, r );
	}

	for( const auto& k: r  ) {
		cout << k.first << '\t' << source.getSequenceName( k.second ) << endl;
	}
}

void Trie::_find1( const string& s ) const {
	cout << "==========================" << endl;
	cout << s << endl;

	_find( s );

	string b = s;
	for( Depth d = 0 ; d < s.size() ; d++ ) {
		b=s;
		if( asc2nu.at( b[d]) == 0xF ) { // 'N' in position d
			continue;
		}

		b[d] = 'N';
		_find( b );
	}
}

void Trie::find1( const vector<string>& v ) const {
	for( const string& s: v ){
		_find1( s );
	}
}

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
