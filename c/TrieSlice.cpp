#define __TrieSlice_cpp__

#include "TrieSlice.h"

TrieSlice::TrieSlice( Depth d ) : fixed_depth( d ), next_node( 1 ) {
	source.emplace( 0, positionLength( 0, 0 ));
}

void TrieSlice::filterHomolo( Source& src, Depth max_homolo, Symbol sy0, Length h0, Length minim ) {
	lock.lock();
	if( h0 > max_homolo ) {
		eraseChildren( 0 );
	} else {
		_filterHomolo( src, max_homolo, 0, fixed_depth, sy0, h0, minim );
	}
// 	TODO: clear/rehash unsorted maps
	lock.unlock();
}

void TrieSlice::_filterHomolo( Source& src, Depth max_homolo, Node n0, Depth d0, Symbol sy0, Length h0, Length minim ) {
	for( SymbolNode sn: children( n0 )) {
		Depth  h  = h0;
		Symbol sy = sy0;
		Node   n  = node( sn );

		PositionLength pl = source.at( n );
		Position       p  = position( pl );
		Length         l  = length( pl );

		for( Depth d = 0 ; d < l ; d++ ) {

			if( sy == src.symbol( p+d )) {
				if( ++h > max_homolo ) {
					if(( d0+d ) < minim ) {
// 	erase
						eraseChild( n0, symbol( sn ));
					} else {
// 	shrink
						source.at( n ) = positionLength( p, d );
						eraseChildren( n );
					}

					goto next_child;
				}

				continue;
			}

			sy = src.symbol( p+d );
			h = 1;
		}

		_filterHomolo( src, max_homolo, n, d0+l, sy, h, minim );

next_child:
		;
	}
}

void TrieSlice::collectClusters( Trie& trie, Depth d ) {
	_collectClusters( trie, d, 0 ); // call recursive worker function

// 	Cleanup list of occurrences
	occurrences.clear();
	occurrences.rehash( occurrences.size());
}

void TrieSlice::_collectClusters( Trie& trie, Depth d, Node n ) {
	for( auto& c: children( n )) {
		_collectClusters( trie, d+length( source.at( n )), node( c ));
	}

	const auto o = occ( n );
	if( o.empty()) return; // no occurrences

	const set<Sequence> s{ o.begin(), o.end() };
	cluster.emplace( n, trie.source.clusters.at( s )); // add cluster id to this node
}

void TrieSlice::encodeClusters( set<set<Sequence>>& cluster_set, Depth d, Node n ) const {
	for( auto& c: children( n )) {
		encodeClusters( cluster_set, d+length( source.at( n )), node( c ));
	}

	const auto o = occ( n );
	if( o.empty()) return; // no occurrences

	const set<Sequence> s{ o.begin(), o.end() };

	cluster_set.emplace( s );
}

void TrieSlice::collectMatches( Trie& trie, Depth d ) {
	deque<pair<Cluster,PositionDepthLength>> m;

	_collectMatches( trie, m, d, 0 );

	static mutex lock;
	lock.lock(); // protect unique Trie::matches from multithreaded access
	for( const auto& e2ma : m ) {
		trie.matches[ e2ma.first ].push_back( e2ma.second );
	}
	lock.unlock();
}

void TrieSlice::_collectMatches( Trie& trie, deque<pair<Cluster,PositionDepthLength>>& m, Depth d, Node n ) {
	Position p = position( source.at( n ));
	Length l   = length( source.at( n ));

	for( auto& c: children( n )) {
		_collectMatches( trie, m, d+l, node( c ));
	}

	if( !cluster.count( n )) return; // no occurrences

	m.emplace_back( cluster.at( n ), positionDepthLength( p, d, l ));
// 		trie.matches[ trie.source.slice_cluster.at( cluster.at( n ))].push_back( positionDepthLength( p, d, l )); // add this node's cluster id to the list of matches
}

// |ACGTA|                 (p   ;l=5) <-- array<TrieSlice,fixed_depth> from Trie
//        +|AG|            (p+5 ;l=2) <-- TrieSlice::add starts here; calls _add recursively
//             +|CGA|      (p+7 ;l=3) 
//                   +|TT| (p+10;l=2)
/**
 * Add subsequence to a node in the Trie, recursively
 * 
 * \param d is the depth of the top of the node (recursive parameter), since the depth
 * is not stored in the Trie, just inferred from traversal
 */
// TODO: TrieSlice::_add remove debug (LOW)
void TrieSlice::_add( Source& src, Sequence s, Node n0, Position p, Depth d, Length l, Length minim ) {
	Position p0 = position( source.at( n0 ));
	Length l0   = length( source.at( n0 ));

	Depth dd;
	bool mismatch = false;

	for( dd = 0 ; dd < min( l0, l ) ; dd++ ){
		if( src.symbol( p0+dd ) == src.symbol( p+dd )){
			continue;
		}

		mismatch = true;
		break;
	}

	if( !mismatch && ( l0 == l )) {
// (1)  n0 |-----------|E [s0,...]     =>   n0 |-----------|E [s0,...,s]
//      s  |-----------|

// 	NOTE: since the subsequence s, p, d, l has at least "minimum" length,
// 	it is safe to add s to the occurrences
		occurrences.emplace( n0, s );

		if( p<p0 ) {
// 	point the source to the earlied occurrence in the database
			source.at( n0 ) = positionLength( p, l0 );
		}

		return;
	}

	if( !mismatch && ( l0 < l )) {
// (2)  n0 |-----------|E [s0,...]     =>   n0 |-----------|E [s0,...,s]
//      s  |-----------====|                                +|====|   [s]


		Symbol sy = src.symbol( p+dd );
		Node n;

		if(( d+dd ) >= minim ) {
// 	only add occurrences to the nodes below the "minimum" depth
			occurrences.emplace( n0, s );
		}

		if( p<p0 ) {
// 	point the source to the earlied occurrence in the database
			source.at( n0 ) = positionLength( p, l0 );
		}

		if( childAt( n0, sy, n )){
			_add( src, s, n, p+dd, d+dd, Length( l-dd ), minim );
			return;
		}

		n = newChild( n0, sy, p+dd, Length( l-dd ));

		occurrences.emplace( n, s );

		return;
	}

// 	TODO: TrieSlice::_add (3) inherit children (MEDIUM)
	if( !mismatch ) {
// (3)  n0 |-------====|E [s0,...]     =>   n0 |-------|   [s0,...,s]
//      s  |-------|                                    +|====|E [s0,...]

		Symbol sy = src.symbol( p0+dd );
		vector<SymbolNode> ch = children( n0 );

		for( auto& c: ch ) {
			children_hash.erase( symbolNode( symbol( c ), n0 ));
		}
		children_list.erase( n0 );

// 	insert n1
		Node n = newChild( n0, sy, p0+dd, l0-dd );
		for( auto& c: ch ) {
			children_hash.emplace( symbolNode( symbol( c ), n ), node( c ));
			children_list.emplace( n, c );
		}

// 	shrink n0
		if( p<p0 ) {
// 	point the source to the earlied occurrence in the database
			source.at( n0 ) = positionLength( p, dd );
		} else {
			source.at( n0 ) = positionLength( p0, dd );
		}

// 		occurrences of n
		for( Sequence& ss: occ( n0 )) {
			occurrences.emplace( n, ss );
		}

		occurrences.emplace( n0, s );

		return;
	}

// 	TODO: TrieSlice::_add (4) inherit children (MEDIUM)
// (4)  n0 |-------====|E [s0,...]     =>   n0 |-------|   (if length>min) [s0,...,s] else []
//      s  |-------~~~|                              n1 +|====|E [s0,...]
//                                                   n2 +|~~~|  [s]

	Symbol sy1 = src.symbol( p0+dd );
	Symbol sy2 = src.symbol( p +dd );

// 	capture and erase children of node n0
	vector<SymbolNode> ch = children( n0 );
	for( auto& c: ch ) {
		children_hash.erase( symbolNode( symbol( c ), n0 ));
	}
	children_list.erase( n0 );

// 	insert n1
	Node n1 = newChild( n0, sy1, p0+dd, l0-dd );
	for( auto& c: ch ) {
		children_hash.emplace( symbolNode( symbol( c ), n1 ), node( c ));
		children_list.emplace( n1, c );
	}

// 	insert n2
	Node n2 = newChild( n0, sy2, p+dd, l-dd );

// 	shrink n0
	if( p<p0 ) {
// 	point the source to the earlied occurrence in the database
		source.at( n0 ) = positionLength( p, dd );
	} else {
		source.at( n0 ) = positionLength( p0, dd );
	}

// 	occurrences
// 		occurrences of n2
	occurrences.emplace( n2, s );

// 		occurrences of n1
	for( Sequence& ss: occ( n0 )) {
		occurrences.emplace( n1, ss );
	}

	if(( d+dd ) >= minim ) {
// 	keep n0 occurrences
		occurrences.emplace( n0, s );
	} else {
// 	remove n0 occurrences
		occurrences.erase( n0 );
	}

	return;
}

void TrieSlice::_mark( Source& src, Sequence s, Node n0, Position p, Depth d, Length l, Length minim ) {
	Position p0 = position( source.at( n0 ));
	Length l0   = length( source.at( n0 ));

	Depth dd;
	bool mismatch = false;

	for( dd = 0 ; dd < min( l0, l ) ; dd++ ){
// 	WARNING: matching between ambiguous symbols ("&" and not "=")
		if( src.symbol( p0+dd ) & src.symbol( p+dd )){
			continue;
		}

		mismatch = true;
		break;
	}

// 	cout << dd << " " << l0 << " " << l << endl;

	if( !mismatch && ( l0 == l )) {
// (1)  n0 |-----------|E [s0,...]     =>   n0 |-----------|E [s0,...,s]
//      s  |-----------|

// 	NOTE: since the subsequence s, p, d, l has at least "minimum" length,
// 	it is safe to add s to the occurrences
		occurrences.emplace( n0, s );
		return;
	}

	if( !mismatch && ( l0 < l )) {
// (2)  n0 |-----------|E [s0,...]     =>   n0 |-----------|E [s0,...,s]
//      s  |-----------====|                                +[mark matching children]

		Symbol sy = src.symbol( p+dd );

		if(( d+dd ) >= minim ) {
// 	only add occurrences to the nodes below the "minimum" depth
			occurrences.emplace( n0, s );
		}

		for( auto& c: children( n0 )) {
			if( sy & symbol( c )) {
				_mark( src, s, node( c ), p+dd, d+dd, Length( l-dd ), minim );
			}
		}

		return;
	}

// 	TODO: TrieSlice::_mark (3) inherit children (MEDIUM)
	if( !mismatch ) {
// (3)  n0 |-------====|E [s0,...]     =>   n0 |-------|   [s0,...,s]
//      s  |-------|                                    +|====|E [s0,...]

		Symbol sy = src.symbol( p0+dd );
		vector<SymbolNode> c0 = children( n0 );

		for( auto& c: c0 ) {
			children_hash.erase( symbolNode( symbol( c ), n0 ));
		}
		children_list.erase( n0 );

// 	insert n1
		Node n = newChild( n0, sy, p0+dd, l0-dd );
		for( auto& c: c0 ) {
			children_hash.emplace( symbolNode( symbol( c ), n ), node( c ));
			children_list.emplace( n, c );
		}

// 	shrink n0
		source.at( n0 ) = positionLength( p0, dd );

// 		occurrences of n
		for( Sequence& ss: occ( n0 )) {
			occurrences.emplace( n, ss );
		}

		occurrences.emplace( n0, s );
		return;
	}

// 	TODO: TrieSlice::_mark (4) inherit children (MEDIUM)
// (4)  n0 |-------====|E [s0,...]     =>   n0 |-------|   (if length>min) [s0,...,s] else []
//      s  |-------~~~|                              n1 +|====|E [s0,...]

	if(( d+dd ) < minim ) {
		return;
	}

// 	only add occurrences to the nodes below the "minimum" depth
	Symbol sy = src.symbol( p0+dd );
	vector<SymbolNode> c0 = children( n0 );

	for( auto& c: c0 ) {
		children_hash.erase( symbolNode( symbol( c ), n0 ));
	}
	children_list.erase( n0 );

// 	insert n1
	Node n = newChild( n0, sy, p0+dd, l0-dd );
	for( auto& c: c0 ) {
		children_hash.emplace( symbolNode( symbol( c ), n ), node( c ));
		children_list.emplace( n, c );
	}

// 	shrink n0
	source.at( n0 ) = positionLength( p0, dd );

// 		occurrences of n
	for( Sequence& ss: occ( n0 )) {
		occurrences.emplace( n, ss );
	}

	occurrences.emplace( n0, s );
	return;
}

/**
 * Remove from a Trie occurrences of subsequences with exactly one bp difference compared
 * with a source subsequence, recursively. Exactly one bp difference has been detected before
 * the current symbol in the source subsequence.
 * 
 * \param n0 current node (target subsequence)
 * 
 * \param p position of the current symbol in the source subsequence
 * \param l remaining length of the current subsequence
 * 
 * \param d depth of the current node and the current symbol in the source subsequence
 * \param minim minimum length of an oligonucleotide
 */
void TrieSlice::_smallDiff1( Source& src, Node n0, Position p, Depth d, Length l, Sequence s, Length minim ) {
	PositionLength pl0 = source.at( n0 );
	Position p0 = position( pl0 );
	Length l0   = length( pl0 );

	Depth dd = 0;

	for( dd = 0 ; dd < min( l0, l ) ; dd ++ ) {
		if( src.symbol( p0+dd ) & src.symbol( p+dd )) {
			continue;
		}

// 	HERE: more than one difference
		if(( dd > 0 ) && ( d+dd ) >= minim ) {

// 	split if bottom lower than minim
			splitNode( n0, src.symbol( p0+dd ), dd );
			eraseOccurrencesUnlessOwn( n0, s );
		}

		return;
	}

// 	HERE: no new differences; exactly one difference
	if( l0 == l ) {

// 	HERE: source subsequence ends at the end of the node
		eraseOccurrencesUnlessOwn( n0, s );
		return;
	}

	if( l0 > l ) {

// 	HERE: source subsequence ends before current node
		splitNode( n0, src.symbol( p0+l ), l );
		eraseOccurrencesUnlessOwn( n0, s );
		return;
	}

	assert( l0 == dd );

	if(( d+l0 ) >= minim ) {
		eraseOccurrencesUnlessOwn( n0, s );
	}

// 	HERE: more to compare; regression
	Symbol sy = src.symbol( p+l0 );

	for( SymbolNode sn: children( n0 )) {
		if( sy & symbol( sn )) {
			_smallDiff1( src, node( sn ), p+l0, d+l0, l-l0, s, minim );
		}
	}
}

/**
 * Remove from a Trie occurrences of subsequences with exactly one bp difference compared
 * with a source subsequence, recursively. Exactly zero differences have been detected before
 * the current symbol in the source subsequence.
 * 
 * \param n0 current node (target subsequence)
 * 
 * \param p position of the current symbol in the source subsequence
 * \param l remaining length of the current subsequence
 * 
 * \param d depth of the current node and the current symbol in the source subsequence
 * \param minim minimum length of an oligonucleotide
 */
void TrieSlice::_smallDiff0( Source& src, Node n0, Position p, Depth d, Length l, Sequence s, Length minim ) {
	PositionLength pl0 = source.at( n0 );
	Position p0 = position( pl0 );
	Length l0   = length( pl0 );

	Depth dd = 0;
	Length diff = 0;

	for( dd = 0 ; dd < min( l0, l ) ; dd ++ ) {
		if( src.symbol( p0+dd ) & src.symbol( p+dd )) {
			continue;
		}

		if( ++diff == 1 ) {
// 	HERE: first difference found
// 	split if below minimum depth
			if(( dd > 0 ) && ( d+dd ) >= minim ) {

				Node n1 = splitNode( n0, src.symbol( p0+dd ), dd );
				eraseOccurrencesUnlessOwn( n0, s );

				_smallDiff1( src, n1, p+dd, d+dd, l-dd, s, minim );
				return;
			}

			continue;
		}

// 	HERE: second difference found
		assert( diff > 1 );
		assert( dd > 0 );

		if(( d+dd ) >= minim ) {

			splitNode( n0, src.symbol( p0+dd ), dd );
			eraseOccurrencesUnlessOwn( n0, s );
		}

		return;
	}

	assert( diff <= 1 );

	if( diff == 0 ) {
// 	HERE: no differences found
		if( l <= l0 ) {
			return;
		}

		for( SymbolNode sn: children( n0 )) {
			_smallDiff0( src, node( sn ), p+dd, d+dd, l-dd, s, minim );
		}
		return;
	}

	assert( diff == 1 );

	if( l < l0 ) {
		assert(( d+l ) >= minim );

		splitNode( n0, src.symbol( d+l ), l );
		eraseOccurrencesUnlessOwn( n0, s );
		return;
	}

	eraseOccurrencesUnlessOwn( n0, s );

	if( l > l0 ) {
		assert( dd == l0 );

		for( SymbolNode sn: children( n0 )) {
			_smallDiff1( src, node( sn ), p+dd, d+dd, l-dd, s, minim );
		}
	}
}

void TrieSlice::_confirm( Trie& trie, const Source& src, const string& s, Node n0, Sequence re, Position p, Depth d, Length l, Length minim ) {
	Position p0 = position( source.at( n0 ));
	Length l0   = length( source.at( n0 ));

	Depth dd;
	bool mismatch = false;

	bool has_cluster = cluster.count( n0 );

	for( dd = 0 ; dd < min( l0, l ) ; dd++ ){
// 	WARNING: matching between ambiguous symbols ("&" and not "=")
		if( src.symbol( p0+dd ) & s.at( p+dd )){
			continue;
		}

		mismatch = true;
		break;
	}

	if( !mismatch && ( l0 == l )) {
// (1)  n0 |-----------|E [s0,...]     =>   n0 |-----------|E [s0,...,s]
//      s  |-----------|

// 	NOTE: since the subsequence s, p, d, l has at least "minimum" length,
// 	it is safe to add s to the occurrences

		if( !has_cluster ) return; // no cluster == no occurrences
		if( src.commonSpecies( re, cluster.at( n0 ))) return; // keep this oligo; one of its species matches the species of the reference sequence
		cluster.erase( n0 );
		return;
	}

	if( !mismatch && ( l0 < l )) {
// (2)  n0 |-----------|E [s0,...]     =>   n0 |-----------|E [s0,...] & s?
//      s  |-----------====|                                +[mark matching children]

		Symbol sy = s.at( p+dd );

		if(( d+dd ) >= minim ) {
			if( !has_cluster ) return; // no cluster == no occurrences
			if( src.commonSpecies( re, cluster.at( n0 ))) return; // keep this oligo; one of its species matches the species of the reference sequence
			cluster.erase( n0 );
		}

		for( auto& c: children( n0 )) {
			if( sy & symbol( c )) {
				_confirm( trie, src, s, node( c ), re, p+dd, d+dd, Length( l-dd ), minim );
			}
		}

		return;
	}

	if( !mismatch ) {
// (3)  n0 |-------====|E [s0,...]     =>   n0 |-------| 
//      s  |-------|                                    +|====|E [s0,...]

		if( !has_cluster ) return; // no cluster == no occurrences
		if( src.commonSpecies( re, cluster.at( n0 ))) return; // keep this oligo; one of its species matches the species of the reference sequence

		Symbol sy = src.symbol( p0+dd );
		vector<SymbolNode> c0 = children( n0 );

		for( auto& c: c0 ) {
			children_hash.erase( symbolNode( symbol( c ), n0 ));
		}
		children_list.erase( n0 );

// 	insert n1
		Node n = newChild( n0, sy, p0+dd, l0-dd );
		for( auto& c: c0 ) {
			children_hash.emplace( symbolNode( symbol( c ), n ), node( c ));
			children_list.emplace( n, c );
		}

// 	shrink n0
		source.at( n0 ) = positionLength( p0, dd );

// 		occurrences of n
		cluster.emplace( n, cluster.at( n0 )); // copy the occurrences

		return;
	}

// (4)  n0 |-------====|E [s0,...]     =>   n0 |-------|   (if length>min) [s0,...,s] else []
//      s  |-------~~~|                              n1 +|====|E [s0,...]

	if(( d+dd ) < minim ) {
		return;
	}

	if( !has_cluster ) return; // no cluster == no occurrences
	if( src.commonSpecies( re, cluster.at( n0 ))) return; // keep this oligo; one of its species matches the species of the reference sequence

// 	only add occurrences to the nodes below the "minimum" depth
	Symbol sy = src.symbol( p0+dd );
	vector<SymbolNode> c0 = children( n0 );

	for( auto& c: c0 ) {
		children_hash.erase( symbolNode( symbol( c ), n0 ));
	}
	children_list.erase( n0 );

// 	insert n1
	Node n = newChild( n0, sy, p0+dd, l0-dd );
	for( auto& c: c0 ) {
		children_hash.emplace( symbolNode( symbol( c ), n ), node( c ));
		children_list.emplace( n, c );
	}

// 	shrink n0
	source.at( n0 ) = positionLength( p0, dd );

// 		occurrences of n
	cluster.emplace( n, cluster.at( n0 ));
	cluster.erase( n0 );

	return;
}

// ====================================================================================================================
// TEST
// ====================================================================================================================
/**
 * Print the portion of TrieSlice below a node using bracket notation, recursively
 * - children in round-brackets "()"
 * - occurrences in squiggly brackets "{}"
 */
// TODO: TrieSlice::show
void TrieSlice::show( const Source& src, Node n, Depth d ) const
{
// 	if( n ) {
// // 		print the contents of a node that is not the TrieSlice root
// 		cout << src.subsequence( position( source.at( n )), Length( length( source.at( n ))));
// 	}
// 
// // 	print occurrences
// // 	NOTE: no worries about printing occurrences in the root before its contents; the root has no occurrences
// 	auto o = occurrences.equal_range( n );
// 	if( o.first != o.second ) {
// 		cout << "{";
// 		bool bo = true;
// 		for( auto& oo = o.first ; oo != o.second ; oo++ ) {
// 			if(bo){bo = false;} else {cout << ",";}
// 			cout << oo->second;
// 		}
// 		cout << "}";
// 	}
// 
// 	auto ch = children( n );
// 	if( ch.empty()) {
// 		return;
// 	}
// 
// 	if( !n ) {
// // 	print the contents of the TrieSlice root from the portion accessible from its first child
// 		cout << src.subsequence( position( source.at( node( ch.front()))) - d, d );
// 	}
// 
// // 	print children, recursively
// 	cout << "(";
// 	bool bc = true;
// 	for( auto& c:ch ) {
// 		if(bc){bc = false;} else {cout << " ";}
// 		show( t, node( c ), Depth( d+length( source.at( n ))));
// 	}
// 	cout << ")";
}

/**
 * Print all paths to terminal nodes, "|"-delimited, recursively
 */
// TODO TrieSlice::print
void TrieSlice::print( const Source& src, Node n, Depth d, string root ) const {
}

/**
 * Recursively measure a TrieSlice
 * 
 * \param n0 start at this node (recursive parameter)
 * \param d depth of the node n0 (recursive parameter)
 * 
 * \param t trie
 * 
 * \param n accumulator for the total number of nodes in the TrieSlice
 * \param l accumulator for the total length of all nodes in the TrieSlice
 * \param o accumulator for the total number of occurrences in the TrieSlice
 * 
 * \param ld accumulator for the distribution of lengths in the TrieSlice
 * \param od accumulator for the distribution of occurrences in the TrieSlice
 * \param dd accumulator for the distribution of depths of nodes in the TrieSlice
 */
void TrieSlice:: measure(
	Node n0, Depth d,
	int& nodes, int& leaves, int& l, int& o,
	Distribution& depth_distribution, Distribution& length_distribution, Distribution& occurrence_distribution
) const {
	nodes++;
	if( children_hash.find( n0 ) == children_hash.end()) {
		leaves++;
	}

	l += length( source.at( n0 ));
	o += occurrences.count( n0 );

	depth_distribution[d]++; 
	length_distribution[length( source.at( n0 ) )]++;
	occurrence_distribution[occurrences.count( n0 )]++;

	for( auto& c: children( n0 )) {
		measure(
			node( c ), Depth( d+length( source.at( n0 ))),
			nodes, leaves, l, o,
			depth_distribution, length_distribution, occurrence_distribution );
	}
}

/*
 * Basic, recursive health check of the children of a node
 * 
 * TODO: TrieSlice::verify catch non-discoverable nodes (LOW)
 */
void TrieSlice::verify( const Source& src, Node n0 ) const {
	Node n;

	for( auto& c: children( n0 )) {
		if( !childAt( n0, symbol( c ), n )) {
			error("1");
			error( string{"child in list ('"}+nu2asc.at( symbol( c ))+"') but not in hash" );
		}

		if( n != node( c )) {
			error("2");
			error( string{"different child found in hash ('"}+nu2asc.at( symbol( c ))+"') than in list" );
		}

		if( symbol( c ) != src.symbol( position( source.at( node( c ))))) {
			error(
				string{"list ('"}
				+nu2asc.at( symbol( c ))
				+"') vs. source ('"
				+nu2asc.at(src.symbol( position( source.at( node( c )))))
				+"')"
			);
		}

		verify( src, node( c ));
	}
}

unordered_multimap<string,Sequence>& TrieSlice::find( const Source& so, const string& s, unordered_multimap<string,Sequence>& a ) const {
	return find( so, s, 0, fixed_depth, a );
}

unordered_multimap<string,Sequence>& TrieSlice::find( const Source& so, const string& s, Node n0, Depth d, unordered_multimap<string,Sequence>& a ) const {
	PositionLength pl0 = source.at( n0 );
	Position p0 = position( pl0 );
	Length l0 = length( pl0 );

	Length l = s.size();

	for( Depth dd=0 ; dd<l0 ; dd++ ) {
		if(( d+dd ) >= l ) {
			for( Sequence s: occ( n0 )){
				a.emplace( so.printableSubsequence( p0-d, l ), s );
			}
			return a;
		}

		if(!( so.symbol( p0+dd ) & s.at( d+dd ))) {
			return a;
		}
	}

	assert(( d+l0 ) <= l );

	if(( d+l0 ) == l ) {
		for( Sequence s: occ( n0 )){
			a.emplace( so.printableSubsequence( p0-d, l ), s );
		}
		return a;
	}

	for( const SymbolNode sn: children( n0 )) {
		find( so, s, node( sn ), d+l0, a );
	}

	return a;
}

Cluster TrieSlice::getClusterId( Source& src, Position p, Length l ) const { 
	assert( l > fixed_depth );

	return _getClusterId( src, 0, p+fixed_depth, fixed_depth, l-fixed_depth );
}

Cluster TrieSlice::_getClusterId( Source& src, Node n0, Position p, Depth d, Length l ) const {
	assert( l > 0 );

	Length l0   = length( source.at( n0 ));

	if( l <= l0 ) {
		assert( cluster.count( n0 ) > 0 );

		return cluster.at( n0 );
	}

	Node n;
	bool b = childAt( n0, src.symbol( p+l0 ), n );

	assert( b );

	return _getClusterId( src, n, p+l0, d+l0, l-l0 );
}

const bool TrieSlice::getCluster( const Trie& trie, const string& s, Position p, Length l, Cluster& clu ) const {
	assert( l > fixed_depth );

	return _getCluster( trie, 0, fixed_depth, s, p+fixed_depth, l-fixed_depth, clu );
}

const bool TrieSlice::_getCluster( const Trie& trie, Node n0, Depth d, const string& s, Position p, Length l, Cluster& clu ) const {
	assert( l > 0 );

	PositionLength pl0 = source.at( n0 );
	Position p0 = position( pl0 );
	Length l0 = length( pl0 );

	for( Depth dd = 0 ; dd < min( l0, l ) ; dd++ ){
		if( trie.source.symbol( p0+dd ) == s.at( p+dd )){
			continue;
		}

		return false; // if a mismatch is encountered, the target sequence is not found in the TrieSlice
	}

	if( l <= l0 ) { // the target sequence ends before the current node
		assert( cluster.count( n0 ) > 0 ); // WARNING: make sure that the length matches --oligo-size

		clu = cluster.at( n0 );
		return true;
	}

// 	here, the target sequence continues beyond the current node
	assert( l > l0 );

	Node n;
	bool b = childAt( n0, s.at( p+l0 ), n );

	if( !b ) {
		return false; // cannot continue in any of the children of the current node
	}

	return _getCluster( trie, n, d+l0, s, p+l0, l-l0, clu );
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
