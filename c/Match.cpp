#define __Match_cpp__

#include "Match.h"

void Match::onFragment( const string& na, const string& fn, const Cover<Position>& amb, const bool rc ){
	match_sequences.emplace_back( na, amb );

	if( match_sequences.size() >= match_buffer_size ) processFragments();
};

void Match::onFinish() {
	this->ParserFasta::onFinish();
	processFragments();
}

/**
 * Process a number of fragments read by the parser.
 * 
 * This can happen when match_buffer_size fragments have been read or at EOF input
 */
void Match::processFragments() {
	const ios_base::fmtflags floatflag = out.setf( ios::fixed, ios::floatfield ); // remember + set float fixed
	const unsigned int prec = out.precision( 1 );  // remember + set precision to 99.1

	spin( threads, *this, &Match::sequenceLoop ); // start a number of worker threads and wait until they finish

// restore float and precision
	out.setf( floatflag, ios::floatfield );
	out.precision( prec );

	match_sequences.clear();
	content.clear(); // forget/reuse the storage for the current reference sequence
}

/**
 * Worker function for processing sequences in a loop
 * 
 * Can be executed and synchronized between multiple threads
 */
void Match::sequenceLoop( bool& first ) {
	static mutex lock; // initialization lock

	static atomic<unsigned int> i( 0 ); // progress counter
	const unsigned int n = match_sequences.size();

	Alignment al;

	lock.lock();
	if( first ) { // initialization
		first = false;
		i = 0;
	}
	lock.unlock();

	while( true ) { // sequence loop
		const unsigned int ii = i++; // WARNING: will go over n in multiple threads
		if( ii >= n ) break;
		const auto& pa = match_sequences.at( ii );
		onSequence( pa.first, pa.second, al ); // call the worker for each sequence
	}
}

/**
 * Callback for processing a sequence from the input
 * 
 * Writes to output one line associated with the sequence
 * 
 * \param na FASTA name of the sequence
 * \param amb ambiguous Cover of the sequence
 * 
 * \see Cover
 */
void Match::onSequence( const string& na, const Cover<Position>& amb, Alignment& al ) {
	const LLength target_sequence_length = amb.range().size();

	const double min_cluster_area_ratio = 0.75; // ignore target sequences with less area explained by clusters
	assert( min_cluster_area_ratio > 0 );

	const LLength min_cluster_area = min_cluster_area_ratio * target_sequence_length;

	const double size_factor = 2.0; // multiplier for the signature size vs. minimum target sequence length

	set<Sequence> max_set;

	if( target_sequence_length < size_factor * le ) { // ignore short sequences
		print( out,
			na, "",
			0.0, target_sequence_length,
			target_sequence_length, 0, max_set.size()
		);

		return;
	}

	if( target_sequence_length >= Alignment::max_length_smallest_sequence )
		error(
			"target sequence too long (", na, "): ",
			target_sequence_length, " >= ", Alignment::max_length_smallest_sequence
		);

// 	(1) Calculate list of clusters matching the sequence
	set<Cluster> set_clusters;
	deque<pair<Position,Cluster>> po_cluster;

	Position last_position( 0 ); // last position covered by any cluster
	LLength cluster_area( 0 );   // total area covered by any cluster

	for( const Range<Position>& ra: -amb ) { // O(1)
		for( Position p = ra.lo() ; Length l = ra.cover( p, le, le ) ; p++ ) { // O(s) for each cluster matching the target sequence
			Cluster clu( Cluster_invalid );
			if( !trie.getCluster( content, p, l, clu )) continue; // next if there is no cluster at the current position

			assert( clu != Cluster_invalid );

			if( p <= last_position ) cluster_area -= ( last_position - p );

			last_position = p + l;
			cluster_area += l;

			if( !contains( set_clusters, clu )) max_set += trie.source.clusters.at( clu );

			set_clusters.emplace( clu );
			po_cluster.emplace_back( p, clu );

		}
	}

	if( cluster_area < min_cluster_area ) { // skip target sequences with too many positions unexplained by any cluster
		print( out,
			na, "",
			100. * cluster_area / target_sequence_length, target_sequence_length,
			target_sequence_length, 0
		);

		return;
	}

// 	(2) build the minimum set of sequences
	struct ClusterSizeCompare{
		bool operator() ( const Cluster& c1, const Cluster& c2 ) {
			const size_t s1 = clusters.at( c1 ).size();
			const size_t s2 = clusters.at( c2 ).size();

			if( s1 < s2 ) return true;
			if( s1 > s2 ) return false;

			return ( c1 < c2 );
		}
		ClusterSizeCompare( const One2One<set<Sequence>,Cluster>& clus ) : clusters( clus ) {};
		const One2One<set<Sequence>,Cluster>& clusters;
	} clusterSizeCompare( trie.source.clusters );

// 	set of clusters sorted by size
	set<Cluster,ClusterSizeCompare> sorted_clusters( set_clusters.begin(), set_clusters.end(), clusterSizeCompare );
	set<Sequence> min_set;

	for( Cluster cl1 : sorted_clusters ) {
		set<Sequence> s1 = trie.source.clusters.at( cl1 );
		if(!( min_set & s1 ).empty()) continue;

		last_position = 0;
		cluster_area = 0;

		for( const auto& e2pocl: po_cluster ) {
			const Position& p = e2pocl.first;
			const Cluster& cl2 = e2pocl.second;

			set<Sequence> s2 = trie.source.clusters.at( cl2 );
			s2 &= s1;

			if( s2.empty()) continue; // the current cluster has no overlap with the minimum cluster

			s1.swap( s2 );

			if( p <= last_position ) cluster_area -= ( last_position - p );

			last_position = p + le;
			cluster_area += le;
		}

		if( cluster_area < min_cluster_area ) continue;

		min_set += s1;
	}

	if( !min_set.size()) { // skip target sequences with too many positions unexplained by any cluster
		print( out,
			na, "",
			100. * min_cluster_area_ratio, target_sequence_length,
			target_sequence_length, 0, max_set.size()
		);

		return;
	}


// 	(3) Align the target sequence to [all fragments of] sequences of the minumum set
	deque<pair<Sequence,pair<double,LLength>>> aligned_sequences; // sequence by match percentage and overlap length

	for( const Sequence& se: min_set ) {
		for( auto it2fr = trie.source.instance_fragments.from.equal_range( se ) ; it2fr.first != it2fr.second ; ++it2fr.first ) {
			const Fragment& fr = trie.source.fragments.at( it2fr.first->second );
			const auto a2 = al.align(
				content, amb.range().lo(), amb.range().size(),
				trie.source.getSource(), fr.getRange().lo(), fr.getRange().size()
			);

			assert( a2.L > 0 );
			aligned_sequences.emplace_back( se, make_pair( 100.0 * a2.M / a2.L, a2.L ));
		}
	}

	assert( aligned_sequences.size() > 0 );

// 	(4) pick source sequences with highest match percentage
	struct { // reverse lexicographic order on match percentage
		bool operator() ( const pair<Sequence,pair<double,LLength>> a, const pair<Sequence,pair<double,LLength>> b ) const { 
			return ( a.second.first < b.second.first );
		};
	} sequenceAlignmentCompare;

// 	maximum match percentage after alignment
	const double max_score = max_element( aligned_sequences.begin(), aligned_sequences.end(), sequenceAlignmentCompare )->second.first;

	for( const auto& as: aligned_sequences ) {
		if( as.second.first != max_score ) continue;
		
		print( out,
			na, trie.source.instances.at( as.first ),
			as.second.first, as.second.second,
			target_sequence_length, min_set.size(), max_set.size()
		);
	}
}
/**
 * Display the result of a match against the target sequence on one line
 */
ostream& Match::print(
	ostream& out,

	const string& target_sequence_name,
	const string& source_sequence_name,

	const double match_percentage,
	const LLength overlap,

	const LLength target_sequence_length,
	const LLength min_set_size,
	const LLength max_set_size
) {
	static mutex lock;

	lock.lock();
	out
		<< target_sequence_name // name of the target sequence
		<< '\t' << ( !source_sequence_name.size() ? "-" : source_sequence_name.c_str()) // name of the best matching source sequence

		<< '\t' << (( !source_sequence_name.size() && ( match_percentage > 0 )) ? "<":"" ) << match_percentage << '%' // matching percentage
		<< '\t' << overlap

		<< '\t' << target_sequence_length
		<< '\t' << min_set_size
		<< '\t' << max_set_size

		<< endl;
	lock.unlock();

	return out;
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
