#ifndef __Source_h__
#define __Source_h__

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <stack>
#include <map>
#include <unordered_map>

#include <atomic>

#include "Tree.h"

#include "Types.h"
#include "Cover.h"
#include "Relation.h"
#include "Fragment.h"
#include "Thermo.h"

#include "Fold.h"

#include "ParserFasta.h"

#include "util.h"

using namespace std;

extern int  newick_parse( Source& );
extern void newick_restart( FILE* );
extern void newickClear();

extern int  tax_parse( Source& );
extern void tax_restart( FILE* );

/**
 * Source database of sequences
 */
class Source : public ParserFasta
{
	friend class Reference;
private:
	Length minim;
	Length maxim;

	Position max_ambiguities;
	Position max_crowded_ambiguities;
	Length max_homolo;

	ostream& fold_output;
public:
	One2One<string,Sequence> instances;
	One2Many<Sequence,TypeFragment> instance_fragments;
	One2One<TypeFragment,Position> fragment_position;

	One2One<TypeFragment,Fragment> fragments; // fragment id to fragment body

	One2One<string,Sequence,unordered_map<string,Sequence>> reference; // lookup: reference-sequence-name -> reference-sequence-id

	One2One<string,Species,unordered_map<string,Species>> species; // list of species names from the taxonomy file (temporary)
	One2Many<Species,string,unordered_multimap<Species,string>> species_reference; // lookup: species-name -> reference-sequence-name (temporary)
	unordered_map<Sequence,Species> instance_species; // lookup: instance (many) --> species-name (one) (temporary)

	One2One<set<Sequence>,Cluster> clusters;      // clusters are groups of Sequences
	map<Cluster,set<Species>> cluster_species;    // 

	/**
	 * Map of "sets of sequences" (set of Sequence id's) to names for targets of oligo signature search
	 * NOTE: a set with just one sequence element represents oligo signatures for the sequence itself
	 */
	One2One<set<Sequence>,string> targets;

	/**
	 * Array of "melting lengths" corresponding to each position in the Source content
	 * 
	 * Elements of type Length (unsigned char) 
	 *  - 0 : subsequence cannot start at given position (end of the Cover)
	 *  - minim - 1 : minimum length subsequence starting at position is solid at melting temperature
	 *  - minim <= mele <= maxim : maximum length of an unmelted subsequence
	 */
	deque<Length> max_length_at;

	Source( Length m, Length M, Position ma, Position mca, Length mh, ostream& fo = onull, const bool rc = false ) :
		ParserFasta( rc ),
		minim( m ),
		maxim( M ),
		max_ambiguities( ma ),
		max_crowded_ambiguities( mca ),
		max_homolo( mh ),
		fold_output( fo ),
		has_tree( false )
		{};

	virtual void onFragment( const string& na, const string& fn, const Cover<Position>& amb, const bool rc = false ) {
		static Sequence seq = 0; // sequence id
		static TypeFragment fra = 0; // fragment id

		max_length_at.resize( max_length_at.size() + amb.range().size(), 0 ); // extend the array of lengths; default = 0

		for( const Range<Position>& r: -amb ) { // initialize array of lengths with "cover" values
			for( Position p = r.lo() ; ; p++ ) {
				const Length le = r.cover( p, minim, maxim );
				if( !le ) break;

				max_length_at.at( p ) = le;
			}
		}

		if( filterAmbiguous( na, fn, amb )) {
			return;
		}

		if( filterAmbiguousCrowded( na, fn, amb )) {
			return;
		}

		Sequence se = seq;

		int r = instances.emplace( na, seq );
		if( !r ) {
			targets.emplace( { seq }, na ); // add sequence to list of targets
			++seq;
		}

		switch( r ) {
			case 0:
				instance_fragments.emplace( se, fra );
				fragment_position.emplace( fra, amb.range().hi());
				fragments.emplace( fra++, Fragment{ fn, amb, maxim, rc });
				break;
			case 1:
				instance_fragments.emplace( instances.at( na ), fra );
				fragment_position.emplace( fra, amb.range().hi());
				fragments.emplace( fra++, Fragment{ fn, amb, maxim, rc });
				break;
			case 2:
			default:
				error( "A sequence with this id already exists: ", na, "/", seq );
		}
	};

	/**
	 * Write names of excluded fragments to a file
	 */
	void printExcluded( const string& file_name = "excluded.fasta" ) {
		if( excluded.empty()) {
			return;
		}

		fstream sink_excluded;

		try{
			sink_excluded.open( file_name, fstream::out );
		}catch( ios_base::failure f ){
			error( "cannot open ambiguities exclusion file (", file_name, "); ", f.what());
		}

		for( const string& e: excluded ) {
			sink_excluded << e << endl;
		}

		sink_excluded.close();
	};

	string last_newick_token;
	/**
	 * Interface with the Newick flex/bison parser
	 * 
	 * The parser will call the event listeners setTree and setLastTreeLabel
	 */
	void parseNewick( const string& path ) {
		assert( !has_tree );

		if( !path.size()) {
			return;
		}

		FILE* newick_in = fopen( path.c_str(), "rt" );

		if( !newick_in ) {
			error( "cannot open Newick tree file ( ", path, " )" );
		}

// 	connect flex to open input file
		newick_restart( newick_in );

// 	bison parser call
		if( newick_parse( *this )) {
// HERE: trap parsing errors

// 	NOTE: yyparse will return 0 when everything is ok; 1 on YYABORT; 2 on memory exhaustion
// 	reference: http://www.gnu.org/software/bison/manual/html_node/Parser-Function.html
// 	WARNING: yyparse may not return 0 even if yyerror was NOT called
			if( last_newick_token == "__no_token__" ) error(
				"cannot parse Newick tree file (", path, ")\n",
				"** failed before reading the first node label\n",
				"** Newick tree format reference: \n",
				"   http://evolution.genetics.washington.edu/phylip/newick_doc.html"
			);

			error(
				"cannot parse Newick tree file (", path, ")\n",
				"** last label read:", last_newick_token, "\n",
				"** Newick tree format reference: \n",
				"   http://evolution.genetics.washington.edu/phylip/newick_doc.html"
			);
		}

		fclose( newick_in );

		newickClear(); // clear temporary data structures

		tree.label( "Node" );

// 	Add groups to targets
		for( const auto& g: tree.getGroups()){
			if( g.second.size() < 2 )
// 	ignore groups with only one sequence
// 	NOTE: individual sequence are added above
				continue;

			set<Sequence> sg;

			for( const string& sn: g.second ){
				if( !instances.has( sn )) {
// 	NOTE: sequences that are found in the phylogeny tree, but not found in the data source are ignored
// 	      This is necessary when sequences are excluded because of ambiguities (--max-ambiguities or --max-crowded-ambiguities)
					continue;
				}

				sg.emplace( instances.at( sn ));
			}

			targets.emplace( sg, g.first );
		}
	};

	inline void setTree( Tree& t ) { tree = t; has_tree = true; };
	inline const Tree& getTree() { return tree; };

// parsing Taxononomy files
	void parseTaxonomy( const string& path ) {
// 	Extract Genus/species from names that look like: XX_999999_Genus_species...
		for( const auto& e2in: instances.from ) {
			auto v = split( e2in.first, "_" ); // split instance name
			if( v.size() < 4 ) continue; // name does not seem to be correctly encoded; ignore it!

			string  spe = lower( v.at( 2 )) + "_" + lower( v.at( 3 )); // canonical species name: "genus_species"
			static atomic<Species> _spe; // zero-initialized species counter

			if( !species.has( spe )) {
				species.emplace( spe, _spe++ );
			}

			instance_species.emplace( e2in.second, species.at( spe ));
		}

// 	Build set of species for each Cluster
		for( const auto& e2cl: clusters.from ) {
			set<Species> ss;
			for( Sequence se: e2cl.first ) {
				if( !instance_species.count( se )) continue; // ignore sequences without a species
				ss.emplace( instance_species.at( se ));
			}
			cluster_species.emplace( e2cl.second, ss );
		}

// 	Actual parsing starts here
		FILE* tax_in = fopen( path.c_str(), "rt" );

		if( !tax_in ) {
			error( "cannot open sequence file ( ", path, " )" );
		}

// 	connect flex to open input file
		tax_restart( tax_in );

// 	bison parser call
		if( tax_parse( *this )) {
// 	NOTE: yyparse will return 0 when everything is ok; 1 on YYABORT; 2 on memory exhaustion
// 	reference: http://www.gnu.org/software/bison/manual/html_node/Parser-Function.html
			error( "cannot parse sequence file ( ", path, " )" );
		}

		fclose( tax_in );

// 	Clean temporary structures
		species.from.clear(); species.from.rehash( species.from.size());
		species.to.clear();

		species_reference.from.clear(); species_reference.from.rehash( species_reference.from.size());
		species_reference.to.clear();

		instance_species.clear(); instance_species.rehash( instance_species.size());
	};

	void onTaxonomyEntry( const string& id, const string& sp ) {
		static atomic<Sequence> re;

		if( reference.has( id )) {
			error( "taxonomy file: multiple entries with the same name (", id, ")" );
		}

		reference.emplace( id, re++ );

		auto v = split( sp, "_" );
		if( v.size() < 2 ) return; // does not encode a species
		string spe = lower( v.at( 0 )) + "_" + lower( v.at( 1 ));

		if( species.has( spe )) {
			species_reference.emplace( species.at( spe ), id );
		}
	};

	bool filterAmbiguous( const string& na, const string& fn, const Cover<Position>& amb );
	bool filterAmbiguousCrowded( const string& na, const string& fn, const Cover<Position>& amb );

	void filterOutgroup( const string& outgroup_file_name );
	void readIsolationList( const string& isolation_file_name );

	void filterMelting( const unsigned int threads, const double max_melting, const double strand_concentration, const double salt_concentration );
	void _filterMelting( bool& first, const double& max_melting, const double& strand_concentration, const double& salt_concentration );

	inline const string& getSource() const { return content; };
	inline deque<Length>& getMaxLengthAt() { return max_length_at; };
	inline const string& getSequenceName( Sequence i ) const { return instances.at( i ); };

	inline const map<set<Sequence>,string>& getTargets() { return targets.from; };

	inline string printableSubsequence( Position p, Position l ) const {
		return convertNu2Asc( content.substr( p, l ));
	};
	inline string printableSubsequence( const Range<Position>& r ) const {
		return printableSubsequence( r.lo(), r.hi()-r.lo());
	};
	inline string printableSubsequence( const Cover<Position>& c ) const {
		return printableSubsequence( c.range());
	};

	inline Symbol symbol( Position p ) const {
		return content.at( p );
	};

	inline Position length() const {
		return content.size();
	};

	inline TypeFragment getFragmentAtPosition( const Position p ) {
		auto s = fragment_position.to.lower_bound( p );
		if( s == fragment_position.to.end()) {
			error( "cannot find any fragment at position", p, "in the sequence database" );
		}

		return s->second;
	};

	inline bool commonSpecies( Sequence re, Cluster cl ) const {
		return
			species_reference.has( reference.at( re ))
			&& cluster_species.at( cl ).count( species_reference.at( reference.at( re )));
	};

// 	TEST
	void show( ostream& out ) const;
	void showSequence( ostream& out, Sequence s ) const;

	void find( const vector<string>& v ) const {
		for( const string& s: v ) {
			find( s );
		}
	};
	void find( const string& s ) const {
		for( auto i: instances.to ) {
			findInSequence( i.first, convertAsc2Nu( s ));
		}
	};
	void findInSequence( Sequence se, const string& s ) const {
		for( auto f = instance_fragments.from.equal_range( se ) ; f.first != f.second ; ++( f.first )) {
			findInFragment( se, fragments.at( f.first->second ), s );
		}
	};
	void findInFragment( Sequence se, const Fragment& fr, const string& s ) const {
		assert( s.size() < 256 );
		Length l = s.size();

		Position lo = fr.getRange().lo();
		Position hi = fr.getRange().hi() - s.size();

		for( Position p = lo ; p < hi ; ++p ) {
			Length i=0;

			for( i = 0 ; i < l ; ++i ) {
				if( s.at( i ) & symbol( p+i )) {
					continue;
				}

				break;
			}

			if( i == l ) {
				cout << convertNu2Asc( s ) << '\t' << printableSubsequence( p, l ) << '\t' << instances.at( se ) << endl;
			}
		}
	};

private:
	/**
	 * Target names that should be the only ones included in the result of the oligo signature search
	 */
	set<string> isolation;
	/**
	 * FASTA names of sequences excluded from the analysis because of ambiguities
	 */
	vector<string> excluded;
	/**
	 * Whether there is a phylogeny tree associated with the sequence database
	 */
	bool has_tree;
	/**
	 * The actual phylogeny tree
	 * When there is no phylogeny tree, the tree will be empty
	 */
	Tree tree;
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
