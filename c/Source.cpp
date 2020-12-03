#include "Source.h"

/**
 * Determine if a fragment has mode ambiguities than max_ambiguities
 */
bool Source::filterAmbiguous( const string& na, const string& fn, const Cover<Position>& amb )
{
	if(( max_ambiguities > 0 ) && ( amb.length() > max_ambiguities )) {
		excluded.emplace_back( na + '\t' + fn );
		return true;
	}

	return false;
}

/**
 * Count ambiguities "crowded" inside a window of size "maxim"
 * if the "number of crowded ambiguities" > max_ambiguities:
 * the fragment will not be added to the Source and add its name will be added to excluded.fasta
 */
bool Source::filterAmbiguousCrowded( const string& na, const string& fn, const Cover<Position>& amb )
{
	if(( max_crowded_ambiguities > 0 ) && ( amb.window( maxim ) > max_ambiguities )) {
		excluded.emplace_back( na + '\t' + fn );
		return true;
	}

	return false;
}

void Source::filterOutgroup( const string& outgroup_file_name )
{
	if( !outgroup_file_name.size())
		return;

	fstream f;
	string item;

	try{
		f.open( outgroup_file_name, fstream::in );
	}catch( ios_base::failure f ){
		error( "cannot open outgroup file (", outgroup_file_name, "): ", f.what());
	}

	while( true ){
		f >> item;
		if( f.eof() )
			break;

		bool found = false;

// 	Lookup item in the list of sequences
		for( const auto& e2in: instances.from ){
			if( e2in.first.find( item ) == string::npos )
				continue;

			found = true;

			targets.erase( e2in.first );
		}

		if( !found ){
			error( "cannot find sequence specified in outgroup file (", item, ")" );
		}
	}
}

void Source::readIsolationList( const string& isolation_file_name )
{
	if( !isolation_file_name.size())
		return;

	fstream f;
	string item;

	try{
		f.open( isolation_file_name, fstream::in );
	}catch( ios_base::failure f ){
		error( "cannot open isolation file (", isolation_file_name, "): ",  f.what());
	}

	while( true ){
		f >> item;
		if( f.eof() )
			break;

		bool found = false;

// 	Lookup item in the list of sequences
		for( const auto& s: instances.from ){
			if( s.first.find( item ) == string::npos )
				continue;

			found = true;
			isolation.insert( s.first );
		}

		if( !found ){
			error( "cannot find sequence specified in isolation file (", item, ")" );
		}
	}

	for( const auto& e2in: instances.from ) {
		if( !isolation.count( e2in.first )) {
			targets.erase( e2in.first );
		}
	}
}

void Source::filterMelting( const unsigned int threads, const double max_melting, const double strand_concentration, const double salt_concentration ) {
	if( &fold_output != &onull ) fold_output << fixed;
	spin( threads, *this, &Source::_filterMelting, max_melting, strand_concentration, salt_concentration );
}

/**
 * Remove occurrences of subsequences that have melting temperatures higher than \param max_melting
 */
void Source::_filterMelting( bool& first, const double& max_melting, const double& strand_concentration, const double& salt_concentration ) {
	const auto e = fragments.from.end();
	static One2One<TypeFragment,Fragment>::type_from::iterator i;
	static mutex lock;
	static double t0 = max_melting + Thermo::K; // convert from Celsius to Kelvin

	static Thermo th( t0, strand_concentration, salt_concentration );

	lock.lock();

	while( true ) {
		if( first ) {
			i = fragments.from.begin();
			first = false;
		}

		if( i == e ) break;
		auto ii = i;

		lock.unlock();
		for( const Range<Position>& r: ii->second.getAmbigCompl()) {
			if( r.size() < minim ) continue; // skip ranges that are too small

// 	WARNING: allocating stack storage for the Fold (like so: "Fold h") fails on some systems (clusters)
// 	Possible explanation: stack overflow for on stack storage
// 	Solution: Allocate the Fold on heap storage (new Fold)
			Fold* h = new Fold( content, r.lo(), r.size(), max_length_at, minim, min( r.size(), Position( maxim )), th, fold_output );
			h->fold();
			delete h; // make sure to delete the Fold !
		}
		lock.lock();

		i++;
	}
	lock.unlock();

	return;
}

/**
 * Print all fragments from all sequences (EXPERIMENTAL)
 * 
 * \see Source::showSequence
 */
void Source::show( ostream& out ) const {
	for( const auto& e: instances.to ) {
		showSequence( out, e.first );
	}
}

/**
 * Print all fragments from a sequence (EXPERIMENTAL)
 * 
 * For each fragment, print one line for each of the following: 
 *  - sequence name
 *  - actual fragment
 *  - positions of ambiguous bases in the fragment
 *  - the actual ambiguous bases
 */
void Source::showSequence( ostream& out, Sequence s ) const {
	for( auto i = instance_fragments.from.equal_range( s ) ; i.first != i.second ; ++( i.first )) { // for each fragment in the sequence
		const Fragment& fr = fragments.at( i.first->second );

		out << instances.at( s ) << endl; // print sequence name
		out << printableSubsequence( fr.getRange()) << endl; // print actual fragment

		out << fr.getAmbig() << endl; // print the positions of ambiguous bases in the fragment

		for( const Range<Position>& r: fr.getAmbig()) {
			out << printableSubsequence( r );  // print all ambiguous bases
		}

		out << endl;
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
