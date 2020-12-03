#ifndef __util_h__
#define __util_h__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include <cassert>

#include <sys/types.h>
#include <sys/stat.h>
#include <wordexp.h>

#include <thread>

#include "Error.h"

using namespace std;

/**
 * Lowers all characters of a string \param s in place
 * 
 * \returns the lowered string itself
 * 
 * WARNING: assumes single-byte-per-character; does NOT support multibyte encodings (Unicode)
 */
inline string& lower( string& s )
{
	transform( s.begin(), s.end(), s.begin(), ::tolower ); // single character encoding
	return s;
}

/**
 * \returns a copy of string s reversed in place
 * 
 * WARNING: assumes single-byte-per-character; does NOT support multibyte encodings (Unicode)
 */
inline string operator- ( const string& s ) {
	string _s = s;
	reverse( _s.begin(), _s.end());
	return _s;
}

/**
 * True if string s starts with prefix
 * 
 * \param suffix the remaining part after stripping away the prefix
 */
inline bool startsWith( const string& s, const string& prefix, string& suffix )
{
	suffix = "";

	if( s.find( prefix ) != 0 )
		return false;

	suffix = s.substr( prefix.length());
	return true;
}

/**
 * String to integer
 */
inline bool convert( const string& source, int& destination )
{
	stringstream ss;

	ss << source;
	ss >> destination;

	return !ss.fail() && !ss.bad();
}

/**
 * String to bool
 */
inline bool convert( string source, bool& destination )
{
	lower( source );

	if(( source == "yes" ) || ( source == "y" ) || ( source == "true" ) || ( source == "1" )){
		destination = true;
		return true;
	}

	if(( source == "no" ) || ( source == "n" ) || ( source == "false" ) || ( source == "0" )){
		destination = false;
		return true;
	}

	return false;
}

/**
 * String to double
 */
inline bool convert( const string& source, double& destination )
{
	stringstream ss;

	ss << source;
	ss >> destination;

	return !ss.fail() && !ss.bad();
}

/**
 * Split string source into parts, as separated by separator
 */
inline void split( const string& source, const string& separator, string& r1, string& r2 )
{
	size_t f = source.find( separator );

	switch( f ){
		case 0:
		case string::npos:
			r1 = source;
			r2 = "";
			return;
		default:
			r1 = source.substr( 0, f );
			r2 = source.substr( f+separator.length());
			return;
	}
}

/**
 * Test whether file can be open for reading
 * 
 * NOTE: reading from the file may still fail
 */
inline bool canOpenFileRead( const string& file )
{
	struct stat sb;

	if( stat( file.c_str(),&sb ) != 0 ) return false; // file does not exist
	if( !S_ISREG( sb.st_mode )) return false;         // not a regular file

	ifstream f( file );          // open for reading; will be closed when going out of scope
	if( !f.good()) return false; // file cannot be read

	return true;
}

/**
 * Test whether file can be open for reading
 * 
 * NOTE: reading from the file may still fail
 */
inline bool canOpenFileWrite( const string& fn )
{
	struct stat sb;
	const int sf = stat( fn.c_str(), &sb );

	if( sf == 0 ) { // file already exists
		if( !S_ISREG( sb.st_mode )) return false; // not a regular file; bail!

		ofstream f( fn, ios_base::app ); // do not destroy existing content
		if( !f.good()) return false;     // file exists; but cannot be written

		return true;
	}

	assert( sf == -1 ); // invalid return from stat (unistd.h)

	ofstream f( fn );   // try to open file for writing
	bool ok = f.good();
	f.close();
	remove( fn.c_str()); // cleanup; remove file

	return ok;
}

/**
 * Add file names (including wildcards) from an expression (source) to a list
 */
inline bool addFiles( const string& source, vector<string>& list )
{
	wordexp_t p;
	if( wordexp( source.c_str(), &p, 0 )){
		wordfree( &p );
		error( "cannot expand expression (", source, ")" );
	}

	if( !p.we_wordc ){
		wordfree( &p );
		return false;
	}

	for( unsigned int i=0 ; i< p.we_wordc ; i++ ){
		if( !canOpenFileRead( p.we_wordv[i])){
			error( "cannot open file (", p.we_wordv[i], ")" );
		}

		list.emplace_back( p.we_wordv[i] );
	}

	wordfree( &p );
	return true;
}

/**
 * Split a string \param s into chunks separated by \param sep
 */
inline vector<string> split( const string& s, const string& sep ) {
	vector<string> r;
	size_t p1 = 0;
	size_t p2 = 0;

	while( true ) {
		p2 = s.find( sep, p1 );

		if( p2 == string::npos ) {
			r.emplace_back( s.substr( p1, string::npos ));
			return r;
		}

		r.emplace_back( s.substr( p1, p2-p1 ));
		p1 = p2+1;
	}
}

/**
 * Set intersection
 */
template<typename T> inline set<T> operator& ( const set<T>& s1, const set<T>& s2 )
{
	set<T> result;
	set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), inserter( result, result.begin()));
	return result;
}

/**
 * Set intersection
 */
template<typename T> inline set<T>& operator&= ( set<T>& s1, const set<T>& s2 )
{
	set<T> result;
	set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), inserter( result, result.begin()));
	s1.swap( result );
	return s1;
}

/**
 * Set union
 */
template<typename T> set<T> operator+ ( const set<T>& s1, const set<T>& s2 )
{
	set<T> result;
	set_union( s1.begin(), s1.end(), s2.begin(), s2.end(), inserter( result, result.begin()));
	return result;
}

/**
 * Set union
 */
template<typename T> inline set<T>& operator+= ( set<T>& s1, const set<T>& s2 )
{
	set<T> result;
	set_union( s1.begin(), s1.end(), s2.begin(), s2.end(), inserter( result, result.begin()));
	s1.swap( result );
	return s1;
}

/**
 * Set difference
 */
template<typename T> set<T> operator- ( const set<T>& s1, const set<T>& s2 )
{
	set<T> result;
	set_difference( s1.begin(), s1.end(), s2.begin(), s2.end(), inserter( result, result.begin()));
	return result;
}

/**
 * Set union
 */
template<typename T> inline set<T>& operator-= ( set<T>& s1, const set<T>& s2 )
{
	set<T> result;
	set_difference( s1.begin(), s1.end(), s2.begin(), s2.end(), inserter( result, result.begin()));
	s1.swap( result );
	return s1;
}

/**
 * Container "element containment" operator
 */
template<typename C, typename T> inline bool contains( const C& s, const T& e )
{
	return s.find( e ) != s.end();
}

/**
 * Set "less" operator
 */
template<typename T> inline bool setLess( const set<T>& a, const set<T>& b)
{
	return ( a.size() < b.size());
}

/**
 * Set "inclusion" operator
 */
template<typename T> inline bool contains( const set<T>& a, const set<T>& b )
{
	return includes( a.begin(), a.end(), b.begin(), b.end());
}

/**
 * Start and wait for the finishing of \param t threads, using function \param f of object \param o
 * and passing to each the parameter list \param args
 * 
 * Each thread function is responsible with managing its own workload (e.g. using static static variables)
 * and not interfering with other threds (e.g. using mutexes)
 * 
 * Each thread function expects a boolean parameter ("first") that indicates whether this is the first
 * thread in the pool (to support "unique initialization
 * 
 * WARNING: derived functors cannot be passed as parameters
 */
template<class C, typename... Args>
void spin( unsigned int t, C& o, void ( C::* f )( bool&, Args&... ), Args&... args ) {
	vector<thread> v;  // thread pool
	static bool first; // shared flag indicating whether this is the first thread executing

	first = true;

	for( unsigned int i=0 ; i<t ; i++ ){ // start all threads
// 	equivalent to thread t( f, &o, ref( first ), args... );
		v.emplace_back( f, &o, ref( first ), args... );
// 	NOTE: "first" must be passed by reference, so that "unique initialization" can be shared by all spinned threads
	}

	for( thread& th: v ){ // wait for all threads to finish
		th.join();
	}
};


// TEST
template<typename T> ostream& operator<< ( ostream& o, const set<T>& s )
{
	o << '{';
	for( const T& e : s ) {
		o << e << ' ';
	}
	o << '}';
	return o;
}

#ifndef __util_cpp__
extern ostream onull;
#endif

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
