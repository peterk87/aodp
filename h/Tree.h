#ifndef __Tree_h__
#define __Tree_h__

#include <string>
#include <vector>
#include <set>

#include "util.h"

using namespace std;

class Source;

/**
 * A phylogeny tree
 * Leaves have a name and a length
 * 
 * TODO: Tree::private members + get/set
 */
class Tree{
	public:
		Tree( const string& n = "", const vector<Tree>& c = vector<Tree>{}, const string& l = "" ):  len(l), name(n), children(c){};

		string& show( string& result, bool is_top = true ) const;
		vector<string>& showLineage( vector<string>& result ) const;
		vector<pair<string,set<string>>> getGroups() const;

		void label( string root_name = "Node" );
		void label( unsigned int &l, const string& root_name );

// 	NOTE: the bison parser needs access to the internal structure (private members) of the tree
		friend int newick_parse( Source& source );

		/**
		 * Print a Newick tree, recursively
		 * 
		 * WARNING: the final ';' needs to be added after this
		 */
		friend ostream& operator<< ( ostream& o, const Tree& t ) {
			o << t.children << t.name;

			if( t.len.size()){
				o << ':' << t.len;
			}
			return o;
		};

		/**
		 * Return a new tree, with \param suffix appended to each node whose name is in \param m
		 */
		Tree mark( const set<string>& m, const string& suffix = "*" ) const {
			Tree result = *this; // deep, recursive copy
			result._mark( m, suffix );
			return result;
		};

	private:
		/**
		 * Append \param suffix to each node whose name is in \param m (recursive worker function)
		 */
		void _mark( const set<string>& m, const string& suffix = "*" ) {
			if( m.count( name )) {
				name += suffix;
			}
			for( Tree& c: children ) {
				c._mark( m, suffix );
			}
		};

		friend ostream& operator<< ( ostream& o, const vector<Tree>& c ) {
			if( c.empty()) return o;

			o << '(';

			bool first = true;
			for( const Tree& child: c ) {
				if( first ) {
					first = false;
				}
				else {
					o << ',';
				}
				o << child;
			}

			o << ')';
			return o;
		};

		vector<pair<string,set<string>>>& getGroups( vector<pair<string,set<string>>>& result ) const;

		string len;
		string name;
		vector<Tree> children;
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
