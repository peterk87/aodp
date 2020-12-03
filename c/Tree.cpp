#include "Tree.h"

/**
 * Print a phylogeny tree, recursively
 */
string& Tree::show( string& result, bool is_top ) const
{
	if( children.empty()){
		result += name;
		if( len.size()){
			result += ":";
			result += len;
		}
		return result;
	}

	result += "(";

	bool first = true;

	for( const Tree& child: children ){
		if( first )
			first = false;
		else
			result += ",";

		child.show( result, false );
	}

	result += ")";
	result += name;
	if( len.size()){
		result += ":";
		result += len;
	}

	if( is_top ){
		result += ";";
	}

	return result;
};

/**
 * Calculate the lineage (list of ancestors for every leaf node) of a phylogeny tree
 * 
 * \return reference to the result (input parameter)
 */
vector<string>& Tree::showLineage( vector<string>& result ) const
{
	static vector<string> lineage;

	if( children.empty()){
		string element;

		for( string parent: lineage ){
			element += parent;
			element += ":";
		}

		element += name;

		result.push_back( element );
		return result;
	}

	lineage.push_back( name );
	for( const Tree& child: children ){
		child.showLineage( result );
	}
	lineage.pop_back();

	return result;
}

/**
 * Enumerate the groups (list of leaf nodes) in a phylogeny tree
 * 
 * \return reference to the result (input parameter)
 */
vector<pair<string,set<string>>>& Tree::getGroups( vector<pair<string,set<string>>>& result ) const
{
	if( children.empty()){
		result.emplace_back(pair<string,set<string>>( name,set<string>{name}));
		return result;
	}

	set<string> current = {};

	for( auto& child: children ){
		child.getGroups( result );
		current += result.back().second;
	}

	result.push_back( pair<string,set<string>>( name, current ));

	return result;
};

vector<pair<string,set<string>>> Tree::getGroups() const
{
	vector<pair<string,set<string>>> result;
	return getGroups( result );
}


/**
 * Label all nodes in a tree using a numeric label, in pre-order, recursively
 * 
 * The format of the label is root_name(default:"Node") + numeric label
 * 
 * The last_label parameter must be initialized to the label of the first node
 * during the top call.
 * 
 * NOTE: leaf nodes are not labelled, but they are counted in the pre-order
 */
void Tree::label( string root_name ) {
	unsigned int l = 0;
	this->label( l, root_name );
}

/**
 * Recursive node labelling function
 */
void Tree::label( unsigned int &last_label, const string& root_name ) {
	last_label++;

	if( children.empty()) {
		return;
	}

	name = root_name + to_string( last_label );

	for( auto& child: children ) {
		child.label( last_label, root_name );
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
