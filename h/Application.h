#ifndef __Application_h__
#define __Application_h__

#include <cstring>
#include <cstdio>
#include <cstdlib>

#include <unistd.h>

#include <sstream>
#include <string>
#include <cctype>

#include "Error.h"
#include "Clock.h"
#include "util.h"

#include "Source.h"
#include "TrieSlice.h"

#include "Trie.h"
#include "TrieAmbig.h"
#include "Reference.h"
#include "Match.h"

class Application
{
public:
	Application( int iArgc, char* asArgv[] );
	~Application();

	int run();
	int _run( Trie& trie );

protected:
	static const int I_oligo_size_min;
	static const int I_oligo_size_max;

private:
	/**
	 * aodp version; defined in c/version.cpp
	 * WARNING: the version constant is automatically overwritten
	 * when making a release ('make release' in the development root)
	 */
	static const string version;

	static const string dash;

	/**
	 * man page content, defined in aodp.help.cpp
	 * WARNING: the content of the aodp.help.cpp is generated from the POD file man/pod/aodp
	 * when making a release ('make release' in the development root)
	 */
	static const string _help;

	void readArguments();

	void readSequences();

	void help();
	void printVersion();

	void printOligoStrings( ostream* o, Trie& trie ) const;

	void printOligoPositions( ostream* o, Trie& trie ) const;
	void printOligoRanges( ostream* o, Trie& trie ) const;

	void printFasta( ostream* o, Trie& trie ) const;
	void printGff( ostream* o, Trie& trie ) const;
	void printTab( ostream* o, Trie& trie ) const;

	void printNewick( ostream* o, Trie& trie ) const;
	void printNodeList( ostream* o, Trie& trie ) const;
	void printLineage( ostream* o, Trie& trie ) const;
	void printCladogram( Trie& trie ) const;

// 	Clusters
	void printClusterList( ostream* o, Trie& trie ) const;
	void printClusterOligos( ostream* o, Trie& trie ) const;
	void printSequenceClusters( ostream* o, Trie& trie ) const; // EXPERIMENTAL

	void printClusterShape( ostream* o, Trie& trie ) const; // EXPERIMENTAL

// 	TEST
	void printMetrics( ostream* o, Trie& trie );
	void printSource( ostream* o, Trie& trie );

	bool testDriver();

// 	TrieT processor;

	pair<map<string, map<Sequence, Cover<Position>>>,map<Sequence, Cover<Position>>> calculateRanges( Trie& trie ) const;
	vector<unsigned long> fillRange( unsigned long lo, unsigned long hi, int first_site_gap, int inter_site_gap );

	string s_first_input_file_name;

// 	command line arguments
	deque<string> args;

// 	bool options
	map<string,bool> flags;

//	numeric integer options
	map<string,int> integers;

//	numeric floating point (double) options
	map<string,double> floats;

// 	numeric range options
	map<string,pair<int,int>> ranges;

// 	names, keys in the name_options map
	const set<string> names;
	map<string,string> name_options;

// 	output streams
	map<string,ostream*> output;
	set<string> output_files;

// 	input files
	map<string, string> input;

//	everything else is sequence files
	vector<string> sequence_files;

//	Measure duration of steps
	Clock timer;
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
