#define __Application_cpp__

#include "Application.h"

// version number; automatically generated from the latest git tag
// overwrites stuff between "" in the line below magic
// see the Makefile in the project root
// magic (do not modify): iMMmwyWGTTGeJ6TYQg2myA
const string Application::version = "2.5.0.1"; // do not modify

const int Application :: I_oligo_size_min = 12;
const int Application :: I_oligo_size_max = 100;

const string Application :: dash = "--";

/**
 * Call the aodp manpage, trying first the man repository, then
 * (if not found) look for it in the MAN_TARGET_DIR location.
 * 
 * Terminate the program after showing the manpage
 * 
 * NOTE: MAN_TARGET_DIR is read in "./configure" and passed as a
 * macro definition in the compile phase ("make")
 * 
 * WARNING: calls "system" to execute the manpage (potential attack vector);
 * the wordexp library is used to confirm that we are dealing with an actual
 * file that is the manpage
 */
void Application::help()
{
	cout << _help;
	exit( 0 );
}

/**
 * Display with "aodp --version", then terminate the program
 */
void Application::printVersion()
{
	cout << "aodp (Automated Oligonucleotide Design Pipeline) " << version << endl;
	cout << endl;
	cout << "(C) HER MAJESTY THE QUEEN IN RIGHT OF CANADA (2014-2018)" << endl;
	cout << "(C) Manuel Zahariev mz@alumni.sfu.ca (2000-2008,2014-2018)" << endl;
	cout << endl;
	cout << "License GPLv3: GNU GPL version 3 <http://gnu.org/licenses/gpl.html>" << endl;
	cout << "This is free software: you are free to change and redistribute it." << endl;
	cout << "There is NO WARRANTY, to the extent permitted by law." << endl;
	cout << endl;
	cout << "author: Manuel Zahariev, mz@alumni.sfu.ca" << endl;

	exit( 0 );
}

/**
 * Constructor taking the parameters of "int main(...)"
 */
Application :: Application( int iArgc, char* asArgv[] ) :
// 	arguments
	args( asArgv, asArgv+iArgc ),

// 	bool options
	flags({
		{ "crowded", false },
		{ "ambiguous-sources", true },
		{ "ambiguous-oligos", false },
		{ "ignore-snp", false },
		{ "reverse-complement", false },
	}),

//	integer options
	integers({
		{ "max-homolo", 4 },
		{ "max-ambiguities", -1 },
		{ "max-crowded-ambiguities", -1 },
		{ "first-site-gap", 5 },
		{ "inter-site-gap", 5 },
		{ "threads",
// 	the default number of threads is the "number of processors - 1" or "1" for single processor systems
			max( thread::hardware_concurrency(), unsigned( 2 )) - 1
		},

		{ "cluster-shape", 0 } // EXPERIMENTAL
	}),

// 	floating point options
	floats({
		{ "max-melting", -Thermo::K },
		{ "salt", 1 },
		{ "strand", 0.1 }
	}),

// 	numeric ranges
	ranges({
		{ "oligo-size", { 16, 16 }}
	}),

// 	string options, read
	names({
		{ "basename" },
		{ "cladogram" },
		{ "ignore-SNP" },
		{ "clusters" },
	}),

// 	output
	output({
		{ "help", &cout },
		{ "strings", &onull },
		{ "positions", &onull },
		{ "ranges", &onull },
		{ "fasta", &onull },
		{ "gff", &onull },
		{ "tab", &onull },
		{ "newick", &onull },
		{ "node-list", &onull },
		{ "lineage", &onull },
		{ "options", &onull },
		{ "time", &onull },
		{ "fold", &onull },

		{ "cluster-list", &onull },
		{ "cluster-oligos", &onull },

		{ "sequence-clusters", &onull },

		{ "metrics", &onull },
		{ "source", &onull },

		{ "match-output", &cout },
	}),

// 	input
	input({
		{ "tree-file", "" },
		{ "outgroup-file", "" },
		{ "isolation-file", "" },
		{ "names-file", "" },
		{ "database", "" },
		{ "taxonomy", "" },
		{ "match", "" },
	})

// //	everything else is sequence files
// 	,sequence_files({})
{
	testDriver();
}

/**
 * Destructor
 * 
 * Side effect: Delete (close) all open output streams
 */
Application :: ~Application()
{
	timer.check( "cleanup" );
	timer.stop();

	for( auto& o: output ){
		if(( o.second == &cout ) || ( o.second == &onull ))
			continue;

		delete o.second;
	}
}

/**
 * Read + process command line arguments.
 * 
 * In general, they look like: "--argument=value"
 * 
 * Fill in the maps: flags, integers, ranges, input and output and the sequence_files vector
 */
void Application :: readArguments()
{
	unordered_set<string> option_names;

// 	first item is the program name
	if( args.size() < 2 ) help();

	args.pop_front();

	for( string& arg: args ){
		if( arg == ( dash+"help" ))
			help();

		if( arg == ( dash+"version" ))
			printVersion();

		string dash_option;
		string prefix, suffix;

		if( !startsWith( arg, dash, dash_option )){
			addFiles( arg, sequence_files );
			continue;
		}

		split( dash_option, "=", prefix, suffix );
		lower( prefix );

		if( option_names.count( prefix ) > 0 ) { // make sure option has not been encountered before
			error( string{"duplicate option ( --"} + prefix + " )" );
		}
		option_names.emplace( prefix );

		if( flags.count( prefix )){
			bool value = true;
			if(( suffix.size() > 0 ) // zero-suffix meands "yes", like in --ignore-SNP
				&& !convert( suffix, value ))
				error( "cannot understand command line option (", arg , ")\n*** expecting yes/no" );

			flags[prefix] = value;
			continue;
		}

		if( integers.count( prefix )){
			int value;

			if( !convert( suffix, value ))
				error( "cannot understand command line option (", arg , ")\n*** expecting number" );

			if( value < 0 )
				error( "cannot understand command line option (", arg , ")\n*** expecting positive number" );

			for( const char c : suffix ) {
				if( !isdigit( c ))
					error( "cannot understand command line option (", arg , ")\n*** expecting number" );
			}

			integers[prefix] = value;
			continue;
		}

		if( floats.count( prefix )){
			double value;

			if( !convert( suffix, value ))
				error( "cannot understand command line option (", arg , ")\n*** expecting floating point number" );

			floats[prefix] = value;
			continue;
		}

		if( ranges.count( prefix )){
			string lo, hi;

			split( suffix, "-", lo, hi );

			if( hi == "" )
				hi = lo;

			int vlo, vhi;
			if( !convert( lo, vlo ))
				error( "invalid value for low range (", lo, ")" );

			if( !convert( hi, vhi ))
				error( "invalid value for high range (", hi, ")" );

			if( vlo < 0 )
				error( "cannot understand command line option ( ", arg , ")\n*** low range must be positive number" );

			if( vlo > vhi )
				error( "invalid numeric range (", lo, "-", hi, ")" );

			ranges[prefix] = { vlo, vhi };
			continue;
		}

		if( names.count( prefix )) {
			name_options.emplace( prefix, suffix );
			continue;
		}

		if( input.count( prefix )){
			input[ prefix ] = suffix;

			if( suffix == "" ) error( "empty input option (", dash+prefix, " )" );

			if( !canOpenFileRead( suffix )) error( "cannot open file (", dash+prefix+"="+suffix, " )" );
			continue;
		}

		if( output.find( prefix ) != output.end()){
			if(( suffix == "" ) || (suffix == "-" )){
				output[prefix] = &cout;
				continue;
			}

			if( contains( output_files, suffix ))
				error( "output file has already been used ( ", dash+prefix+"="+suffix, " )" );

			output_files.insert( suffix );

			output[prefix] = new ofstream( suffix );
			if( output[prefix]->good()) continue;

			error( "cannot open output file (", dash+prefix+"="+suffix, ")" );
		}

		error( "cannot understand command line option (", dash+prefix, ")" );
	}

// 	Argument tests go here

//	Whether any sequence files have been specified
	if( !sequence_files.size())
		error( "no sequence files specified. Nothing to do." );

	if( ranges.at( "oligo-size" ).first < 8 ) // issue #75
		error(
			"invalid value for --oligo-size (", ranges.at( "oligo-size" ).first, ")\n",
			" ** expecting at least 8"
		);

	if(( input.at( "taxonomy" ).size() > 0 ) && !input.at( "database" ).size()) // issue #74
		error( "--taxonomy option, but no --database specified" );

	if(( input.at( "database" ).size() > 0 ) && !input.at( "taxonomy" ).size()) // issue #74
		error( "--database option, but no --taxonomy specified" );
 
// 	If both max-crowded-ambiguities and max-ambiguities have not been specified
	if(( integers[ "max-ambiguities" ] < 0 ) && ( integers[ "max-crowded-ambiguities" ] < 0 ))
		integers[ "max-ambiguities" ] = 5; // set default for max-ambiguities

// 	Generate + validat options against "--basename"
	if( name_options.count( "basename" )){
// 		Validate generated options
		for( string s: { "strings", "fasta", "gff", "tab", "positions", "ranges", "newick", "node-list", "lineage" }){
			if( output.at( s ) != &onull )
				error( "incompatible options --basename and --"+s );
		}

		for( string s: { "cladogram" }){
			if( name_options.count( s ))
				error( "incompatible options --basename and --"+s );
		}

// 	Generate options
		for( string s: { "strings", "fasta", "gff", "tab", "positions", "ranges" }){
			output[ s ] = new ofstream( name_options[ "basename" ] + ".oligo." + s );
			if( !output[ s ]->good())
				error( "cannot open output file ( ", name_options[ "basename" ] + ".oligo." + s, ")" );
		}

		if( input["tree-file"].size()){
// 		HERE, there is an input tree-file; generate the tree-related output
			for( string s: { "newick", "node-list", "lineage" }){
				output[ s ] = new ofstream( name_options[ "basename" ] + "." + s );
				if( !output[ s ]->good())
					error( "cannot open output file ( ", name_options[ "basename" ] + s, ")" );
			}

			for( string s: { "cladogram" }){
				name_options[ s ] = name_options[ "basename" ] + "." + s + ".eps";
			}
		}
	}

// 	Generate + validat options against "--basename"
	if( name_options.count( "clusters" )){
// 		Validate generated options
		for( string s: { "cluster-list", "cluster-oligos" }){
			if( output.at( s ) != &onull )
				error( "incompatible options --clusters and --"+s );

			output[ s ] = new ofstream( name_options[ "clusters" ] + "." + s );
			if( !output[ s ]->good())
				error( "cannot open output file ( ", name_options[ "clusters" ] + "." + s, ")" );
		}
	}

	if( name_options.count( "cladogram" )) {
		if( !input[ "tree-file" ].size()) error( "--cladogram option but no --tree-file specified" );
		if( !canOpenFileWrite( name_options.at( "cladogram" ))) error( "cannot write to cladogram file (", name_options.at( "cladogram" ));
	}

	if( !integers["inter-site-gap"])
		error( "invalid value for parameter --inter-site-gap=0\n*** expecting positive number");

	if(( floats["max-melting"] > -Thermo::K ) && flags["ambiguous-oligos"] )
		error( "incompatible options --ambiguous-oligos and --max-melting" );

	if(( floats["salt"] > 1.1 ) || ( floats["salt"] < 0.05 ))
		error( "invalid value for option --salt (", floats["salt"], ")\n*** expecting value between 0.05 and 1.1" );

	if(( floats["strand"] > 100 ) || ( floats["strand"] < 0.01 ))
		error( "invalid value for option --strand (", floats["strand"], ")\n*** expecting value between 0.01 and 100" );

	for( const string s: { "newick", "node-list", "lineage" }) {
		if(( output.at( s ) != &onull ) && ( input.at( "tree-file" ).size() < 1 )) error( dash+s, "option but no --tree-file specified" );
	}

	if(( integers.at( "cluster-shape" ) < 0 ) || ( integers.at( "cluster-shape" ) > 3 ))
		error( "invalid value for experimental option --cluster-shape (", integers.at( "cluster-shape" ), ")\n*** expecting integer value between 1 and 3" );

// 	Whether any output has been specified
	bool has_output = false;
	for( auto& o: output ){
		if( o.second != &onull ){
			has_output = true;
			break;
		}
	}

	if( !has_output )
		error( "no output specified. Nothing to do." );
}

/**
 * Encapsulate the runtime processing of Application
 */
int Application :: run()
{

// 	Read command-line arguments
	readArguments();

	Source source(
		ranges["oligo-size"].first,
		ranges["oligo-size"].second,
		integers["max-ambiguities"],
		integers["max-crowded-ambiguities"],
		integers["max-homolo"],
		*output["fold"],
		flags["reverse-complement"]
	);

	timer.setOutput( output["time" ]);

// 	Read input sequences
	for( string& fn: sequence_files ){
		source.parse( fn.c_str());
	}

// 	Read phylogeny
	source.parseNewick( input["tree-file"] );

// 	Read limiters
	source.filterOutgroup( input["outgroup-file"]);
	source.readIsolationList( input["isolation-file"]);

	timer.check( "read" );

// 	Create trie, populate it, then print
	if( flags["ambiguous-oligos"] ) {
// 	Here, will create an ambiguous trie
		TrieAmbig trie_ambig( source, ranges["oligo-size"].first, ranges["oligo-size"].second );
		return _run( trie_ambig );
	}

// 	Here, will create an unambiguous trie
	Trie trie( source, ranges["oligo-size"].first, ranges["oligo-size"].second );

	return _run( trie );
}

/**
 * Populate Trie and print result
 * 
 * NOTE: same behaviour for ambiguous/non-ambiguous trie
 */
int Application::_run( Trie& trie )
{
	trie.buildSlices();
	timer.check( "prepare" );

	if(( floats["max-melting"] > -Thermo::K ) // filter on maximum melting temperature
		|| ( output["fold"] != &onull )) {    // display the secondary structure and melting temperature
		trie.source.filterMelting( integers["threads"], floats["max-melting"], floats["strand"], floats["salt"]);
		timer.check( "melt" );
	}

	trie.cover(integers["threads"]);
	timer.check( "cover" );

	if( integers["max-homolo"]>0 ){
		trie.filterHomolo( integers["threads"], integers["max-homolo"]);
		timer.check( "homolo" );
	}

	trie.touch(integers["threads"]);
	timer.check( "touch" );

	if( flags[ "ignore-snp" ]) {
		trie.smallDiff( integers["threads"]);
		timer.check( "snp" );
	}

	trie.encodeClusters();
	timer.check( "encode" );

	trie.collectClusters(integers["threads"]);
	timer.check( "clusters" );

// 	Read taxonomy file
	if( input["taxonomy"].size()) {
		trie.source.parseTaxonomy( input["taxonomy"]);
		timer.check( "taxonomy" );
	}

	if( input["database"].size()) {
		Reference reference( trie, integers["threads"] );
		reference.parse( input["database"]);
		timer.check( "reference" );
	}

	trie.collectMatches( integers["threads"]);
	timer.check( "collect" );

	trie.sortMatches(integers["threads"]);
	timer.check( "sort" );

	trie.source.printExcluded( "excluded.fasta" );

// 	Read matching file
	if( input.at( "match" ).size()) {
		Match match( *output.at( "match-output" ), integers.at( "threads" ), trie, ranges.at( "oligo-size" ).first );
		match.parse( input.at( "match" ));
		timer.check( "match" );
	}

	printOligoStrings(    output["strings"          ], trie );
	printOligoPositions(  output["positions"        ], trie );
	printOligoRanges(     output["ranges"           ], trie );
	printFasta(           output["fasta"            ], trie );
	printGff(             output["gff"              ], trie );
	printTab(             output["tab"              ], trie );
	printNewick(          output["newick"           ], trie );
	printNodeList(        output["node-list"        ], trie );
	printLineage(         output["lineage"          ], trie );

	printClusterList(     output["cluster-list"     ], trie );
	printClusterOligos(   output["cluster-oligos"   ], trie );
	printSequenceClusters(output["sequence-clusters"], trie );

	printMetrics( output["metrics"], trie );
	printSource( output["source"], trie );

	printCladogram( trie );

	if( integers.at( "cluster-shape" ) > 0 ) printClusterShape( &cout, trie ); // EXPERIMENTAL

// 	return value transfered to main
	timer.check( "print" );

	return 0;
}

/**
 * Print calculated oligo signatures as strings
 */
void Application :: printOligoStrings( ostream* o, Trie& trie ) const
{
	if( o == &onull ) return;

	for( auto& t: trie.source.getTargets()) {
		*o << "------------------------" << endl;
		*o << t.second << endl;
		*o << "------------------------" << endl;

		if( !trie.source.clusters.has( t.first )) {
			continue;
		}

		auto cluster = trie.matches.find( trie.source.clusters.at( t.first ));
		if( cluster == trie.matches.end()) {
			continue;
		}

		for( PositionDepthLength pdl: cluster->second ) {
			Position p = pdlPosition( pdl );
			Depth    d = pdlDepth( pdl );
			Length   l = pdlLength( pdl );

			for( Length x=1 ; x<=l ; x++ ){
				if(( x+d ) < trie.minim )
					continue;

				*o << trie.source.printableSubsequence( p-d, d+x );

				*o << endl;
			}
		}
	}
}

/**
 * Print calculated positions of oligo signatures
 */
void Application :: printOligoPositions( ostream* out, Trie& trie ) const
{
	if( out == &onull )
		return;

	const map<string, map<TypeFragment, Cover<Position>>>& ranges_by_target_fragment = calculateRanges( trie ).first;

	*out <<
		"Filename" << "\t" <<
		"Accession" << "\t" <<
		"Sites" << endl;

	for( const auto& e2rf : ranges_by_target_fragment ) {
			for( const auto& e2fc : e2rf.second ) {
			Position start = trie.source.fragments.at( e2fc.first ).getRange().lo();

			*out << e2rf.first << '\t' << trie.source.fragments.at( e2fc.first ).file_name;

			for( const Range<Position>& ra : e2fc.second ) {
// 	NOTE: "+1" to convert from zero-based indexes (Range)to 1-based indexes (sequence display)
				Position lo = ra.lo() - start + 1;
				Position hi = ra.hi() - start + 1;

				if( flags.at( "crowded" ))
				{
					for( Position l: ra.fill( integers.at("first-site-gap"), integers.at("inter-site-gap"))){
						*out << '\t' << l - start + 1;
					}
				}
				else{
					*out << '\t' << ( lo+hi )/2;
				}
			}

			*out << endl;
		}
	};
}

/**
 * Print calculated ranges of oligo signatures
 */
// TODO: Application :: printOligoRanges
void Application :: printOligoRanges( ostream* out, Trie& trie ) const
{
	if( out == &onull )
		return;

	map<TypeFragment, Cover<Position>> ranges_by_fragment = calculateRanges( trie ).second;

	for( const auto& e2rf : ranges_by_fragment ) {
		Position start = trie.source.fragments.at( e2rf.first ).getRange().lo();

		*out << "------------------------" << endl;
		*out << trie.source.instances.at( trie.source.instance_fragments.to.at( e2rf.first )) << " (" ;
		*out << trie.source.fragments.at( e2rf.first ).file_name;
		*out << ")" << endl;
		*out << "------------------------" << endl;

		for( const Range<Position>& ra: e2rf.second ) {
			Position lo = ra.lo() - start + 1;
			Position hi = ra.hi() - start + 1;
			Position mid = ( lo+hi )/2;

			*out << "[ ";
				*out << lo;
				*out << " - ";
				*out << hi;
			*out << " ] : ";

			*out << mid;

			if( flags.at("crowded"))
			{
				vector<Position > v = ra.fill( integers.at("first-site-gap"), integers.at("inter-site-gap"));
// 	NOTE: Range::fill always returns a vector with an ODD number of elements
				int half = ( v.size() - 1 ) / 2;

				*out << " {";

				for( int i=0 ; i<half ; i++ ){
					*out << " < " << ( v[i] - start );
				}

				*out << " < " << ( v[half] - start ) << " > ";
		
				for( unsigned int i=half+1 ; i<v.size() ; i++ ){
					*out << ( v[i] - start ) << " > ";
				}

				*out << "}";
			}

			*out << endl;
		}
	}
}

/**
 * Transform the collected oligos into ranges of base-pairs
 * 
 * \return pair:
 *   - map of target names by maps of sequence identifiers by Cover of ranges 
 *       - used for --positions : collect oligo positions (within originator sequence) for all targets
 *   - map of sequence identifiers by Cover of ranges
 *       - used for --ranges : collect oligo ranges (between two base pairs) for each sequence (irrelevant of the target)
 * 
 * \see Range, Cover
 */
pair<map<string, map<TypeFragment, Cover<Position>>>,map<TypeFragment, Cover<Position>>> Application::calculateRanges( Trie& trie ) const
{
	pair<map<string, map<TypeFragment, Cover<Position>>>,map<TypeFragment, Cover<Position>>> result;

	for( const auto& e2ta: trie.source.getTargets()) {
		if( !trie.source.clusters.has( e2ta.first )) continue; // targets with no matches do not have clusters

		auto cluster = trie.matches.find( trie.source.clusters.at( e2ta.first ));
		if( cluster == trie.matches.end()) {
			continue;
		}

		string ta = e2ta.second;

		for( PositionDepthLength pdl: cluster->second ) {
			Position p = pdlPosition( pdl );
			Depth    d = pdlDepth( pdl );
			Length   l = pdlLength( pdl );

			TypeFragment f = trie.source.getFragmentAtPosition( p );

			Range<Position> r{ p-d, Position{d}+l };

			result.first[ ta ][ f ].combineWithRange( r );
			result.second[ f ].combineWithRange( r );
		}
	}

	return result;
}

/**
 * Print calculated oligos into a file in the FASTA format
 */
void Application::printFasta( ostream* out, Trie& trie ) const
{
	if( out == &onull ) return;

	for( const auto& t: trie.source.getTargets() ){
		if( !trie.source.clusters.has( t.first )) {
			continue;
		}

		auto m = trie.matches.find( trie.source.clusters.at( t.first ));

		if( m == trie.matches.end()){
			continue;
		}

		for( const PositionDepthLength v: m->second ){
			const Position p = pdlPosition( v );
			const Depth    d = pdlDepth( v );
			const Length   l = pdlLength( v );

			const TypeFragment fid = trie.source.getFragmentAtPosition( p );
			const Fragment f = trie.source.fragments.at( fid );
			const Position s = f.getRange().lo();

			for( Length x=1 ; x<=l ; x++ ){
				if(( x+d ) < trie.minim )
					continue;

				*out
					<< ">"
					<< t.second
					<< "-len" << (d+x) << "-(s"
					<< (p-s-d)+1 << "e" << (p-s+x) << ")"
					<< f.rcId()
					<< endl;

				*out << trie.source.printableSubsequence( p-d, d+x ) << endl;
			}
		}
	}
}

/**
 * Print calculated oligos into a file in the GFF format
 */
void Application::printGff( ostream* out, Trie& trie ) const
{
	if( out == &onull ) return;

	*out << "##gff-version3" << endl;

	unsigned long id = 1;

	for( const auto& t: trie.source.getTargets() ){
		if( !trie.source.clusters.has( t.first )) {
			continue;
		}

		auto m = trie.matches.find( trie.source.clusters.at( t.first ));

		if( m == trie.matches.end()){
			continue;
		}

		for( const PositionDepthLength v: m->second ){
			const Position p = pdlPosition( v );
			const Depth    d = pdlDepth( v );
			const Length   l = pdlLength( v );

			const TypeFragment fid = trie.source.getFragmentAtPosition( p );
			const Fragment f = trie.source.fragments.at( fid );
			const Position s = f.getRange().lo();

			for( Length x=1 ; x<=l ; x++ ){
				if(( x+d ) < trie.minim )
					continue;

				*out
					<< t.second
					<< "\t.\t" << "len"
					<< "\t" << (p-s-d)+1
					<< "\t" << (p-s)+x
					<< "\t.\t+\t.\tID=" << t.second << "-" << id++
					<< f.rcId()
					<< ":" << trie.source.printableSubsequence( p-d, d+x ) << endl;
			}
		}
	}

	*out << "##FASTA" << endl;

	for( const auto& s: trie.source.instances.from ){
		*out << ">" << s.first << endl;
		for( auto i=trie.source.instance_fragments.from.equal_range( s.second ) ; i.first != i.second ; ++( i.first )) {
			*out << trie.source.printableSubsequence( trie.source.fragments.at( i.first->second ).getRange()) << endl;
		}
	}
}

/**
 * Print calculated oligos into a file in "tab" format
 */
void Application::printTab( ostream* out, Trie& trie ) const
{
	if( out == &onull ) return;

	for( const auto& t: trie.source.getTargets() ){
		if( !trie.source.clusters.has( t.first )) {
			continue;
		}

		auto m = trie.matches.find( trie.source.clusters.at( t.first ));

		if( m == trie.matches.end()){
			continue;
		}

		for( const PositionDepthLength v: m->second ){
			const Position p = pdlPosition( v );
			const Depth    d = pdlDepth( v );
			const Length   l = pdlLength( v );

			const TypeFragment fid = trie.source.getFragmentAtPosition( p );
			const Fragment f = trie.source.fragments.at( fid );
			const Position s = f.getRange().lo();

			for( Length x=1 ; x<=l ; x++ ){
				if(( x+d ) < trie.minim )
					continue;

				*out
					<< t.second
					<< "-len" << (d+x) << "-(s"
					<< (p-s-d)+1 << "e" << (p-s+x) << ")"
					<< f.rcId()
					<< "\t"
					<< trie.source.printableSubsequence( p-d, d+x )
					<< endl;
			}
		}
	}
}

/**
 * Print an annotated phylogeny tree (=including generated labels for internal nodes)
 * in the Newick format, based on the tree that has been read as input (--tree-file)
 * 
 * \see Tree
 */
void Application::printNewick( ostream* out, Trie& trie ) const
{
	if( out == &onull ) return;

	*out << trie.source.getTree() << ';' << endl; // the overloaded operator<< does not print the final ';' -- needs to be added
}

/**
 * Print a cladogram of the phylogeny, marking with '*' nodes with oligos and colouring red sections of the tree with oligos
 * 
 * \see `clado` command-line utility
 * \see Bio::Tree::Draw::Cladogram (BioPerl.org)
 * 
 * NOTE: BIN_TARGET_DIR is read in "./configure" and passed as a macro definition in the compile phase ("make")
 * 
 */
void Application::printCladogram( Trie& trie ) const
{
	if( !name_options.count( "cladogram" )) return;

	const string fn = name_options.at( "cladogram" ) + ".BxtKaD6ARVe2J4M155p16w.tmp";      // temporary file name that receives the annotated Newick tree
//	                                                    |         magic        |

	ofstream o( fn, ios_base::trunc ); // overwrite if existing
	if( !o.good()) error( "cladogram: cannot open temporary annotated treee file (", fn, ")" );

	o << trie.source.getTree().mark( trie.getNodesWithMatches(), "*" ) << ';';             // write an annotated Newick tree
	o.close();

	const string cmd = BIN_TARGET_DIR "/clado " + fn + " " + name_options.at( "cladogram" ); // command-line call for calling utility `clado`

	int result = system( cmd.c_str());
	if( !!remove( fn.c_str())) error( "cannot remove temporary file (", fn, ")");            // remove temporary file

	if( result ) error( "cannot run cladogram utility `clado` (looking in " BIN_TARGET_DIR ")" );
}

/**
 * Print a list of nodes in the philogeny, including generated labels for internal
 * nodes, based on a tree read as input (--tree-file)
 * 
 * \see Tree
 */
void Application::printNodeList( ostream* out, Trie& trie ) const
{
	if( out == &onull ) return;

	for( const auto& g: trie.source.getTree().getGroups()){
		if( g.second.size() <= 1 )
			continue;

		*out << g.first << "\t";

		bool first = true;
		for( string name: g.second ){
			if( !first ) *out << ",";
			else first = false;

			*out << name;
		}

		*out << endl;
	}
}

/**
 * Print the lineage (ordered list of ancestor nodes) of each leaf node
 * in the phylogeny
 * 
 * \see Tree
 */
void Application::printLineage( ostream* out, Trie& trie ) const
{
	if( out == &onull ) return;

	vector<string> lineage;
	for( string& line: trie.source.getTree().showLineage( lineage )){
		*out << line << endl;
	}
}

/**
 * Generate the list of clusters.
 * 
 * Definition: For a given database of sequences, a "cluster" is a grouping of sequences
 *             for which at least one oligonucleotide signature can be found,
 *             where that oligonucleotide signature matches all sequences in the grouping
 *             but does not match any sequences not in the grouping.
 * 
 * The output contains the following tab-separated columns:
 *  - Numeric identifier of the cluster. This is a generated value.
 *  - Space-separated list of all identifiers of sequences contained in the cluster
 *  - If there is a phylogeny tree and the cluster matches exactly a node
 *    in the phylogeny tree (internal or leaf), an additional column with the
 *    name of this node is included
 */
void Application::printClusterList( ostream* out, Trie& trie ) const
{
	for( auto e2cl: trie.source.clusters.from ) {
		*out << e2cl.second << '\t';

		bool first = true;
		for( Sequence s: e2cl.first ) {
			if( !first ) *out << ' ';
			else first = false;

			*out << trie.source.instances.at( s );
		}

		*out << '\t';

		if( trie.source.targets.has( e2cl.first )) {
			*out << trie.source.targets.at( e2cl.first );
		} else {
			*out << "-";
		}

		*out << endl;
	}
}

/**
 * Generate cluster oligonucleotide signatures
 * 
 * Definition: For a given database of sequences, a "cluster oligonucleotide signature"
 *             is an oligonucleotide signature that matches all sequences in a given
 *             cluster and does not match any sequences not in the cluster.
 * 
 */
void Application::printClusterOligos( ostream* out, Trie& trie ) const
{
	for( auto e2ma: trie.matches ) {
		for( const PositionDepthLength& pdl: e2ma.second ) {
			Position p = pdlPosition( pdl );
			Depth    d = pdlDepth( pdl );
			Length   l = pdlLength( pdl );

			*out
				<< e2ma.first << '\t'
				<< trie.source.printableSubsequence( p-d, d+l ) << endl;
		}
	}
}

/**
 * Prints tab-separated list of clusters associated with sequences. Fields: 
 * - numeric sequence identifier
 * - sequence name
 * 
 * - number of sequences with the same "cluster pattern" as this sequence
 * - number of sequences with the same "clade pattern" as this sequence
 * 
 * - number of sequences "hidden" by the "cluster pattern" of the current sequence (incl. current sequence)
 * - number of sequences "hidden" by the "clade pattern" of the current sequence (incl. current sequence)
 * 
 * - space-separated list of clusters containing the sequence ("cluster pattern")
 * - space-separated list of signature clades containing the sequence ("clade pattern") 
 * 
 */
void Application::printSequenceClusters( ostream* out, Trie& trie ) const {
	if( out == &onull ) return;

	unordered_map<Sequence,set<Cluster>> seq2set_clu; // clusters containing sequence
	unordered_map<Sequence,set<Cluster>> seq2set_clade; // signature clades containing sequence

	for( const auto& e2in: trie.source.instances.from ) {
		seq2set_clu.emplace( e2in.second, set<Cluster>{} );
		seq2set_clade.emplace( e2in.second, set<Cluster>{} );
	}

	for( const auto& e2cl: trie.source.clusters.from ) {
		for( const auto& se: e2cl.first ) {
			seq2set_clu.at( se ).emplace( e2cl.second );

			if( !trie.source.targets.has( e2cl.first )) continue; // skip CLUSTERS that are not CLADES
			seq2set_clade.at( se ).emplace( e2cl.second );
		}
	}

	multimap<set<Cluster>, Sequence> set_clu2seq; // set of clusters containing a sequence
	set<set<Cluster>> all_set_clu; // all possible combinations of clusters associated with each sequence

	for( const auto e2ssc: seq2set_clu ) { // first:Sequence second:set<Cluster>
		set_clu2seq.emplace( e2ssc.second, e2ssc.first );
		all_set_clu.emplace( e2ssc.second );
	}

	multimap<set<Cluster>, Sequence> set_clade2seq; // set of signature clades containing a sequence
	set<set<Cluster>> all_set_clade; // all possible combinations of signature clades associated with each sequence

	for( const auto& e2ssn: seq2set_clade ) { // first:Sequence second:set<Cluster>
		set_clade2seq.emplace( e2ssn.second, e2ssn.first );
		all_set_clade.emplace( e2ssn.second );
	}

	map<Sequence,set<Sequence>> cocluster; // for each sequence, the set of sequences "hidden" by the cluster pattern of the sequence:
	                                       // sequences with cluster patterns included in the cluster pattern of the sequence
	map<Sequence,set<Sequence>> coclade;   // for each sequence, the set of sequences "hidden" by the clade pattern of the sequence
	                                       // sequences with clade patterns included in the cluster pattern of the sequence
	
	for( const auto& e2in: trie.source.instances.from ) { // for all sequences
		const Sequence& se1 = e2in.second;
		const set<Cluster>& sc1 = seq2set_clu.at( se1 );
		const set<Cluster>& sn1 = seq2set_clade.at( se1 );

		cocluster.emplace( se1, set<Sequence>{});
		coclade.emplace( se1, set<Sequence>{});

		for( const auto& e2ssc: seq2set_clu ) {
			const Sequence& se2 = e2ssc.first;
			const set<Cluster>& sc2 = e2ssc.second;

			if( contains( sc1, sc2 )) cocluster.at( se1 ).emplace( se2 );
		}

		for( const auto& e2ssn: seq2set_clade ) {
			const Sequence& se2 = e2ssn.first;
			const set<Cluster>& sn2 = e2ssn.second;

			if( contains( sn1, sn2 )) coclade.at( se1 ).emplace( se2 );
		}
	}

// 	PRINT
	for( const auto& e2in: trie.source.instances.from ) { // for each sequence
		const string& name = e2in.first;
		const Sequence& se = e2in.second;

		*out << se << '\t'; // sequence identifier
		*out << name << '\t'; // Sequence name

		*out <<
			set_clu2seq.count( // size of the set of Sequences with the same set of Clusters
				seq2set_clu.at( se ) // set of Clusters containing the Sequence
			) << '\t';

		*out <<
			set_clade2seq.count( // size of the set of Sequences with the same set of Nodes
				seq2set_clade.at( se ) // set of signature clades containing the Sequence
			) << '\t';

		*out << cocluster.at( se ).size() << '\t';  // sequences hidden by the current sequence "cluster pattern"
		*out << coclade.at( se ).size() << '\t'; // sequences hidden by the current sequence "clade pattern"

		bool first = true;
//  - space-separated list of clusters containing the sequence ("cluster pattern")
		for( const Cluster& cl: seq2set_clu.at( se )) { // PRINT space-separated list of clusters
			if( first ) first = false; else *out << ' ';

			*out << cl;
		}
		*out << '\t';

		first = true;
//  - space-separated list of signature clades containing the sequence ("clade pattern") 
		for( const Cluster& cl: seq2set_clade.at( se )) { // PRINT space-separated list of targets (phylogeny-representative clusters)
			if( first ) first = false; else *out << ' ';
			*out << cl;
		}
		if( first ) *out << '-'; // "-" if there are no clade signatures
		*out << '\t';



		*out << endl;
	}
}

/**
 * Print metrics of the Source and Trie
 */
void Application::printMetrics( ostream* out, Trie& trie )
{
	if( out == &onull ) return;

	int nodes       = 0;
	int leaves      = 0;
	int length      = 0;
	int occurrences = 0;
	int clusters    = 0;

	Distribution nucleo_distribution;
	Distribution prefix_distribution;

	Distribution depth_distribution;
	Distribution length_distribution;
	Distribution occurrence_distribution;
	Distribution cluster_distribution;

	trie.measure(
		nodes, leaves, length, occurrences, clusters,
		nucleo_distribution, prefix_distribution, depth_distribution, length_distribution, occurrence_distribution, cluster_distribution
	);

	for( auto& m: trie.matches ) {
		cluster_distribution[ m.second.size() ]++;
	}

	*out << "===================================" << endl;
	*out << "sequences   : " << trie.source.instances.to.size() << endl;
	*out << "database    : " << trie.source.length() << endl;

	*out << "nodes       : " << nodes  << endl;
	*out << "leaves      : " << leaves << endl;
	*out << "length      : " << length << endl;
	*out << "occurrences : " << occurrences << endl;
	*out << "clusters    : " << trie.matches.size() << endl;

	*out << "=============nucleotides===========" << endl;
	for( auto& n: nucleo_distribution ) {
		*out << nu2asc.at( n.first ) << '\t' << n.second << endl;
	}

	*out << "=============depths================" << endl;
	*out << depth_distribution << endl;

	*out << "=============lengths================" << endl;
	*out << length_distribution << endl;

	*out << "=============occurrences============" << endl;
	*out << occurrence_distribution << endl;

	*out << "=============clusters===============" << endl;
	*out << cluster_distribution << endl;

	*out << "=============prefixes==============" << endl;
	for( const auto& ipr : prefix_distribution ) {
		*out << convertNu2Asc( p42nu( Prefix4( ipr.first ))) << '\t' << ipr.second << endl;
	}
}

/**
 * Print metrics of the Source and Trie
 */
void Application::printSource( ostream* out, Trie& trie )
{
	if( out == &onull ) return;
	trie.source.show( *out );
}

/**
 * Print an array of subsequences based on the value of the --cluster-shape indicating:
 *  - 1: cluster_id
 *  - 2: size of the cluster
 *  - 3: '*' if the cluster has exactly one sequence; the Node id, if the cluster matches exactly a node in the phylogeny tree (--tree-file)
 * 
 * EXPERIMENTAL
 */
void Application::printClusterShape( ostream* o, Trie& trie ) const
{
	if( integers.at( "max-homolo" ) > 0 ) // does not work with max-homolo filtering
		error( "--cluster-shape does not support --max-homolo > 0" );

	const int cluster_sequence_option = integers.at( "cluster-shape" );

	assert( cluster_sequence_option > 0 );

	for( const auto& e2in: trie.source.instances.from ) {
		*o << e2in.first << endl;

		for( Length ll = trie.minim ; ll <= trie.maxim ; ll++ ){ // display column headings
			*o << '\t' << Position( ll );
		}
		*o << endl << endl;

		for( auto i2fr = trie.source.instance_fragments.from.equal_range( e2in.second ) ; i2fr.first != i2fr.second ; ++i2fr.first ) {

			for( const Range<Position>& r: trie.source.fragments.at( i2fr.first->second ).getAmbigCompl()) {
				Position p = r.lo();

				while( true ) {
					Length l = r.cover( p, trie.minim, trie.maxim );
					if( l == 0 ) break; // no more elements in the cover
					*o << nu2asc.at( trie.source.symbol( p ));

					for( Length ll = trie.minim ; ll <= l ; ll++ ){
						Cluster c = trie.getClusterId( p, ll );

						unsigned int z = trie.source.clusters.at( c ).size();

						*o << '\t';

						switch( cluster_sequence_option ) {
							case 1:
								*o << '#' << c; // 1: cluster id
								break;
							case 2:
								*o << z; // 2: cluster size
								break;
							case 3:
								if( z==1 ) *o << '*'; // sequence oligo sig
								else if( trie.source.targets.has( trie.source.clusters.at( c )))
									*o << trie.source.targets.at( trie.source.clusters.at( c )).substr( 4 ); 
								break;
						}
					}

					*o << endl;

					++p;
				}
			}

		}
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
