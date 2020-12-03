%{

#include <iostream>
#include <string>
#include <cstdlib>

#include "Source.h"

using namespace std;

extern int yylex();
extern void yyerror( const Source&, const char* );

static string last = "";

%}

%define api.prefix {tax_}

%parse-param { Source& source }

%union {
//	C-string values
//	WARNING: these are allocated in the lexer; the parser is responsible with free-ing them
	char*	_string;
}

//	%token	<_string> NAME
//	%token	<_string> SPECIES

%token TAB ENTER SEMICOLON CLADE
%token <_string> NAME SPECIES

%%

taxonomy	: line
			| taxonomy line
			;

line	: NAME TAB clades SPECIES end {
			last = $1;

			source.onTaxonomyEntry( $1, $4 );
			free( $1 );
			free( $4 );
		}
		;

end	: ENTER
	| SEMICOLON ENTER
	;

clades	: CLADE SEMICOLON
		| clades CLADE SEMICOLON
		;
%%

void yyerror( const Source& source, const char* s )
{
	if( !last.size()) error( "cannot parse taxonomy file\n ** failed on first line" );

	error(
		"cannot parse taxonomy file\n",
		"** last sequence id processed:", last,
		"\n ** The triggered error is likely on the following line.",
		"\n ** Common issues are invalid characters (such as spaces) in lineage names"
	);
}

// This file is part of aodp (the Automated Oligonucleotide Design Pipeline)
// 
// (C)	HER MAJESTY THE QUEEN IN RIGHT OF CANADA (2014,2015)
// (C)	Manuel Zahariev mz@alumni.sfu.ca (2000-2008,2014,2015)
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
