%{

#include <iostream>
#include <string>
#include <cstdlib>

#include "Error.h"
#include "ParserFasta.h"

using namespace std;

string buffer;

extern int yylex();
extern void yyerror( ParserFasta&, const string&, const char* );

%}

%define api.prefix {fasta_}

%parse-param { ParserFasta& parser }
%parse-param { const string& file_name }

%union {
//	single character values
	char	_char;
//	C-string values
//	WARNING: these are allocated in the lexer; the parser is responsible with free-ing them
	char*	_string;
}

%token	<_string> FASTA_NAME
%token	<_char> NUCLEOTIDE
%token	<_char> AMBIG

%%

/**
 * A FASTA file is a sequence of fragments
 */
fasta : fragments {
		parser.onFinish();
	}
	;

fragments : fragment
		| fragments fragment
		;

/**
 * A fragment is a human-readable name followed by a nucleotide sequence body
 */
fragment	: name body
			;

/**
 * The name of a FASTA fragment
 */
name	: FASTA_NAME	{
//	encountering a name indicates that a new fragment is starting
			parser.addFragmentName( $1, file_name );
//	WARNING: the fragment name is allocated with strdup in fasta.l ; free it here!
			free( $1 );
		}
		;

/**
 * The body of a fragment is made from alternating ambiguous and unambiguous portions
 */
body	: nn
		| na
		;

/**
 * A portion of a fragment ending in an ambiguous section
 */
nn	: na ambig	{
//	record the ambiguous section
		parser.addAmbig( buffer );
	}
	| ambig	{
//	record the ambiguous section
		parser.addAmbig( buffer );
	}
	;

/**
 * A portion of a fragment ending in an un-ambiguous section
 */
na	: nn nucleotides	{
//	record the non-ambiguous ambiguous section
		parser.addNucleotides( buffer );
	}
	| nucleotides	{
//	record the non-ambiguous ambiguous section
		parser.addNucleotides( buffer );
	}
	;

/**
 * An unambiguous contiguous section
 */
nucleotides	: nucleotides NUCLEOTIDE	{
//	build a buffer of unambiguous nucleotides
				buffer += $2;
			}
			| NUCLEOTIDE	{
//	start building a buffer of unambiguous nucleotides
				buffer = $1;
			}
			;

/**
 * An ambiguous contiguous section
 */
ambig	:	ambig AMBIG	{
//	build a buffer of ambiguous nucleotides
				buffer += $2;
			}
			| AMBIG	{
//	start building a buffer of unambiguous nucleotides
				buffer = $1;
			}
			;

%%

void yyerror( ParserFasta& parser, const string& f, const char* s )
{
	error( "parsing error ( ", f, " ):", s );
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
