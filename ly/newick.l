%{
	#include <iostream>
	#include <string>
	#include <vector>

	#include "Source.h"

	#include "newick.y.hpp"

	using namespace std;

	extern vector<string> parts;

	void add_part( const char * const text = yytext );	
%}

%option prefix="newick_"
%option noyywrap nodefault nounput

/*
	Source
	http://evolution.genetics.washington.edu/phylip/newick_doc.html

	tree ==> descendant_list [ root_label ] [ : branch_length ] ;

	descendant_list ==> ( subtree { , subtree } )

	subtree ==> descendant_list [internal_node_label] [: branch_length]
			==> leaf_label [: branch_length]

	root_label ==> label
	internal_node_label ==> label
	leaf_label ==> label

	label	==> unquoted_label
			==> quoted_label

	unquoted_label ==> string_of_printing_characters
	quoted_label ==> ' string_of_printing_characters '

	branch_length	==> signed_number
					==> unsigned_number
*/


/*  DECLARATIONS  */

/*	Unquoted labels may not contain blanks, parentheses, square brackets,
	single_quotes, colons, semicolons, or commas.	*/
LABEL	[^ \t()\[\]':;,\n\r]+

SEPARATOR	[ \t]
PUNCTUATION	[':;()\[\],]

/*    Underscore characters in unquoted labels are converted to blanks. */
/*    Single quote characters in a quoted label are represented by two single quotes. */
/*    Blanks or tabs may appear anywhere except within unquoted labels or branch_lengths. */
/*    Newlines may appear anywhere except within labels or branch_lengths. */
/*    Comments are enclosed in square brackets and may appear anywhere newlines are permitted. */

/*  Number  */
DIGIT	[0-9]
DIGIT_1_9	[1-9]

INT	{DIGIT}|{DIGIT_1_9}{DIGIT}+
SIGNED	[+-]?{INT}|{INT}
DECIMAL	{SIGNED}|{SIGNED}.[0-9]+
NUMBER	{DECIMAL}|{DECIMAL}[Ee]{SIGNED}

/*  ACTIONS  */
%%

{NUMBER}	{
/* 	cout << "NUMBER(" << yytext << ")" << endl; */
	add_part();
	return NUMBER;
}

{LABEL}	{
/* 	cout << "LABEL(" << yytext << ")" << endl; */
	add_part();
	return LABEL;
}

{SEPARATOR} 	{
/* 	cout << "SEPARATOR" << endl; */
	add_part();
/* 	 TODO: test separator inside (ok) and outside (fail) of quoted label */
	return SEPARATOR;
}

''	{
/* 	cout << "QUOTE" << endl; */
	add_part( "'" );
/* 	 TODO: test quote inside (ok) and outside (fail) */
	return QUOTE;
}

{PUNCTUATION}	{
/* 	cout << "PUNCTUATION(" << yytext << ")" << endl; */
/* 	 TODO: test punctuation inside quoted label */
	return *yytext;
}

[ \t\n\r]+          /* eat up whitespace TODO: only in allowed areas */

.           {
	error(
		"cannot process Newick tree file\n",
		string{"** encountered invalid character: '"}+yytext+"'\n",
		"** Newick tree format reference: \n",
		"   http://evolution.genetics.washington.edu/phylip/newick_doc.html"
	);
}


%%

void add_part( const char * const text )
{
	parts.emplace_back( text );
	newick_lval = parts.size()-1;
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
