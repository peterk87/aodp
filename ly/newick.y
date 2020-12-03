%{
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include "Source.h"

extern void yyerror(Source& source, const char* msg );
extern int yylex();

using namespace std;

/**
* Accummulator of strings read from input
*/
vector<string> parts;

/**
* Index of the last label 
*/
int last_label = -1;

/**
* Accumulator for the phylogeny
*/
vector<Tree> trees;

int crt_node = 1;
%}

%define api.prefix {newick_}

%parse-param { Source& source }

/*  TYPES  */

%token NUMBER
%token LABEL
%token SEPARATOR
%token PUNCTUATION
%token QUOTE

/*
	http://evolution.genetics.washington.edu/phylip/newick_doc.html

              tree ==> descendant_list [ root_label ] [ : branch_length ] ;

   descendant_list ==> ( subtree { , subtree } )

           subtree ==> descendant_list [internal_node_label] [: branch_length]
                   ==> leaf_label [: branch_length]

            root_label ==> label
   internal_node_label ==> label
            leaf_label ==> label

                 label ==> unquoted_label
                       ==> quoted_label

        unquoted_label ==> string_of_printing_characters
          quoted_label ==> ' string_of_printing_characters '

         branch_length ==> signed_number
                       ==> unsigned_number
*/
/*  GRAMMAR  */
%%

tree:
	descendant_list ';' {
		source.setTree( trees[$1]);
	}
	| descendant_list root_label ';' {
		source.setTree( trees[$1]);
	}
	| descendant_list ':' branch_length ';'	{
//	capture length
		trees[$1].len = parts[$3];
		source.setTree( trees[$1]);
	}
	| descendant_list root_label ':' branch_length ';'	{
//	capture length
		trees[$1].len = parts[$3];
		source.setTree( trees[$1]);
	}
	;

descendant_list:
	'(' subtrees ')'	{
//	this is an internal node
/* 		trees[$2].name = "Node"+to_string( crt_node++ ); */
		$$ = $2;
	}
	;

subtrees:
	subtree	{
//	a new internal node is created here and pushed forward
		trees.emplace_back( "", vector<Tree>{ trees[$1]});
		$$ = trees.size() - 1;
	}
	| subtrees ',' subtree {
//	combine internal groups and pass the result on
		trees[$1].children.emplace_back( trees[$3]);
	}
	;

/*  make a subgroup, pass on the subgroup index */
subtree:
	descendant_list
	| descendant_list internal_node_label
	| descendant_list ':' branch_length	{
//	capture length
		trees[$1].len = parts[$3];
	}
	| descendant_list internal_node_label ':' branch_length		{
//	capture length
		trees[$1].len = parts[$4];
	}
	| leaf_label	{
//	a new leaf node is created here
		trees.emplace_back( parts[$1]);
	}
	| leaf_label ':' branch_length	{
//	a new leaf node is created here
//	capture length
		trees.emplace_back( parts[$1], vector<Tree>{}, parts[$3] );
		$$ = trees.size() - 1;
	}
	;

root_label:
	label
	;

internal_node_label:
	label
	;

leaf_label:
	label
	;

label:
	LABEL { source.last_newick_token = parts.at( $1 ); }
	| NUMBER
	| '\'' phrase '\''	{
		$$ = $2;
	}
	;

phrase:
	LABEL
	| LABEL separator	{
		parts.push_back( parts[$1] + parts[$2] );
		$$ = parts.size()-1;
	}
	| LABEL separator phrase	{
		parts.push_back( parts[$1] + parts[$2] + parts[$3]);
		$$ = parts.size()-1;
	}
	;

separator:
	SEPARATOR
	| QUOTE
	| PUNCTUATION
	;

branch_length:
	NUMBER
	;

%%

/*  CODE  */
void yyerror( Source& source, const char* msg ) { // must be provided, otherwise compiler error
// do nothing here, error processing based on result of newick_parse
}

void newickClear()
{
	trees.clear();
	parts.clear();
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
