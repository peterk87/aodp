%{

#include <iostream>

#include "Source.h"

using namespace std;

#include "tax.y.hpp"

%}

%option prefix="tax_"
%option noyywrap nodefault nounput

%%

[\t]	{
	return TAB;
}

[\n\r]+	{
	return ENTER;
}

;	{
	return SEMICOLON;
}

s__[^;\n\r]*	{
	tax_lval._string = strdup( yytext+3 );
	return SPECIES;
}

[a-z]+__[^;\n\r]*	{
	return CLADE;
}

[A-Za-z0-9_.]+	{
	tax_lval._string = strdup( yytext );
	return NAME;
}

			/* ignore everything else*/
.

%%

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
