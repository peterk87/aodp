/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_FASTA_LY_FASTA_Y_HPP_INCLUDED
# define YY_FASTA_LY_FASTA_Y_HPP_INCLUDED
/* Debug traces.  */
#ifndef FASTA_DEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define FASTA_DEBUG 1
#  else
#   define FASTA_DEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define FASTA_DEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined FASTA_DEBUG */
#if FASTA_DEBUG
extern int fasta_debug;
#endif

/* Token type.  */
#ifndef FASTA_TOKENTYPE
# define FASTA_TOKENTYPE
  enum fasta_tokentype
  {
    FASTA_NAME = 258,
    NUCLEOTIDE = 259,
    AMBIG = 260
  };
#endif

/* Value type.  */
#if ! defined FASTA_STYPE && ! defined FASTA_STYPE_IS_DECLARED
typedef union FASTA_STYPE FASTA_STYPE;
union FASTA_STYPE
{
#line 24 "ly/fasta.y" /* yacc.c:1909  */

//	single character values
	char	_char;
//	C-string values
//	WARNING: these are allocated in the lexer; the parser is responsible with free-ing them
	char*	_string;

#line 76 "ly/fasta.y.hpp" /* yacc.c:1909  */
};
# define FASTA_STYPE_IS_TRIVIAL 1
# define FASTA_STYPE_IS_DECLARED 1
#endif


extern FASTA_STYPE fasta_lval;

int fasta_parse (ParserFasta& parser, const string& file_name);

#endif /* !YY_FASTA_LY_FASTA_Y_HPP_INCLUDED  */
