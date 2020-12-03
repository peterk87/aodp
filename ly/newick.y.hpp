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

#ifndef YY_NEWICK_LY_NEWICK_Y_HPP_INCLUDED
# define YY_NEWICK_LY_NEWICK_Y_HPP_INCLUDED
/* Debug traces.  */
#ifndef NEWICK_DEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define NEWICK_DEBUG 1
#  else
#   define NEWICK_DEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define NEWICK_DEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined NEWICK_DEBUG */
#if NEWICK_DEBUG
extern int newick_debug;
#endif

/* Token type.  */
#ifndef NEWICK_TOKENTYPE
# define NEWICK_TOKENTYPE
  enum newick_tokentype
  {
    NUMBER = 258,
    LABEL = 259,
    SEPARATOR = 260,
    PUNCTUATION = 261,
    QUOTE = 262
  };
#endif

/* Value type.  */
#if ! defined NEWICK_STYPE && ! defined NEWICK_STYPE_IS_DECLARED
typedef int NEWICK_STYPE;
# define NEWICK_STYPE_IS_TRIVIAL 1
# define NEWICK_STYPE_IS_DECLARED 1
#endif


extern NEWICK_STYPE newick_lval;

int newick_parse (Source& source);

#endif /* !YY_NEWICK_LY_NEWICK_Y_HPP_INCLUDED  */
