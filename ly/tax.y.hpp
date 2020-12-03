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

#ifndef YY_TAX_LY_TAX_Y_HPP_INCLUDED
# define YY_TAX_LY_TAX_Y_HPP_INCLUDED
/* Debug traces.  */
#ifndef TAX_DEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define TAX_DEBUG 1
#  else
#   define TAX_DEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define TAX_DEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined TAX_DEBUG */
#if TAX_DEBUG
extern int tax_debug;
#endif

/* Token type.  */
#ifndef TAX_TOKENTYPE
# define TAX_TOKENTYPE
  enum tax_tokentype
  {
    TAB = 258,
    ENTER = 259,
    SEMICOLON = 260,
    CLADE = 261,
    NAME = 262,
    SPECIES = 263
  };
#endif

/* Value type.  */
#if ! defined TAX_STYPE && ! defined TAX_STYPE_IS_DECLARED
typedef union TAX_STYPE TAX_STYPE;
union TAX_STYPE
{
#line 22 "ly/tax.y" /* yacc.c:1909  */

//	C-string values
//	WARNING: these are allocated in the lexer; the parser is responsible with free-ing them
	char*	_string;

#line 77 "ly/tax.y.hpp" /* yacc.c:1909  */
};
# define TAX_STYPE_IS_TRIVIAL 1
# define TAX_STYPE_IS_DECLARED 1
#endif


extern TAX_STYPE tax_lval;

int tax_parse (Source& source);

#endif /* !YY_TAX_LY_TAX_Y_HPP_INCLUDED  */
