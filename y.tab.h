/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     STR = 258,
     IDENT = 259,
     INT = 260,
     REAL = 261,
     LOOKAT = 262,
     PROJ = 263,
     IMAGE = 264,
     RECDEPTH = 265,
     AMBIENT = 266,
     BACKGROUND = 267,
     LIGHTSRC = 268,
     SIN = 269,
     COS = 270,
     TAN = 271,
     SQRT = 272,
     KA = 273,
     KD = 274,
     KS = 275,
     KT = 276,
     NI = 277,
     PHONG = 278,
     COLOR = 279,
     MARBLE = 280,
     CHECKER = 281,
     SPHERE = 282,
     PLANE = 283,
     BEZIER3 = 284,
     TEAPOT = 285,
     ELEVMAP = 286,
     HERMITEFUNC = 287,
     BILINEARFIS = 288,
     HERMITEFIS = 289,
     BILINEARRIS = 290,
     HERMITERIS = 291,
     SUPERELLIPSOID = 292,
     UMINUS = 293
   };
#endif
/* Tokens.  */
#define STR 258
#define IDENT 259
#define INT 260
#define REAL 261
#define LOOKAT 262
#define PROJ 263
#define IMAGE 264
#define RECDEPTH 265
#define AMBIENT 266
#define BACKGROUND 267
#define LIGHTSRC 268
#define SIN 269
#define COS 270
#define TAN 271
#define SQRT 272
#define KA 273
#define KD 274
#define KS 275
#define KT 276
#define NI 277
#define PHONG 278
#define COLOR 279
#define MARBLE 280
#define CHECKER 281
#define SPHERE 282
#define PLANE 283
#define BEZIER3 284
#define TEAPOT 285
#define ELEVMAP 286
#define HERMITEFUNC 287
#define BILINEARFIS 288
#define HERMITEFIS 289
#define BILINEARRIS 290
#define HERMITERIS 291
#define SUPERELLIPSOID 292
#define UMINUS 293




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 89 "ray.y"
{
  char *s;
  int i;
  double d;
  double v[3];
  double *a;
  double **m;
  struct VecList *veclist;
  struct VecListList *vecllist;
  struct NumList *numlist;
  struct Matrix *matrix;
}
/* Line 1529 of yacc.c.  */
#line 138 "y.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

