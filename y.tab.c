/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



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




/* Copy the first part of user declarations.  */
#line 1 "ray.y"

#include <stdlib.h>
#include <math.h>

extern int yylex(void);
void yyerror(char *msg);
void *myalloc(size_t sz);

typedef struct List {  /* base type for all lists */
  struct List *next;
} List;

int length(List *list);
List *insert(List *list, List *node);  /* insert onto tail */

typedef struct VecList {
  struct VecList *next;
  double v[3];
} VecList;

typedef struct VecListList {
  struct VecListList *next;
  VecList *list;
} VecListList;

typedef struct NumList {
  struct NumList *next;
  double num;
} NumList;

typedef struct Matrix {
  struct Matrix *next;
  NumList *row;
} Matrix;

VecList *newVecList(double v[3]);
VecListList *newVecListList(VecList *);
NumList *newNumList(double num);
Matrix *newMatrixRow(NumList *row);

void go(void);
void setLookAt(double eye[3], double lookat[3], double up[3]);
void setProjection(double fovy, double aspect);
void setImage(char *fname, int w, int h, int samples);
void setRecursionDepth(int depth);
void setAmbient(double ambient[3]);
void setBackground(double background[3]);
void setLight(double pos[3], double color[3], 
              double c0, double c1, double c2);
void setVariable(char *name, double val);
double lookupVariable(char *name);
void setMaterialAmbient(double);
void setMaterialDiffuse(double);
void setMaterialSpecular(double);
void setMaterialTransmission(double);
void setMaterialRefraction(double);
void setMaterialPhong(double);
void setMaterialColor(double color[3]);
void setMarbleTexture(VecList *colors, double veinDir[3], int octaves);
void setCheckerTexture(void);
void makeSphere(double center[3], double rad);
void makePlane(double normal[3], double point[3]);
void makeBezier3(VecListList *mesh);
void makeTeapot(double org[3], double up[3], double spout[3], double scale);
void makeElevationMap(double orgxy[], double widthheight[], 
                      double zminmax[], char *fname);
void makeHermiteFunc(double x[], double y[], 
                     int W, int H,
                     double *f, double *fx, double *fy, double *fxy);
void makeBilinearFis(double x[], double y[],
		     int W, int H,
		     double **alpha, double **z);
void makeHermiteFis(double x[], double y[],
		    int W, int H,
		    double **alpha, 
		    double **z, double **zx, double **zy, double **zxy);
void makeBilinearRis(int M, int N, int Sx, int Sy, double x0, double y0,
	             double width, double height, int W, int H,
		     double alpha, Matrix *z);
void makeHermiteRis(int M, int N, int Sx, int Sy, double x0, double y0,
                    double width, double height, int W, int H, double alpha, 
                    Matrix *z, Matrix *zx, Matrix *zy, Matrix *zxy);
void makeSuperellipsoid(double n, double m, 
			double size[3], double z[3], double x[3], double center[3]);
#ifdef YYDEBUG
int yydebug=1;
#endif


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

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
/* Line 193 of yacc.c.  */
#line 274 "y.tab.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 287 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  69
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   567

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  49
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  25
/* YYNRULES -- Number of rules.  */
#define YYNRULES  73
/* YYNRULES -- Number of states.  */
#define YYNSTATES  340

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   293

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      45,    46,    40,    38,    44,    39,     2,    41,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    43,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    47,     2,    48,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    42
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     7,    10,    12,    14,    16,    18,
      20,    22,    24,    26,    28,    30,    38,    44,    54,    62,
      68,    72,    76,    80,    92,    98,   106,   108,   112,   116,
     122,   128,   134,   144,   154,   158,   160,   162,   166,   169,
     172,   176,   180,   184,   188,   193,   198,   203,   208,   210,
     212,   214,   218,   222,   228,   232,   236,   240,   244,   248,
     252,   256,   266,   268,   274,   280,   284,   294,   304,   322,
     336,   356,   384,   424
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      50,     0,    -1,    51,    -1,    52,    -1,    52,    51,    -1,
      53,    -1,    54,    -1,    55,    -1,    56,    -1,    57,    -1,
      58,    -1,    59,    -1,    67,    -1,    72,    -1,    73,    -1,
       7,    43,    60,    44,    60,    44,    60,    -1,     8,    43,
      68,    44,    68,    -1,     9,    43,     3,    44,    68,    44,
      68,    44,    68,    -1,     9,    43,     3,    44,    68,    44,
      68,    -1,     9,    43,     3,    44,    68,    -1,    10,    43,
      68,    -1,    11,    43,    60,    -1,    12,    43,    60,    -1,
      13,    43,    60,    44,    60,    44,    68,    44,    68,    44,
      68,    -1,    13,    43,    60,    44,    60,    -1,    45,    68,
      44,    68,    44,    68,    46,    -1,    60,    -1,    61,    44,
      60,    -1,    47,    61,    48,    -1,    62,    44,    47,    61,
      48,    -1,    45,    68,    44,    68,    46,    -1,    45,    63,
      44,    63,    46,    -1,    45,    68,    44,    68,    44,    68,
      44,    68,    46,    -1,    45,    65,    44,    65,    44,    65,
      44,    65,    46,    -1,     4,    43,    68,    -1,     4,    -1,
      69,    -1,    45,    68,    46,    -1,    39,    68,    -1,    38,
      68,    -1,    68,    38,    68,    -1,    68,    39,    68,    -1,
      68,    40,    68,    -1,    68,    41,    68,    -1,    14,    45,
      68,    46,    -1,    15,    45,    68,    46,    -1,    16,    45,
      68,    46,    -1,    17,    45,    68,    46,    -1,     5,    -1,
       6,    -1,    68,    -1,    70,    44,    68,    -1,    47,    70,
      48,    -1,    71,    44,    47,    70,    48,    -1,    18,    43,
      68,    -1,    19,    43,    68,    -1,    20,    43,    68,    -1,
      21,    43,    68,    -1,    22,    43,    68,    -1,    23,    43,
      68,    -1,    24,    43,    60,    -1,    25,    43,    45,    61,
      46,    44,    60,    44,    68,    -1,    26,    -1,    27,    43,
      60,    44,    68,    -1,    28,    43,    60,    44,    60,    -1,
      29,    43,    62,    -1,    30,    43,    60,    44,    60,    44,
      60,    44,    68,    -1,    31,    43,    63,    44,    63,    44,
      63,    44,     3,    -1,    32,    43,    63,    44,    63,    44,
      68,    44,    68,    44,    64,    44,    64,    44,    64,    44,
      64,    -1,    33,    43,    63,    44,    63,    44,    68,    44,
      68,    44,    66,    44,    66,    -1,    34,    43,    63,    44,
      63,    44,    68,    44,    68,    44,    66,    44,    66,    44,
      66,    44,    66,    44,    66,    -1,    35,    43,    68,    44,
      68,    44,    68,    44,    68,    44,    68,    44,    68,    44,
      68,    44,    68,    44,    68,    44,    68,    44,    68,    44,
      47,    71,    48,    -1,    36,    43,    68,    44,    68,    44,
      68,    44,    68,    44,    68,    44,    68,    44,    68,    44,
      68,    44,    68,    44,    68,    44,    68,    44,    47,    71,
      48,    44,    47,    71,    48,    44,    47,    71,    48,    44,
      47,    71,    48,    -1,    37,    43,    68,    44,    68,    44,
      60,    44,    60,    44,    60,    44,    60,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   124,   124,   127,   128,   131,   132,   133,   134,   135,
     136,   137,   138,   139,   140,   143,   146,   149,   152,   155,
     160,   163,   166,   169,   171,   174,   178,   179,   184,   185,
     191,   194,   199,   205,   211,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   229,   230,
     233,   234,   239,   240,   245,   246,   247,   248,   249,   250,
     251,   252,   254,   257,   258,   259,   260,   262,   265,   270,
     274,   279,   287,   297
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "STR", "IDENT", "INT", "REAL", "LOOKAT",
  "PROJ", "IMAGE", "RECDEPTH", "AMBIENT", "BACKGROUND", "LIGHTSRC", "SIN",
  "COS", "TAN", "SQRT", "KA", "KD", "KS", "KT", "NI", "PHONG", "COLOR",
  "MARBLE", "CHECKER", "SPHERE", "PLANE", "BEZIER3", "TEAPOT", "ELEVMAP",
  "HERMITEFUNC", "BILINEARFIS", "HERMITEFIS", "BILINEARRIS", "HERMITERIS",
  "SUPERELLIPSOID", "'+'", "'-'", "'*'", "'/'", "UMINUS", "'='", "','",
  "'('", "')'", "'{'", "'}'", "$accept", "args", "params", "param",
  "lookat", "projection", "image", "rec_depth", "ambient", "background",
  "light", "vec", "vec_list", "vec_llist", "pair", "mat2x2", "a4",
  "mat4x4", "assign", "expr", "num", "num_list", "matrix", "material",
  "object", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,    43,    45,
      42,    47,   293,    61,    44,    40,    41,   123,   125
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    49,    50,    51,    51,    52,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    53,    54,    55,    55,    55,
      56,    57,    58,    59,    59,    60,    61,    61,    62,    62,
      63,    64,    65,    66,    67,    68,    68,    68,    68,    68,
      68,    68,    68,    68,    68,    68,    68,    68,    69,    69,
      70,    70,    71,    71,    72,    72,    72,    72,    72,    72,
      72,    72,    72,    73,    73,    73,    73,    73,    73,    73,
      73,    73,    73,    73
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     7,     5,     9,     7,     5,
       3,     3,     3,    11,     5,     7,     1,     3,     3,     5,
       5,     5,     9,     9,     3,     1,     1,     3,     2,     2,
       3,     3,     3,     3,     4,     4,     4,     4,     1,     1,
       1,     3,     3,     5,     3,     3,     3,     3,     3,     3,
       3,     9,     1,     5,     5,     3,     9,     9,    17,    13,
      19,    27,    39,    13
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    62,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       2,     3,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    14,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     1,
       4,    35,    48,    49,     0,     0,     0,     0,     0,     0,
       0,    34,    36,     0,     0,     0,     0,    20,    21,    22,
       0,    54,    55,    56,    57,    58,    59,    60,     0,     0,
       0,     0,    65,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    39,    38,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    26,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    37,    40,    41,
      42,    43,     0,     0,    16,    19,    24,     0,     0,    63,
      64,    28,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    44,    45,    46,    47,     0,     0,     0,     0,
      27,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    15,    18,     0,     0,    29,     0,    30,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    25,
      17,     0,    61,    66,    67,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    23,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      69,     0,     0,     0,    73,     0,     0,     0,     0,     0,
       0,     0,    31,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    68,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    33,    70,     0,
       0,    32,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    50,     0,     0,
      71,     0,     0,    52,     0,     0,    51,     0,     0,    53,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    72
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,   128,   129,   102,   105,   240,   250,   242,    39,   317,
      82,   318,   315,    40,    41
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -294
static const yytype_int16 yypact[] =
{
     263,   -36,   -29,   -14,     0,     8,    16,    32,    33,    44,
      59,    64,    65,    66,    82,    86,    88,  -294,    92,    96,
      97,    98,    99,   100,   106,   112,   114,   115,   117,    56,
    -294,   263,  -294,  -294,  -294,  -294,  -294,  -294,  -294,  -294,
    -294,  -294,    89,    51,    89,   129,    89,    51,    51,    51,
      89,    89,    89,    89,    89,    89,    51,   116,    51,    51,
     130,    51,   133,   133,   133,   133,    89,    89,    89,  -294,
    -294,  -294,  -294,  -294,   134,   136,   142,   143,    89,    89,
      89,    28,  -294,    89,   160,    77,   162,    28,  -294,  -294,
     164,    28,    28,    28,    28,    28,    28,  -294,    51,   166,
     174,    51,   176,   178,    89,   179,   180,   182,   183,   161,
     218,   264,    89,    89,    89,    89,  -294,  -294,    73,    89,
      89,    89,    89,   271,    51,    89,    89,    51,  -294,   -16,
      89,    51,   -43,   172,    51,   278,   133,   133,   133,   133,
      89,    89,    89,   113,   124,   145,   157,  -294,    13,    13,
    -294,  -294,    89,   186,    28,   285,   191,    51,   196,    28,
    -294,  -294,    51,   203,    89,   204,   206,   211,   219,   292,
     299,   306,  -294,  -294,  -294,  -294,   313,    51,    89,    89,
    -294,    51,   -25,    51,   175,   133,    89,    89,    89,    89,
      89,    51,    89,  -294,   320,   327,   222,  -294,   224,  -294,
     225,   334,   341,   348,   355,   362,   233,   193,    89,    89,
      89,    89,   234,    89,    89,    89,    89,    89,    51,  -294,
      28,   369,    28,    28,  -294,   376,   383,   390,   397,   404,
     235,    89,   256,   261,   261,    89,    89,    51,    28,   133,
     269,   262,   270,   276,   411,   418,   277,   283,   256,    89,
     284,   261,   261,    89,    89,    51,   133,   290,   425,   262,
    -294,   291,   432,   439,  -294,   232,   256,    89,   297,   261,
      89,    89,  -294,   298,   446,   262,   304,   453,   460,   256,
      89,   305,   261,    89,    89,  -294,   467,   262,   311,   474,
     481,    89,   310,   261,    89,    89,   205,  -294,  -294,   488,
     495,  -294,    89,    89,   502,   509,    89,    89,   516,   523,
     315,   316,   322,   322,    89,   -17,     1,    28,    17,   323,
    -294,   332,    89,  -294,    89,   330,    28,    37,   322,  -294,
      38,   339,   337,   322,    53,   346,   344,   322,    76,  -294
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -294,  -294,   366,  -294,  -294,  -294,  -294,  -294,  -294,  -294,
    -294,    -1,   -98,  -294,   -47,  -227,  -253,  -219,  -294,   -42,
    -294,    74,  -293,  -294,  -294
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint16 yytable[] =
{
      81,   157,    85,   132,    87,   161,   268,    42,    91,    92,
      93,    94,    95,    96,    43,   243,   106,   107,   108,   157,
     316,   257,   281,   197,   109,   110,   111,   319,   157,    44,
     158,   320,   260,   261,   292,   330,   116,   117,   118,   273,
     334,   123,    84,    45,   338,   319,    88,    89,    90,   321,
     276,    46,   285,   121,   122,    97,    69,    99,   100,    47,
     103,   322,   135,   288,   182,   323,   119,   120,   121,   122,
     143,   144,   145,   146,   298,    48,    49,   148,   149,   150,
     151,   322,   319,   154,   155,   329,   331,    50,   159,   165,
     166,   167,   168,    71,    72,    73,    83,   319,   169,   170,
     171,   335,    51,    74,    75,    76,    77,    52,    53,    54,
     176,   119,   120,   121,   122,   119,   120,   121,   122,   147,
     319,   125,   184,   153,   339,    55,   156,    78,    79,    56,
     160,    57,    86,   163,    80,    58,   194,   195,   200,    59,
      60,    61,    62,    63,   201,   202,   203,   204,   205,    64,
     207,   119,   120,   121,   122,    65,   180,    66,    67,   172,
      68,    98,   119,   120,   121,   122,   220,   221,   222,   223,
     173,   225,   226,   227,   228,   229,   193,   101,   104,   112,
     196,   113,   198,   119,   120,   121,   122,   114,   115,   238,
     206,   174,   247,   244,   245,   119,   120,   121,   122,   119,
     120,   121,   122,   175,   124,   140,   126,   258,   127,   265,
     130,   262,   263,   119,   120,   121,   122,   230,   131,   162,
     133,   199,   134,   136,   137,   274,   138,   139,   277,   278,
     177,   119,   120,   121,   122,   179,   246,   224,   286,   219,
     181,   289,   290,   119,   120,   121,   122,   183,   185,   296,
     186,   301,   299,   300,   264,   187,   119,   120,   121,   122,
     304,   305,   141,   188,   308,   309,   210,     1,   211,   212,
       2,     3,     4,     5,     6,     7,     8,   218,   272,   237,
     326,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,   239,   119,   120,   121,   122,   241,   249,   142,   119,
     120,   121,   122,   248,   251,   152,   119,   120,   121,   122,
     252,   255,   164,   119,   120,   121,   122,   256,   259,   178,
     119,   120,   121,   122,   266,   269,   189,   119,   120,   121,
     122,   275,   279,   190,   119,   120,   121,   122,   282,   287,
     191,   119,   120,   121,   122,   293,   297,   192,   119,   120,
     121,   122,   312,   313,   208,   119,   120,   121,   122,   314,
     324,   209,   119,   120,   121,   122,   325,   328,   213,   119,
     120,   121,   122,   332,   333,   214,   119,   120,   121,   122,
     336,   337,   215,   119,   120,   121,   122,    70,   327,   216,
     119,   120,   121,   122,     0,     0,   217,   119,   120,   121,
     122,     0,     0,   231,   119,   120,   121,   122,     0,     0,
     232,   119,   120,   121,   122,     0,     0,   233,   119,   120,
     121,   122,     0,     0,   234,   119,   120,   121,   122,     0,
       0,   235,   119,   120,   121,   122,     0,     0,   236,   119,
     120,   121,   122,     0,     0,   253,   119,   120,   121,   122,
       0,     0,   254,   119,   120,   121,   122,     0,     0,   267,
     119,   120,   121,   122,     0,     0,   270,   119,   120,   121,
     122,     0,     0,   271,   119,   120,   121,   122,     0,     0,
     280,   119,   120,   121,   122,     0,     0,   283,   119,   120,
     121,   122,     0,     0,   284,   119,   120,   121,   122,     0,
       0,   291,   119,   120,   121,   122,     0,     0,   294,   119,
     120,   121,   122,     0,     0,   295,   119,   120,   121,   122,
       0,     0,   302,   119,   120,   121,   122,     0,     0,   303,
     119,   120,   121,   122,     0,     0,   306,   119,   120,   121,
     122,     0,     0,   307,   119,   120,   121,   122,     0,     0,
     310,   119,   120,   121,   122,     0,     0,   311
};

static const yytype_int16 yycheck[] =
{
      42,    44,    44,   101,    46,    48,   259,    43,    50,    51,
      52,    53,    54,    55,    43,   234,    63,    64,    65,    44,
     313,   248,   275,    48,    66,    67,    68,    44,    44,    43,
      46,    48,   251,   252,   287,   328,    78,    79,    80,   266,
     333,    83,    43,    43,   337,    44,    47,    48,    49,    48,
     269,    43,   279,    40,    41,    56,     0,    58,    59,    43,
      61,    44,   104,   282,   162,    48,    38,    39,    40,    41,
     112,   113,   114,   115,   293,    43,    43,   119,   120,   121,
     122,    44,    44,   125,   126,    48,    48,    43,   130,   136,
     137,   138,   139,     4,     5,     6,    45,    44,   140,   141,
     142,    48,    43,    14,    15,    16,    17,    43,    43,    43,
     152,    38,    39,    40,    41,    38,    39,    40,    41,    46,
      44,    44,   164,   124,    48,    43,   127,    38,    39,    43,
     131,    43,     3,   134,    45,    43,   178,   179,   185,    43,
      43,    43,    43,    43,   186,   187,   188,   189,   190,    43,
     192,    38,    39,    40,    41,    43,   157,    43,    43,    46,
      43,    45,    38,    39,    40,    41,   208,   209,   210,   211,
      46,   213,   214,   215,   216,   217,   177,    47,    45,    45,
     181,    45,   183,    38,    39,    40,    41,    45,    45,   231,
     191,    46,   239,   235,   236,    38,    39,    40,    41,    38,
      39,    40,    41,    46,    44,    44,    44,   249,    44,   256,
      44,   253,   254,    38,    39,    40,    41,   218,    44,    47,
      44,    46,    44,    44,    44,   267,    44,    44,   270,   271,
      44,    38,    39,    40,    41,    44,   237,     3,   280,    46,
      44,   283,   284,    38,    39,    40,    41,    44,    44,   291,
      44,    46,   294,   295,   255,    44,    38,    39,    40,    41,
     302,   303,    44,    44,   306,   307,    44,     4,    44,    44,
       7,     8,     9,    10,    11,    12,    13,    44,    46,    44,
     322,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
      37,    45,    38,    39,    40,    41,    45,    45,    44,    38,
      39,    40,    41,    44,    44,    44,    38,    39,    40,    41,
      44,    44,    44,    38,    39,    40,    41,    44,    44,    44,
      38,    39,    40,    41,    44,    44,    44,    38,    39,    40,
      41,    44,    44,    44,    38,    39,    40,    41,    44,    44,
      44,    38,    39,    40,    41,    44,    46,    44,    38,    39,
      40,    41,    47,    47,    44,    38,    39,    40,    41,    47,
      47,    44,    38,    39,    40,    41,    44,    47,    44,    38,
      39,    40,    41,    44,    47,    44,    38,    39,    40,    41,
      44,    47,    44,    38,    39,    40,    41,    31,   324,    44,
      38,    39,    40,    41,    -1,    -1,    44,    38,    39,    40,
      41,    -1,    -1,    44,    38,    39,    40,    41,    -1,    -1,
      44,    38,    39,    40,    41,    -1,    -1,    44,    38,    39,
      40,    41,    -1,    -1,    44,    38,    39,    40,    41,    -1,
      -1,    44,    38,    39,    40,    41,    -1,    -1,    44,    38,
      39,    40,    41,    -1,    -1,    44,    38,    39,    40,    41,
      -1,    -1,    44,    38,    39,    40,    41,    -1,    -1,    44,
      38,    39,    40,    41,    -1,    -1,    44,    38,    39,    40,
      41,    -1,    -1,    44,    38,    39,    40,    41,    -1,    -1,
      44,    38,    39,    40,    41,    -1,    -1,    44,    38,    39,
      40,    41,    -1,    -1,    44,    38,    39,    40,    41,    -1,
      -1,    44,    38,    39,    40,    41,    -1,    -1,    44,    38,
      39,    40,    41,    -1,    -1,    44,    38,    39,    40,    41,
      -1,    -1,    44,    38,    39,    40,    41,    -1,    -1,    44,
      38,    39,    40,    41,    -1,    -1,    44,    38,    39,    40,
      41,    -1,    -1,    44,    38,    39,    40,    41,    -1,    -1,
      44,    38,    39,    40,    41,    -1,    -1,    44
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     4,     7,     8,     9,    10,    11,    12,    13,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    50,
      51,    52,    53,    54,    55,    56,    57,    58,    59,    67,
      72,    73,    43,    43,    43,    43,    43,    43,    43,    43,
      43,    43,    43,    43,    43,    43,    43,    43,    43,    43,
      43,    43,    43,    43,    43,    43,    43,    43,    43,     0,
      51,     4,     5,     6,    14,    15,    16,    17,    38,    39,
      45,    68,    69,    45,    60,    68,     3,    68,    60,    60,
      60,    68,    68,    68,    68,    68,    68,    60,    45,    60,
      60,    47,    62,    60,    45,    63,    63,    63,    63,    68,
      68,    68,    45,    45,    45,    45,    68,    68,    68,    38,
      39,    40,    41,    68,    44,    44,    44,    44,    60,    61,
      44,    44,    61,    44,    44,    68,    44,    44,    44,    44,
      44,    44,    44,    68,    68,    68,    68,    46,    68,    68,
      68,    68,    44,    60,    68,    68,    60,    44,    46,    68,
      60,    48,    47,    60,    44,    63,    63,    63,    63,    68,
      68,    68,    46,    46,    46,    46,    68,    44,    44,    44,
      60,    44,    61,    44,    68,    44,    44,    44,    44,    44,
      44,    44,    44,    60,    68,    68,    60,    48,    60,    46,
      63,    68,    68,    68,    68,    68,    60,    68,    44,    44,
      44,    44,    44,    44,    44,    44,    44,    44,    44,    46,
      68,    68,    68,    68,     3,    68,    68,    68,    68,    68,
      60,    44,    44,    44,    44,    44,    44,    44,    68,    45,
      64,    45,    66,    66,    68,    68,    60,    63,    44,    45,
      65,    44,    44,    44,    44,    44,    44,    64,    68,    44,
      66,    66,    68,    68,    60,    63,    44,    44,    65,    44,
      44,    44,    46,    64,    68,    44,    66,    68,    68,    44,
      44,    65,    44,    44,    44,    64,    68,    44,    66,    68,
      68,    44,    65,    44,    44,    44,    68,    46,    66,    68,
      68,    46,    44,    44,    68,    68,    44,    44,    68,    68,
      44,    44,    47,    47,    47,    71,    71,    68,    70,    44,
      48,    48,    44,    48,    47,    44,    68,    70,    47,    48,
      71,    48,    44,    47,    71,    48,    44,    47,    71,    48
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 124 "ray.y"
    {go();}
    break;

  case 15:
#line 143 "ray.y"
    {setLookAt((yyvsp[(3) - (7)].v), (yyvsp[(5) - (7)].v), (yyvsp[(7) - (7)].v));}
    break;

  case 16:
#line 146 "ray.y"
    {setProjection((yyvsp[(3) - (5)].d), (yyvsp[(5) - (5)].d));}
    break;

  case 17:
#line 150 "ray.y"
    {setImage((yyvsp[(3) - (9)].s), (int) (yyvsp[(5) - (9)].d), (int) (yyvsp[(7) - (9)].d), (int) (yyvsp[(9) - (9)].d)); 
			      free((yyvsp[(3) - (9)].s));}
    break;

  case 18:
#line 153 "ray.y"
    {setImage((yyvsp[(3) - (7)].s), (int) (yyvsp[(5) - (7)].d), (int) (yyvsp[(7) - (7)].d), 1); 
			      free((yyvsp[(3) - (7)].s));}
    break;

  case 19:
#line 156 "ray.y"
    {setImage((yyvsp[(3) - (5)].s), (int) (yyvsp[(5) - (5)].d), -1, 1); 
			      free((yyvsp[(3) - (5)].s));}
    break;

  case 20:
#line 160 "ray.y"
    {setRecursionDepth((int) (yyvsp[(3) - (3)].d));}
    break;

  case 21:
#line 163 "ray.y"
    {setAmbient((yyvsp[(3) - (3)].v));}
    break;

  case 22:
#line 166 "ray.y"
    {setBackground((yyvsp[(3) - (3)].v));}
    break;

  case 23:
#line 170 "ray.y"
    {setLight((yyvsp[(3) - (11)].v), (yyvsp[(5) - (11)].v), (yyvsp[(7) - (11)].d), (yyvsp[(9) - (11)].d), (yyvsp[(11) - (11)].d));}
    break;

  case 24:
#line 171 "ray.y"
    {setLight((yyvsp[(3) - (5)].v), (yyvsp[(5) - (5)].v), 1, 0, 0);}
    break;

  case 25:
#line 175 "ray.y"
    {(yyval.v)[0] = (yyvsp[(2) - (7)].d); (yyval.v)[1] = (yyvsp[(4) - (7)].d); (yyval.v)[2] = (yyvsp[(6) - (7)].d);}
    break;

  case 26:
#line 178 "ray.y"
    {(yyval.veclist) = newVecList((yyvsp[(1) - (1)].v));}
    break;

  case 27:
#line 179 "ray.y"
    {(yyval.veclist) = (VecList *) 
                                            insert((List *) (yyvsp[(1) - (3)].veclist), 
	                                           (List *) newVecList((yyvsp[(3) - (3)].v)));}
    break;

  case 28:
#line 184 "ray.y"
    {(yyval.vecllist) = newVecListList((yyvsp[(2) - (3)].veclist));}
    break;

  case 29:
#line 185 "ray.y"
    {(yyval.vecllist) = (VecListList *)
                                                  insert((List *) (yyvsp[(1) - (5)].vecllist),
                                                         (List *)
                                                          newVecListList((yyvsp[(4) - (5)].veclist)));}
    break;

  case 30:
#line 191 "ray.y"
    {(yyval.v)[0] = (yyvsp[(2) - (5)].d); (yyval.v)[1] = (yyvsp[(4) - (5)].d);}
    break;

  case 31:
#line 194 "ray.y"
    {(yyval.a) = myalloc(4*sizeof(double));
                                         (yyval.a)[0] = (yyvsp[(2) - (5)].v)[0]; (yyval.a)[1] = (yyvsp[(2) - (5)].v)[1];
                                         (yyval.a)[2] = (yyvsp[(4) - (5)].v)[0]; (yyval.a)[3] = (yyvsp[(4) - (5)].v)[1];}
    break;

  case 32:
#line 200 "ray.y"
    {(yyval.a) = myalloc(4*sizeof(double));
                                         (yyval.a)[0] = (yyvsp[(2) - (9)].d); (yyval.a)[1] = (yyvsp[(4) - (9)].d);
                                         (yyval.a)[2] = (yyvsp[(6) - (9)].d); (yyval.a)[3] = (yyvsp[(6) - (9)].d);}
    break;

  case 33:
#line 206 "ray.y"
    {(yyval.m) = myalloc(4*sizeof(double *));
                                         (yyval.m)[0] = (yyvsp[(2) - (9)].a); (yyval.m)[1] = (yyvsp[(4) - (9)].a);
                                         (yyval.m)[2] = (yyvsp[(6) - (9)].a); (yyval.m)[3] = (yyvsp[(6) - (9)].a);}
    break;

  case 34:
#line 211 "ray.y"
    {setVariable((yyvsp[(1) - (3)].s), (yyvsp[(3) - (3)].d)); free((yyvsp[(1) - (3)].s));}
    break;

  case 35:
#line 214 "ray.y"
    {(yyval.d) = lookupVariable((yyvsp[(1) - (1)].s)); free((yyvsp[(1) - (1)].s));}
    break;

  case 37:
#line 216 "ray.y"
    {(yyval.d) = (yyvsp[(2) - (3)].d);}
    break;

  case 38:
#line 217 "ray.y"
    {(yyval.d) = -(yyvsp[(2) - (2)].d);}
    break;

  case 39:
#line 218 "ray.y"
    {(yyval.d) = (yyvsp[(2) - (2)].d);}
    break;

  case 40:
#line 219 "ray.y"
    {(yyval.d) = (yyvsp[(1) - (3)].d) + (yyvsp[(3) - (3)].d);}
    break;

  case 41:
#line 220 "ray.y"
    {(yyval.d) = (yyvsp[(1) - (3)].d) - (yyvsp[(3) - (3)].d);}
    break;

  case 42:
#line 221 "ray.y"
    {(yyval.d) = (yyvsp[(1) - (3)].d) * (yyvsp[(3) - (3)].d);}
    break;

  case 43:
#line 222 "ray.y"
    {(yyval.d) = (yyvsp[(1) - (3)].d) / (yyvsp[(3) - (3)].d);}
    break;

  case 44:
#line 223 "ray.y"
    {(yyval.d) = sin((yyvsp[(3) - (4)].d));}
    break;

  case 45:
#line 224 "ray.y"
    {(yyval.d) = cos((yyvsp[(3) - (4)].d));}
    break;

  case 46:
#line 225 "ray.y"
    {(yyval.d) = tan((yyvsp[(3) - (4)].d));}
    break;

  case 47:
#line 226 "ray.y"
    {(yyval.d) = sqrt((yyvsp[(3) - (4)].d));}
    break;

  case 48:
#line 229 "ray.y"
    {(yyval.d) = (double) (yyvsp[(1) - (1)].i);}
    break;

  case 50:
#line 233 "ray.y"
    {(yyval.numlist) = newNumList((yyvsp[(1) - (1)].d));}
    break;

  case 51:
#line 234 "ray.y"
    {(yyval.numlist) = (NumList *)
                                          insert((List *) (yyvsp[(1) - (3)].numlist),
                                                 (List *) newNumList((yyvsp[(3) - (3)].d)));}
    break;

  case 52:
#line 239 "ray.y"
    {(yyval.matrix) = newMatrixRow((yyvsp[(2) - (3)].numlist));}
    break;

  case 53:
#line 240 "ray.y"
    {(yyval.matrix) = (Matrix *)
                                             insert((List *) (yyvsp[(1) - (5)].matrix),
                                                  (List *) newMatrixRow((yyvsp[(4) - (5)].numlist)));}
    break;

  case 54:
#line 245 "ray.y"
    {setMaterialAmbient((yyvsp[(3) - (3)].d));}
    break;

  case 55:
#line 246 "ray.y"
    {setMaterialDiffuse((yyvsp[(3) - (3)].d));}
    break;

  case 56:
#line 247 "ray.y"
    {setMaterialSpecular((yyvsp[(3) - (3)].d));}
    break;

  case 57:
#line 248 "ray.y"
    {setMaterialTransmission((yyvsp[(3) - (3)].d));}
    break;

  case 58:
#line 249 "ray.y"
    {setMaterialRefraction((yyvsp[(3) - (3)].d));}
    break;

  case 59:
#line 250 "ray.y"
    {setMaterialPhong((yyvsp[(3) - (3)].d));}
    break;

  case 60:
#line 251 "ray.y"
    {setMaterialColor((yyvsp[(3) - (3)].v));}
    break;

  case 61:
#line 253 "ray.y"
    {setMarbleTexture((yyvsp[(4) - (9)].veclist), (yyvsp[(7) - (9)].v), (yyvsp[(9) - (9)].d));}
    break;

  case 62:
#line 254 "ray.y"
    {setCheckerTexture();}
    break;

  case 63:
#line 257 "ray.y"
    {makeSphere((yyvsp[(3) - (5)].v), (yyvsp[(5) - (5)].d));}
    break;

  case 64:
#line 258 "ray.y"
    {makePlane((yyvsp[(3) - (5)].v), (yyvsp[(5) - (5)].v));}
    break;

  case 65:
#line 259 "ray.y"
    {makeBezier3((yyvsp[(3) - (3)].vecllist));}
    break;

  case 66:
#line 261 "ray.y"
    {makeTeapot((yyvsp[(3) - (9)].v), (yyvsp[(5) - (9)].v), (yyvsp[(7) - (9)].v), (yyvsp[(9) - (9)].d));}
    break;

  case 67:
#line 263 "ray.y"
    {makeElevationMap((yyvsp[(3) - (9)].v), (yyvsp[(5) - (9)].v), 
                                                          (yyvsp[(7) - (9)].v), (yyvsp[(9) - (9)].s));}
    break;

  case 68:
#line 267 "ray.y"
    {makeHermiteFunc((yyvsp[(3) - (17)].v), (yyvsp[(5) - (17)].v), 
                                                         (int) (yyvsp[(7) - (17)].d), (int) (yyvsp[(9) - (17)].d), 
                                                         (yyvsp[(11) - (17)].a), (yyvsp[(13) - (17)].a), (yyvsp[(15) - (17)].a), (yyvsp[(17) - (17)].a));}
    break;

  case 69:
#line 272 "ray.y"
    {makeBilinearFis((yyvsp[(3) - (13)].v), (yyvsp[(5) - (13)].v), (yyvsp[(7) - (13)].d), (yyvsp[(9) - (13)].d),
							 (yyvsp[(11) - (13)].m), (yyvsp[(13) - (13)].m));}
    break;

  case 70:
#line 276 "ray.y"
    {makeHermiteFis((yyvsp[(3) - (19)].v), (yyvsp[(5) - (19)].v), (yyvsp[(7) - (19)].d), (yyvsp[(9) - (19)].d),
                                                        (yyvsp[(11) - (19)].m), 
                                                        (yyvsp[(13) - (19)].m), (yyvsp[(15) - (19)].m), (yyvsp[(17) - (19)].m), (yyvsp[(19) - (19)].m));}
    break;

  case 71:
#line 282 "ray.y"
    {makeBilinearRis((yyvsp[(3) - (27)].d), (yyvsp[(5) - (27)].d), (yyvsp[(7) - (27)].d), (yyvsp[(9) - (27)].d),
                                                         (yyvsp[(11) - (27)].d), (yyvsp[(13) - (27)].d), (yyvsp[(15) - (27)].d), (yyvsp[(17) - (27)].d),
                                                         (yyvsp[(19) - (27)].d), (yyvsp[(21) - (27)].d),
                                                         (yyvsp[(23) - (27)].d), 
							 (yyvsp[(26) - (27)].matrix));}
    break;

  case 72:
#line 291 "ray.y"
    {makeHermiteRis((yyvsp[(3) - (39)].d), (yyvsp[(5) - (39)].d), (yyvsp[(7) - (39)].d), (yyvsp[(9) - (39)].d),
                                                        (yyvsp[(11) - (39)].d), (yyvsp[(13) - (39)].d), (yyvsp[(15) - (39)].d), (yyvsp[(17) - (39)].d),
                                                        (yyvsp[(19) - (39)].d), (yyvsp[(21) - (39)].d),
							(yyvsp[(23) - (39)].d),
							(yyvsp[(26) - (39)].matrix), (yyvsp[(30) - (39)].matrix), (yyvsp[(34) - (39)].matrix), (yyvsp[(38) - (39)].matrix));}
    break;

  case 73:
#line 298 "ray.y"
    {makeSuperellipsoid((yyvsp[(3) - (13)].d),(yyvsp[(5) - (13)].d),(yyvsp[(7) - (13)].v),(yyvsp[(9) - (13)].v),(yyvsp[(11) - (13)].v),(yyvsp[(13) - (13)].v));}
    break;


/* Line 1267 of yacc.c.  */
#line 2091 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 300 "ray.y"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <sys/param.h>
#include <assert.h>
#include <unistd.h>
#include "pnmio.h"
#include "raytrace.h"
#include "marble.h"
#include "checker.h"
#include "sphere.h"
#include "plane.h"
#include "bezier3.h"
#include "teapot.h"
#include "func2.h"
#include "rectgrid.h"
#include "hermite.h"
#include "bilinear.h"
#include "rectfis.h"
#include "bilinearfis.h"
#include "hermitefis.h"
#include "ris.h"
#include "superellipsoid.h"

/*
 * Bitmaps to track which parameters were set at least once.
 */
enum {
  SET_LOOKAT        = 0x00001,
  SET_PROJECTION    = 0x00002,
  SET_IMAGE         = 0x00004,
  SET_RECDEPTH      = 0x00008,
  SET_AMBIENT       = 0x00010,
  SET_BACKGROUND    = 0x00020,
  SET_MAT_AMBIENT   = 0x00080,
  SET_MAT_DIFFUSE   = 0x00100,
  SET_MAT_SPECULAR  = 0x00200,
  SET_MAT_TRANS     = 0x00400,
  SET_MAT_REFRACT   = 0x00800,
  SET_MAT_PHONG     = 0x01000,  
  SET_MAT_COLOR     = 0x02000  
};

static int setParams = 0;

/*
 * Arguments to raytrace().
 */
typedef struct {
  double eyePos[3];      /* position of eye */
  double eyeDir[3];      /* view direction */
  double eyeUp[3];       /* eye up vector */
  double projWidth;      /* width of projection plane */
  double projHeight;     /* height of projection plane */
  double projDist;       /* distannce from eye to projection plane */
  int imageWidth;        /* width (in pixels) of image */
  int imageHeight;       /* height of image */
  int samples;           /* sqrt of samples per pixel */
  int maxDepth;          /* maximum recursion depth (1,2,..) */
  double ambient[3];     /* RGB intensities of ambient light */
  double background[3];  /* RGB background light intensities */
  int numLights;         /* number of point light sources */
  LIGHT *lights;         /* array of point light sources */
  int numObjects;        /* number of objects in world */
  OBJECT **objects;      /* array of objects in world */
} RAYTRACE_ARGS;

static RAYTRACE_ARGS args;

#define MAX_LIGHTS 100
static int numLights = 0;
static LIGHT lights[MAX_LIGHTS];

#define MAX_OBJECTS 200
static int numObjects = 0;
static OBJECT *objects[MAX_OBJECTS];

static char *filename;
static char *imagename;
extern int lineno;      /* lineno declared in ray.l */

void yyerror(char *msg) {
  fprintf(stderr, "%s [%d] : %s\n", filename, lineno, msg);
  exit(-1);
}

void fatalError(char *msg) {
  fprintf(stderr, "%s : %s\n", filename, msg);
  exit(-1);
}

void warning(char *msg) {
  fprintf(stderr, "%s : %s\n", filename, msg);
}

static void checkParams(void) {
  if ((setParams & SET_LOOKAT) == 0)
    fatalError("Look-at transformation not specified!");
  if ((setParams & SET_PROJECTION) == 0)
    fatalError("Perspective projection not specified!");
  if ((setParams & SET_IMAGE) == 0)
    fatalError("Output image not specified!");
  if (numLights <= 0)
    fatalError("No lights specified!");
  if (numObjects <= 0)
    fatalError("No objects specified!");

  if ((setParams & SET_RECDEPTH) == 0)
    warning("Ray tracing recursion depth not specified!");
  if ((setParams & SET_AMBIENT) == 0)
    warning("Ambient light in scene not specified!");
  if ((setParams & SET_BACKGROUND) == 0)
    warning("Background color not specified!");
  if ((setParams & SET_MAT_AMBIENT) == 0)
    warning("Ambient material coefficient never specified!");
  if ((setParams & SET_MAT_DIFFUSE) == 0)
    warning("Diffuse material coefficient never specified!");
  if ((setParams & SET_MAT_SPECULAR) == 0)
    warning("Specular material coefficient never specified!");
  if ((setParams & SET_MAT_TRANS) == 0)
    warning("Transmission material coefficient never specified!");
  if ((setParams & SET_MAT_REFRACT) == 0)
    warning("Index of refraction material coefficient never specified!");
  if ((setParams & SET_MAT_PHONG) == 0)
    warning("Phong material exponent never specified!");
}

void go(void) {
  FILE *f;
  pnm_image *image;

#ifdef VERBOSE
  printf("go()\n");
#endif

  checkParams();

  if ((f = fopen(imagename, "wb")) == NULL) {
    perror(imagename);
    exit(-1);
  }

  if (args.imageHeight < 1) {  /* match aspect ratio of proj screen */
    double aspect = args.projWidth/args.projHeight;
    args.imageHeight = (int) (args.imageWidth/aspect + 0.5);
  }

  args.numLights = numLights;
  args.lights = lights;
  args.numObjects = numObjects;
  args.objects = objects;

  image = raytrace(args.eyePos, args.eyeDir, args.eyeUp,
		   args.projWidth, args.projHeight, args.projDist,
		   args.imageWidth, args.imageHeight, args.samples,
		   args.maxDepth,
		   args.ambient, args.background,
		   numLights, lights,
		   numObjects, objects);

  write_pnm_image(image, f);

  fclose(f);
}

void setLookAt(double eye[3], double lookat[3], double up[3]) {
#ifdef VERBOSE
  printf("setLookAt((%f,%f,%f), (%f,%f,%f), (%f,%f,%f))\n",
	 eye[0], eye[1], eye[2],
	 lookat[0], lookat[1], lookat[2],
	 up[0], up[1], up[2]);
#endif

  args.eyePos[0] = eye[0];
  args.eyePos[1] = eye[1];
  args.eyePos[2] = eye[2];
  args.eyeDir[0] = lookat[0] - eye[0];
  args.eyeDir[1] = lookat[1] - eye[1];
  args.eyeDir[2] = lookat[2] - eye[2];
  args.eyeUp[0] = up[0];
  args.eyeUp[1] = up[1];
  args.eyeUp[2] = up[2];

  if (fabs(args.eyeDir[0]) <= EPSILON  &&
      fabs(args.eyeDir[1]) <= EPSILON  &&
      fabs(args.eyeDir[2]) <= EPSILON)
    fatalError("Eye position and lookat coincide!");

  if (fabs(args.eyeDir[1]*up[2] - args.eyeDir[2]*up[1]) <= EPSILON &&
      fabs(args.eyeDir[2]*up[0] - args.eyeDir[0]*up[2]) <= EPSILON &&
      fabs(args.eyeDir[0]*up[1] - args.eyeDir[1]*up[0]) <= EPSILON)
    fatalError("Eye direction and up coincide!");

  if ((setParams & SET_LOOKAT) != 0)
    warning("Lookat specified more than once!");
  
  setParams |= SET_LOOKAT;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define D2R(d) ((d)*M_PI/180.0)

void setProjection(double fovy, double aspect) {

#ifdef VERBOSE
  printf("setProjection(%f,%f)\n", fovy, aspect);
#endif

  if (fovy <= EPSILON)
    fatalError("Bogus field of view!");

  if (aspect <= EPSILON)
    fatalError("Bogus projection aspect ratio!");

  args.projHeight = 2*tan(D2R(fovy)/2);
  args.projWidth = aspect*args.projHeight;
  args.projDist = 1.0;

  if ((setParams & SET_PROJECTION) != 0)
    warning("Projection specified more than once!");
  
  setParams |= SET_PROJECTION;
}

void setImage(char *fname, int w, int h, int samples) {
  struct stat sb;
  int i = 0;

#ifdef VERBOSE
  printf("setImage(%s,%d,%d,%d)\n", fname, w, h, samples);
#endif

  if (setParams & SET_IMAGE)
    fatalError("Output image specified more than once!");

  imagename = strdup(fname);

  while (stat(imagename, &sb) >= 0 && i < 999) {
    static char buf[MAXPATHLEN+1];
    sprintf(buf, "%s-%03d", fname, ++i);
    imagename = buf;
  }

  if (i >= 999)
    fatalError("Too many output images with the same name!");

  if (errno != ENOENT) {
    perror(imagename);
    exit(-1);
  }

  if (w < 1)
    fatalError("Bogus image dimensions!");
  /* h < 1 is ok
   * We'll adjust according to aspect ratio of
   * projetion screen later on.
   */

  if (samples < 0 || samples > 4)
    fatalError("Bogus pixel sampling!");

  args.imageWidth = w;
  args.imageHeight = h;
  args.samples = samples;

  setParams |= SET_IMAGE;
}

void setRecursionDepth(int depth) {

#ifdef VERBOSE
  printf("setRecursionDepth(%d)\n", depth);
#endif

  if (depth < 0 || depth > 12)
    fatalError("Bogus recusions depth!");

  args.maxDepth = depth;

  if (setParams & SET_RECDEPTH)
    warning("Recursion depth set more than once!");

  setParams |= SET_RECDEPTH;
}

void setAmbient(double ambient[3]) {
  args.ambient[0] = ambient[0];
  args.ambient[1] = ambient[1];
  args.ambient[2] = ambient[2];

  if (setParams & SET_AMBIENT)
    warning("Ambient color set more than once!");

  setParams |= SET_AMBIENT;
}

void setBackground(double background[3]) {

#ifdef VERBOSE
  printf("setBackground((%f,%f,%f)\n", 
	 background[0], background[1], background[2]);
#endif

  args.background[0] = background[0];
  args.background[1] = background[1];
  args.background[2] = background[2];

  if (setParams & SET_BACKGROUND)
    warning("Background color set more than once!");

  setParams |= SET_BACKGROUND;
}

int length(List *list) {
  int n;
  List *node;
  for (n = 0, node = list; node != NULL; n++, node = node->next)
    ;
  return n;
}

List *insert(List *list, List *node) {
  List *tail;
  if (list == NULL)
    return node;
  for (tail = list; tail->next != NULL; tail = tail->next)
    ;
  tail->next = node;
  node->next = NULL;
  return list;
}

VecList *newVecList(double v[3]) {
  VecList *node = (VecList *) myalloc(sizeof(VecList));
  node->v[0] = v[0];
  node->v[1] = v[1];
  node->v[2] = v[2];
  node->next = NULL;
  return node;
}

VecListList *newVecListList(VecList *vlist) {
  VecListList *node = (VecListList *) myalloc(sizeof(VecListList));
  node->list = vlist;
  node->next = NULL;
  return node;
}

NumList *newNumList(double num) {
  NumList *list = (NumList *) myalloc(sizeof(NumList));
  list->num = num;
  list->next = NULL;
  return list;
}

Matrix *newMatrixRow(NumList *row) {
  Matrix *m = (Matrix *) myalloc(sizeof(Matrix));
  m->row = row;
  m->next = NULL;
  return m;
}

void setLight(double pos[3], double color[3], 
              double c0, double c1, double c2) {

#ifdef VERBOSE
  printf("setLight((%f,%f,%f), (%f,%f,%f), %f, %f, %f)\n", 
	 pos[0], pos[1], pos[2],
	 color[0], color[1], color[2],
	 c0, c1, c2);
#endif

  if (numLights >= MAX_LIGHTS)
    fatalError("Too many friggin' lights!");

  lights[numLights].pos[0] = pos[0];
  lights[numLights].pos[1] = pos[1];
  lights[numLights].pos[2] = pos[2];

  lights[numLights].color[0] = color[0];
  lights[numLights].color[1] = color[1];
  lights[numLights].color[2] = color[2];
  
  lights[numLights].c[0] = c0;
  lights[numLights].c[1] = c1;
  lights[numLights].c[2] = c2;

  numLights++;
}

void *myalloc(size_t sz) {
  void *p = malloc(sz);
  if (p == NULL) {
    perror("malloc()");
    exit(-1);
  }
  return p;
}

typedef struct Var {
  char *name;
  double val;
  struct Var *next;
} Var;

Var *vars = NULL;

void setVariable(char *name, double val) {
  Var *v;

#ifdef VERBOSE
  printf("setVariable(%s, %f)\n", name, val);
#endif

  for (v = vars; v != NULL; v = v->next)
    if (strcmp(name, v->name) == 0) {
      v->val = val;
      return;
    }

  v = (Var *) myalloc(sizeof(Var));

  v->name = strdup(name);

  v->val = val;
  v->next = vars;
  vars = v;
}

double lookupVariable(char *name) {
  Var *v;

#ifdef VERBOSE
  printf("lookupVariable(%s)\n", name);
#endif

  for (v = vars; v != NULL; v = v->next)
    if (strcmp(name, v->name) == 0)
      return v->val;

  yyerror("Use of unknown variable!");

  return 0.0;
}

typedef struct {
  double ka;
  double kd;
  double ks;
  double kt;
  double ni;
  double phong;
  double color[3];
} Material;

Material material = {  /* current material settings */
  0.1, 0.7, 0.2, 0.55, 1.52, 6.0, {0.2, 0.6, 1.0}
};

void setMaterialAmbient(double ka) {
#ifdef VERBOSE
  printf("setMaterialAmbient(%f)\n", ka);
#endif

  material.ka = ka;
  setParams |= SET_MAT_AMBIENT;
}

void setMaterialDiffuse(double kd) {
#ifdef VERBOSE
  printf("setMaterialDiffuse(%f)\n", kd);
#endif

  material.kd = kd;
  setParams |= SET_MAT_DIFFUSE;
}

void setMaterialSpecular(double ks) {
#ifdef VERBOSE
  printf("setMaterialSpecular(%f)\n", ks);
#endif

  material.ks = ks;
  setParams |= SET_MAT_SPECULAR;
}

void setMaterialTransmission(double kt) {
#ifdef VERBOSE
  printf("setMaterialTransmission(%f)\n", kt);
#endif

  material.kt = kt;
  setParams |= SET_MAT_TRANS;
}

void setMaterialRefraction(double ni) {
#ifdef VERBOSE
  printf("setMaterialRefraction(%f)\n", ni);
#endif

  material.ni = ni;
  setParams |= SET_MAT_REFRACT;
}

void setMaterialPhong(double phong) {
#ifdef VERBOSE
  printf("setMaterialPhong(%f)\n", phong);
#endif

  material.phong = phong;
  setParams |= SET_MAT_PHONG;
}

enum {
  USE_SOLID_COLOR,
  USE_MARBLE,
  USE_CHECKER
} colorMode = USE_SOLID_COLOR;

void setMaterialColor(double color[3]) {
#ifdef VERBOSE
  printf("setMaterialColor((%f,%f,%f))\n", color[0], color[1], color[2]);
#endif

  material.color[0] = color[0];
  material.color[1] = color[1];
  material.color[2] = color[2];
  setParams |= SET_MAT_COLOR;

  colorMode = USE_SOLID_COLOR;
}

struct {
  int n;              /* number of colors */
  double (*rgb)[3];   /* n rgb colors */
  double veinDir[3];  /* direction of marble "grain" */
  int octaves;        /* number of octaves of noise */
} marble = {0, NULL, {0,0,0}, 0};

void setMarbleTexture(VecList *colors, double veinDir[3], int octaves) {
  int i, n = length((List *) colors);
  double (*rgb)[3];
  VecList *list;
  
  rgb = (double (*)[3]) myalloc(3*n*sizeof(double));
  for (i = 0, list = colors; i < n; i++, list = list->next) {
    rgb[i][0] = list->v[0];
    rgb[i][1] = list->v[1];
    rgb[i][2] = list->v[2];
  }

  marble.n = n;
  marble.rgb = rgb;
  marble.veinDir[0] = veinDir[0];
  marble.veinDir[1] = veinDir[1];
  marble.veinDir[2] = veinDir[2];
  marble.octaves = octaves;

  colorMode = USE_MARBLE;
}

void setCheckerTexture(void) {
  colorMode = USE_CHECKER;
}

OBJECT *setMaterial(OBJECT *object) {
  object->ka = material.ka;
  object->kd = material.kd;
  object->ks = material.ks;
  object->kt = material.kt;
  object->ni = material.ni;
  object->setColor(object, material.color);

  if (colorMode == USE_MARBLE) {
    object = marbleObject(object, marble.n, marble.rgb, 
			  marble.veinDir, marble.octaves);
  } else if (colorMode == USE_CHECKER) {
    object = checkerObject(object);
  }

  return object;
}

void makeSphere(double center[3], double rad) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makeSphere((%f,%f,%f), %f)\n", 
	 center[0], center[1], center[2], rad);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");

  object = createSphereObject(center, rad);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeSuperellipsoid(double n, double m, 
			double size[3], double z[3], double x[3], double center[3]) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makeSuperellipsoid(%f,%f, (%f,%f,%f), (%f,%f,%f), (%f,%f,%f), (%f,%f,%f))\n",
	 n, m, size[0],size[1],size[2], z[0],z[1],z[2], x[0],x[1],x[2],
	 center[0],center[1],center[2]);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");

  object = createSuperellipsoidObject(n,m, size, z, x, center);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makePlane(double normal[3], double point[3]) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makePlane((%f,%f,%f), (%f,%f,%f))\n",
	 normal[0], normal[1], normal[2], point[0], point[1], point[2]);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");
    
  object = createPlaneObjectFromNormalandPoint(normal, point);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeBezier3(VecListList *mesh) {
  OBJECT *object;
  VecListList *vllist;
  POINT3 P[4][4];
  int r,c;

  if (length((List *) mesh) != 4)
    fatalError("Cubic Bezier mesh must have exactly 4 rows!");
  for (vllist = mesh; vllist != NULL; vllist = vllist->next)
    if (length((List *) vllist->list) != 4)
      fatalError("Cubic Bezier mesh must have exactly 4 columns per row!");

  for (r = 0, vllist = mesh; r < 4; r++, vllist = vllist->next) {
    VecList *vlist;
    for (c = 0, vlist = vllist->list; c < 4; c++, vlist = vlist->next) {
      P[r][c].x = vlist->v[0];
      P[r][c].y = vlist->v[1];
      P[r][c].z = vlist->v[2];
    }
  }

  object = (OBJECT *) createBezier3Object(P);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeTeapot(double org[3], double up[3], double spout[3], double scale) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makeTeapot((%f,%f,%f), (%f,%f,%f), (%f,%f,%f), %f)\n",
	 org[0], org[1], org[2], 
	 up[0], up[1], up[2], 
	 spout[0], spout[1], spout[2], scale);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");
    
  object = createTeapotObject(org, up, spout, scale);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeElevationMap(double orgxy[], double widthheight[], 
                      double zminmax[], char *fname) {
  OBJECT *object;
  pnm_image *image = read_pnm_image_from_file(fname);
  double zscale, zshift;
  double *z;
  int r,c;

  if (!IS_PGM(image))
    fatalError("Elevation map is not a PGM image!");

  zscale = (zminmax[1] - zminmax[0])/PNM_MAXVAL(image);
  zshift = zminmax[0];
  z = (double *) myalloc(PNM_NR(image)*PNM_NC(image)*sizeof(double));
  for (r = 0; r < PNM_NR(image); r++)
    for (c = 0; c < PNM_NC(image); c++)
      z[r*PNM_NC(image) + c] = PGM_PIXEL(image, r,c)*zscale + zshift;

  object = 
    createSurfaceSupportedByRectangularGrid(PNM_NC(image), PNM_NR(image),
					    orgxy[0], orgxy[1],
					    widthheight[0], widthheight[1],
					    z, PNM_NC(image));
  object = setMaterial(object);
  objects[numObjects++] = object;
}

void makeHermiteFunc(double x[], double y[], 
                     int W, int H,
                     double *f, double *fx, double *fy, double *fxy) {  
  OBJECT *object;
  double *z = (double *) myalloc(W*H*sizeof(double));
  Func2 *func = createGeneralHermiteFunc2(x[0], x[1], y[0], y[1]);
  setGeneralHermiteFunc2(func, (double (*)[2]) f,
			 (double (*)[2]) fx, (double (*)[2]) fy,
			 (double (*)[2]) fxy);
  getFunc2Samples(func, W, H, x[0], y[0], x[1] - x[0], y[1] - y[0], z, W);
  object =  createSurfaceSupportedByRectangularGrid(W, H, x[0], y[0],
						    x[1] - x[0], y[1] - y[0],
						    z, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;
}

void makeBilinearFis(double x[], double y[],
		     int W, int H,
		     double **alpha, double **z) {
  double X[4], Y[4];
  RectFis *fis;
  double *buf = (double *) myalloc(W*H*sizeof(double));
  Func2 *func;
  double dx = (x[1] - x[0])/3, dy = (y[1] - y[0])/3;
  OBJECT *object;
  X[0] = x[0];
  X[1] = x[0] + dx;
  X[2] = x[0] + 2*dx;
  X[3] = x[1];
  Y[0] = y[0];
  Y[1] = y[0] + dy;
  Y[2] = y[0] + 2*dy;
  Y[3] = y[1];
  fis = createBilinearFis(3, 3, X, Y);
  setBilinearFis(fis, z, alpha);
  func = rectFisToFunc2(fis);
  getFunc2Samples(func, W, H, x[0], y[0], x[1] - x[0], y[1] - y[0], buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x[0], y[0],
						   x[1] - x[0], y[1] - y[0],
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;  
}


void makeHermiteFis(double x[], double y[],
		    int W, int H,
		    double **alpha, 
		    double **z, double **zx, double **zy, double **zxy) {
  double X[4], Y[4];
  RectFis *fis;
  double *buf = (double *) myalloc(W*H*sizeof(double));
  Func2 *func;
  double dx = (x[1] - x[0])/3, dy = (y[1] - y[0])/3;
  OBJECT *object;
  X[0] = x[0];
  X[1] = x[0] + dx;
  X[2] = x[0] + 2*dx;
  X[3] = x[1];
  Y[0] = y[0];
  Y[1] = y[0] + dy;
  Y[2] = y[0] + 2*dy;
  Y[3] = y[1];
  fis = createHermiteFis(3, 3, X, Y);
  setHermiteFis(fis, alpha, z, zx, zy, zxy);
  func = rectFisToFunc2(fis);
  getFunc2Samples(func, W, H, x[0], y[0], x[1] - x[0], y[1] - y[0], buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x[0], y[0],
						   x[1] - x[0], y[1] - y[0],
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;  
}

void makeBilinearRis(int M, int N, int Sx, int Sy, double x0, double y0,
	             double width, double height, int W, int H,
		     double alpha, Matrix *z) {
  Func2 *alphaFunc = createConstantFunc2(alpha);
  Ris *ris;
  Func2 *func;
  double *buf;
  int i,j;
  OBJECT *object;

  assert(W > 10 && H > 10);
  buf = (double *) myalloc(W*H*sizeof(double));

  assert(M >= 2 && N >= 2 &&
	 Sx >= 2 && (M % Sx) == 0 && 
	 Sy >= 2 && (N % Sy) == 0);

  ris = createBilinearRis(M, N, Sx, Sy, x0, y0, width, height, alphaFunc);
  for (j = 0; z != NULL; z = z->next, j++) {
    NumList *n = z->row;
    assert(j <= N);
    for (i = 0; n != NULL; n = n->next, i++) {
      assert(i <= M);
      ris->z[j][i] = n->num;
    }
  }

  assert(i > M && j > N);
  computeBilinearRisMaps(ris);
  func = risToFunc2(ris);

  getFunc2Samples(func, W, H, x0, y0, width, height, buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x0, y0,
						   width, height, 
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;  
}

void makeHermiteRis(int M, int N, int Sx, int Sy, double x0, double y0,
                    double width, double height, int W, int H, double alpha, 
                    Matrix *z, Matrix *zx, Matrix *zy, Matrix *zxy) {
  
  Func2 *alphaFunc = createConstantFunc2(alpha);
  Ris *ris;
  Func2 *func;
  double *buf;
  int i,j;
  OBJECT *object;

  assert(W > 10 && H > 10);
  buf = (double *) myalloc(W*H*sizeof(double));

  assert(M >= 2 && N >= 2 &&
	 Sx >= 2 && (M % Sx) == 0 && 
	 Sy >= 2 && (N % Sy) == 0);

  ris = createHermiteRis(M, N, Sx, Sy, x0, y0, width, height, alphaFunc);
  for (j = 0; z != NULL; j++,
       z = z->next, zx = zx->next, zy = zy->next, zxy = zxy->next) {
    NumList *n = z->row, *nx, *ny, *nxy;
    assert(j <= N && zx != NULL && zy != NULL && zxy != NULL);
    nx = zx->row; ny = zy->row; nxy = zxy->row;
    for (i = 0; n != NULL; i++, 
	   n = n->next, nx = nx->next, ny = ny->next, nxy = nxy->next) {
      assert(i <= M && nx != NULL && ny != NULL && nxy != NULL);
      ris->z[j][i] = n->num;
      ris->zx[j][i] = nx->num;
      ris->zy[j][i] = ny->num;
      ris->zxy[j][i] = nxy->num;
    }
  }

  assert(i > M && j > N);
  computeHermiteRisMaps(ris);
  func = risToFunc2(ris);

  getFunc2Samples(func, W, H, x0, y0, width, height, buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x0, y0,
						   width, height, 
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;    
}

extern FILE *yyin;  /* from ray.l */
extern int yylex(void);
extern int yyparse(void);

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s <scene>\n", argv[0]);
    exit(-1);
  }

  if ((yyin = fopen(filename = argv[1], "r")) == NULL) {
    perror(filename);
    exit(-1);
  }

  setVariable("pi", M_PI);

  yyparse();

  fclose(yyin);

  return 0;
}

