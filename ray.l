%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "y.tab.h"  /* created via yacc -d */
int lineno;         /* lineno and yyin initialized in ray.yl */
char *stripq(char *s);
%}
%option nounput
%%
[0-9]+                              {yylval.i = atoi(yytext); return INT;}
[0-9]+[\.][0-9]*                    {yylval.d = atof(yytext); return REAL;}
[0-9]+[\.]?[0-9]*[eE][\+\-]?[0-9]+  {yylval.d = atof(yytext); return REAL;}
LOOKAT                              {return LOOKAT;}
PROJECTION                          {return PROJ;}
IMAGE                               {return IMAGE;}
RECURSEDEPTH                        {return RECDEPTH;}
AMBIENT                             {return AMBIENT;}
BACKGROUND                          {return BACKGROUND;}
LIGHT                               {return LIGHTSRC;}
sin                                 {return SIN;}
cos                                 {return COS;}
tan                                 {return TAN;}
sqrt                                {return SQRT;}
KA                                  {return KA;}
KD                                  {return KD;}
KS                                  {return KS;}
KT                                  {return KT;}
NI                                  {return NI;}
PHONG                               {return PHONG;}
COLOR                               {return COLOR;}
MARBLE                              {return MARBLE;}
CHECKER                             {return CHECKER;}
SPHERE                              {return SPHERE;}
PLANE                               {return PLANE;}
BEZIER3                             {return BEZIER3;}
TEAPOT                              {return TEAPOT;}
ELEVMAP                             {return ELEVMAP;}
HERMITEFUNC                         {return HERMITEFUNC;}
BILINEARFIS                         {return BILINEARFIS;}
HERMITEFIS                          {return HERMITEFIS;}
BILINEARRIS                         {return BILINEARRIS;}
HERMITERIS                          {return HERMITERIS;}
SUPERELLIPSOID                      {return SUPERELLIPSOID;}
[,=\(\)\{\}\+\-\*\/]                {return yytext[0];}
\".*\"                              {yylval.s = stripq(yytext); return STR;}
[a-zA-Z\_][a-zA-Z0-9\_]*            {yylval.s = strdup(yytext); return IDENT;}
[ \t\r]                             { ; }
\n                                  {lineno++;}
\#     {
         int c;
         while ((c = input()) != EOF && c != '\n')
           ;
         lineno++;
       }
%%

char *stripq(char *s) {
  s = strdup(s+1);
  s[strlen(s)-1] = '\0';
  return s;
}
