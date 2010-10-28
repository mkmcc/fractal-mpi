#ifndef UTILS_C
#define UTILS_C
/* ========================================================================== *
 * utils.c
 * Time-stamp: <2010-10-28 15:11:20 (mkmcc)>
 *
 * Some basic utilities.  Taken (almost) verbatim from the athena
 * source code.
 *
 * Begun by Mike McCourt Sept. 2010
 * ========================================================================== */

/* Copyright information for the original program: */
/*
 ATHENA 

 AUTHORS:
     James M. Stone, Princeton Univeristy
     Thomas A. Gardiner, Sandia National Laboratory
     Peter J. Teuben, University of Maryland
     John F. Hawley, University of Virginia

     This software was produced with support from NSF Grant AST-0113571 

 Copyright 2003 James M. Stone, Thomas A. Gardiner,
                Peter J. Teuben, John F. Hawley

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


FILE *atherr_fp(void);
void ath_error(char *fmt, ...);
char *ath_strdup(const char *in);


FILE *atherr_fp(void)
{
  return (stderr);
}

void ath_error(char *fmt, ...)
{
  va_list ap;
  FILE *atherr = atherr_fp();

  fprintf(atherr,"### Fatal error: ");   /* prefix */
  va_start(ap, fmt);              /* ap starts with string 'fmt' */
  vfprintf(atherr, fmt, ap);      /* print out on atherr */
  fflush(atherr);                 /* flush it NOW */
  va_end(ap);                     /* end varargs */

  exit(EXIT_FAILURE);
}

char *ath_strdup(const char *in)
{
  char *out = (char *)malloc((1+strlen(in))*sizeof(char));
  if(out == NULL) {
    printf("ath_strdup: failed to alloc %d\n", (int)(1+strlen(in)));
    return NULL; /* malloc failed */
  }
  return strcpy(out,in);
}

#endif
