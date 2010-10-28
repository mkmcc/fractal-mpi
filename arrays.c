/* ========================================================================== *
 * arrays.c
 * Time-stamp: <2010-10-28 15:10:23 (mkmcc)>
 *
 * Functions to allocate and de-allocate arrays.  Taken verbatim from
 * the athena source code.
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
#include "utils.c"

static void free_1d_array(void *array);
static void free_2d_array(void **array);
static void free_3d_array(void *array);
static void   *calloc_1d_array(size_t nc, size_t size);
static void  **calloc_2d_array(size_t nr, size_t nc, size_t size);
static void ***calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size);

/* Functions ripped from athena
 * -------------------------------------------------------------------------- */
static void free_1d_array(void *array)
{
  free(array);
}

static void free_2d_array(void **array)
{
  free(array[0]);
  free(array);
}

static void free_3d_array(void *array)
{
  void ***ta = (void ***)array;

  free(ta[0][0]);
  free(ta[0]);
  free(array);
}

static void *calloc_1d_array(size_t nc, size_t size)
{
  void *array;
  
  if ((array = (void *)calloc(nc,size)) == NULL) {
    ath_error("[calloc_1d] failed to allocate memory (%d of size %d)\n",
              (int)nc,(int)size);
    return NULL;
  }
  return array;
}

static void **calloc_2d_array(size_t nr, size_t nc, size_t size)
{
  void **array;
  size_t i;

  if((array = (void **)calloc(nr,sizeof(void*))) == NULL){
    ath_error("[calloc_2d] failed to allocate memory for %d pointers\n",(int)nr);
    return NULL;
  }

  if((array[0] = (void *)calloc(nr*nc,size)) == NULL){
    ath_error("[calloc_2d] failed to allocate memory (%d X %d of size %d)\n",
              (int)nr,(int)nc,(int)size);
    free((void *)array);
    return NULL;
  }

  for(i=1; i<nr; i++){
    array[i] = (void *)((unsigned char *)array[0] + i*nc*size);
  }

  return array;
}

static void ***calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size)
{
  void ***array;
  size_t i,j;

  if((array = (void ***)calloc(nt,sizeof(void**))) == NULL){
    ath_error("[calloc_3d] failed to allocate memory for %d 1st-pointers\n",
              (int)nt);
    return NULL;
  }

  if((array[0] = (void **)calloc(nt*nr,sizeof(void*))) == NULL){
    ath_error("[calloc_3d] failed to allocate memory for %d 2nd-pointers\n",
              (int)(nt*nr));
    free((void *)array);
    return NULL;
  }

  for(i=1; i<nt; i++){
    array[i] = (void **)((unsigned char *)array[0] + i*nr*sizeof(void*));
  }

  if((array[0][0] = (void *)calloc(nt*nr*nc,size)) == NULL){
    ath_error("[calloc_3d] failed to alloc. memory (%d X %d X %d of size %d)\n",
              (int)nt,(int)nr,(int)nc,(int)size);
    free((void *)array[0]);
    free((void *)array);
    return NULL;
  }

  for(j=1; j<nr; j++){
    array[0][j] = (void **)((unsigned char *)array[0][j-1] + nc*size);
  }

  for(i=1; i<nt; i++){
    array[i][0] = (void **)((unsigned char *)array[i-1][0] + nr*nc*size);
    for(j=1; j<nr; j++){
      array[i][j] = (void **)((unsigned char *)array[i][j-1] + nc*size);
    }
  }

  return array;
}
