/* ========================================================================== *
 * color_conversion.c
 * Time-stamp: <2010-10-28 15:08:11 (mkmcc)>
 *
 * Functions to convert between RGB and HSV color triplets.  I think
 * they work, but haven't tested them rigorously.
 *
 * Source: http://en.literateprograms.org/RGB_to_HSV_color_space_conversion_(C)
 * (MIT/X11 License)
 *
 * Notes:
 *   1. Colors are stored in double[3] arrays.
 *   2. r,g,b are elements of [0,1]
 *   3. h is an element of [0,360]; s, v are elements of [0,1]
 *
 *
 * This file was begun by Mike McCourt Sept. 2010
 *
 * ========================================================================== */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Color conversion
 * -------------------------------------------------------------------------- */
void rgb_to_hsv(double *rgb, double *hsv); /* overwrites hsv */
void hsv_to_rgb(double *hsv, double *rgb); /* overwrites rgb */

#define MIN3(c)  ((c[1]) <= (c[2]) ? ((c[0]) <= (c[1]) ? (c[0]) : (c[1])) \
                                   : ((c[0]) <= (c[2]) ? (c[0]) : (c[2])))

#define MAX3(c)  ((c[1]) >= (c[2]) ? ((c[0]) >= (c[1]) ? (c[0]) : (c[1])) \
		                   : ((c[0]) >= (c[2]) ? (c[0]) : (c[2])))



/* RGB to HSV conversion function
 * -------------------------------------------------------------------------- */
void rgb_to_hsv(double *rgb1, double *hsv)
{
  double rgb[3];
  rgb[0] = rgb1[0];
  rgb[1] = rgb1[1];
  rgb[2] = rgb1[2];

  double rgb_min = MIN3(rgb), rgb_max = MAX3(rgb);
  hsv[2] = rgb_max;

  if (hsv[2] == 0) {		/* Black */
    hsv[0] = hsv[1] = 0;
    return;
  }

  hsv[1] = rgb_max - rgb_min;
  if (hsv[1] == 0) {
    hsv[0] = 0;
    return;
  }

  /* Normalize saturation to 1 */
  rgb[0] = (rgb[0] - rgb_min)/(rgb_max - rgb_min);
  rgb[1] = (rgb[1] - rgb_min)/(rgb_max - rgb_min);
  rgb[2] = (rgb[2] - rgb_min)/(rgb_max - rgb_min);
  rgb_min = MIN3(rgb);
  rgb_max = MAX3(rgb);

  /* Compute hue */
  if (rgb_max == rgb[0]) {
    hsv[0] = 0.0 + 60.0*(rgb[1] - rgb[2]);
    if (hsv[0] < 0.0) {
      hsv[0] += 360.0;
    }
  } else if (rgb_max == rgb[1]) {
    hsv[0] = 120.0 + 60.0*(rgb[2] - rgb[0]);
  } else {
    hsv[0] = 240.0 + 60.0*(rgb[0] - rgb[1]);
  }

  return;
}
/* End RGB to HSV conversion
 * -------------------------------------------------------------------------- */


/* HSV to RGB conversion
 * -------------------------------------------------------------------------- */
void hsv_to_rgb(double *hsv1, double *rgb)
{
  double hsv[3];
  hsv[0] = hsv1[0];
  hsv[1] = hsv1[1];
  hsv[2] = hsv1[2];

  int i;
  double f, p, q, t;

  if(hsv[1] == 0) {			/* gray */
    rgb[0] = rgb[1] = rgb[2] = hsv[2];
    return;
  }

  hsv[0] /= 60.0;
  i = (int) floor(hsv[0]);
  f = hsv[0] - i;

  p = hsv[2] * (1 - hsv[1]);
  q = hsv[2] * (1 - hsv[1] * f);
  t = hsv[2] * (1 - hsv[1] * (1 - f));

  switch(i){
  case 0:
    rgb[0] = hsv[2];
    rgb[1] = t;
    rgb[2] = p;
    break;
  case 1:
    rgb[0] = q;
    rgb[1] = hsv[2];
    rgb[2] = p;
    break;
  case 2:
    rgb[0] = p;
    rgb[1] = hsv[2];
    rgb[2] = t;
    break;
  case 3:
    rgb[0] = p;
    rgb[1] = q;
    rgb[2] = hsv[2];
    break;
  case 4:
    rgb[0] = t;
    rgb[1] = p;
    rgb[2] = hsv[2];
    break;
  default:
    rgb[0] = hsv[2];
    rgb[1] = p;
    rgb[2] = q;
  }

  return;
}
/* End HSV to RGB conversion
 * -------------------------------------------------------------------------- */
