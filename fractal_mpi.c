/*============================================================================*
 * fractal.c
 * Time-stamp: <2010-10-28 15:15:52 (mkmcc)>
 *
 * Draws Julia set for F = z^2 + c using DEM/J algorithm.
 *
 * Compile: "gcc -W -Wall -std=c99 -pedantic -O3 \
 *               -o fractal.x fractal.c -lm"
 *
 * Usage: Takes no input; a PPM image, which you can convert to other
 *        formats using imagemagick.
 *
 * Notes:
 *
 *   1. The constant C[2] defines the Julia set; this notation is
 *      standard.
 *
 *   2. center and [xy]m(in||ax) define the plot range.  Julia sets
 *      for this function typically fit within [-2:2]^2
 *
 *   3. I think the results should be pretty independent of iter_max
 *      and M... tolerance does influence the output, though.
 *
 *   4. I'm antialiasing the image by simply supersampling... this is
 *      a bad idea for fractals, which have structure on arbitrarily
 *      small scales, but it seems to work.  A better (but not that
 *      much better) thing to do would be to randomly sample N points
 *      in each cell and average the results.
 *
 *   5. As far as I can tell, MPI_Send and MPI_Recv only work with 1D
 *      arrays, so I had to store color in a rather awkward way.
 *      Maybe there's a better way to do this?
 *
 *   6. This code is parallelized, with pretty good scaling.  results
 *      below:
 *
 *          Parallelization (1920 x 1200 image):
 *           1 processor:  157.536 s
 *           2 processors:  85.501 s (1.84)
 *           4 processors:  45.518 s (3.46)
 *           8 processors:  25.917 s (6.08)
 *
 * TODO:
 *
 *   1. There is an awkward fence-post bug here... there are identical
 *      rows of pixels above and below a processor boundary.  While
 *      this is technically "wrong," it is very hard to detect in the
 *      finished image.
 *
 *   2. This code is a mess!  It should really be rewritten.
 *
 *   3. I should write a usage() function.
 *
 *   4. I should rethink the antialiasing.
 *
 * Begun by Mike McCourt Aug. 2010
 *
 * Parts of this code came from the Wikipedia entry:
 *  http://en.wikipedia.org/wiki/File:Demj.jpg  (CC license)
 *
 *============================================================================*/
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.c"
#include "arrays.c"
#include "color_constants.h"
#include "color_conversion.c"
#include "par.c"

#define MPI_PARALLEL
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif /* MPI_PARALLEL */


/* Macro definitions to for efficiency.
 * -------------------------------------------------------------------------- */
#define one_log2 1.44269504088896

#define MAX(x, y) (x > y) ? x : y
#define MIN(x, y) (x < y) ? x : y


/* Data structures
 * -------------------------------------------------------------------------- */
typedef struct FractalPoint_s {
  double inside;		/* Est. fraction of pixel inside the fractal */
  double outside;		/* "" outside the fractal */
}FractalPointS;

static int ins;			/* insertion point into a circular buffer */


/* Functions defining the image
 * -------------------------------------------------------------------------- */
static inline double jdist(double **Z, int *iter); /* Distance estimator */
static void render(int s, FractalPointS **sample_data);


/* Constants for the fractal iteration
 * -------------------------------------------------------------------------- */
static int iter_max = 4096;
static double C[2];
static double M  = 1e30; /* "infinity"           */
static double M2 = 1e60; /* "infinity" (squared) */


/* Image parameters
 * -------------------------------------------------------------------------- */
static int Nx, Ny, Ny_g;
static double xmin, xmax, ymin, ymax, Lx, Ly;
static double center[2];

static double **color;
static int nstripes;
static double expo, *s, base, tol;

static char *filename;


/* Other
 * -------------------------------------------------------------------------- */
static void process_input(char *athinput, int argc, char *argv[]);
static void rasterize(FractalPointS ***data, char *image);


/* MPI
 * -------------------------------------------------------------------------- */
#ifdef MPI_PARALLEL
static int rank, size;
static double sync_max (double my_var);
#endif /* MPI_PARALLEL */





/* Main function
 * -------------------------------------------------------------------------- */
int main (int argc, char *argv[])
{
  int i, j, c;

  /* image stuff */
  FILE *fp;
  char *comment  = "# ";

  FractalPointS ***array;
  char *color_image;

  char *definput = "input.frac";  /* default input filename */
  char *athinput = definput;


#ifdef MPI_PARALLEL
  /* MPI stuff */
  int r;
  MPI_Status *status;
  MPI_Init (&argc, &argv);	                /* starts MPI              */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id  */
  MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */

  status = (MPI_Status*) calloc_1d_array(size, sizeof(MPI_Status));
#endif /* MPI_PARALLEL */


  /* parse command line options */
  for (i=1; i<argc; i++) {
    switch(*(argv[i]+1)) {
    case 'i':                      /* -i <file>   */
      athinput = argv[++i];
      break;
    default:
      break;
    }
  }

  /* process the input file */
  process_input(athinput, argc, argv);


  /* render the fractal -- do this once for each stripe */
  array = (FractalPointS***) calloc_3d_array(nstripes, Ny, Nx,
					     sizeof(FractalPointS));
  for (i=0; i<nstripes; i++)
    render(s[i], array[i]);


  /* rasterize to 24-bit colors.
     this combines all the stripes into one image. */
  color_image = (char*) calloc_1d_array(Ny*Nx*3, sizeof(char));
  rasterize(array, color_image);
  free_3d_array((void***) array);


#ifdef MPI_PARALLEL
  if (rank != 0) {
    /* Send image data to the root process with a blocking send */
    MPI_Send(color_image, 3*Nx*Ny, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  } else {
#endif /* MPI_PARALLEL */

    /* Write to a file */
    fp = fopen(filename, "wb");
    fprintf(fp, "P6\n%s\n%d %d %d\n", comment, Nx, Ny_g, 255);

#ifdef MPI_PARALLEL
    for (r=0; r<size; r++) {
      if (r != 0)		/* Use a blocking recieve */
      	MPI_Recv(color_image, 3*Nx*Ny, MPI_CHAR,
		 r, 0, MPI_COMM_WORLD, &status[r]);
#endif /* MPI_PARALLEL */

      for (j=0; j<Ny; j++) {
      	for (i=0; i<Nx; i++) {
	  for (c=0; c<3; c++)
	    fwrite(&color_image[c + 3*(i + Nx*j)], 1, 1, fp);
      	}
      }
#ifdef MPI_PARALLEL
    }
    fclose(fp);			/* MPI_PARALLEL and r == 0 */
  }
#endif /* MPI_PARALLEL */


  /* Finish up */
  free_1d_array((void*) color_image);
  free_1d_array((void*) s);
  free_2d_array((void**) color);
#ifdef MPI_PARALLEL
  MPI_Finalize();
#else
  fclose(fp);
#endif /* MPI_PARALLEL */

  return 0;
}
/* End main
 * -------------------------------------------------------------------------- */



/* Process the input file
 * -------------------------------------------------------------------------- */
static void process_input(char *athinput, int argc, char *argv[])
{
  int js, je;
  double ymin_g, ymax_g;
  int i;
  char label[80];

  /* first comes first
     -------------------------------------------------- */
  par_open(athinput);
  par_cmdline(argc,argv);

  filename = par_gets_def("image", "file", "demj.ppm");

  Nx = par_geti("image", "Nx");
  Ny = par_geti("image", "Ny");
  Ny_g = Ny;



  /* fractal properties
     -------------------------------------------------- */
  center[0] = par_getd("fractal", "center_x");
  center[1] = par_getd("fractal", "center_y");

  Lx = par_getd("fractal", "Lx");
  Ly = Lx * Ny / Nx;

  xmin = center[0] - Lx/2;  xmax = center[0] + Lx/2;
  ymin = center[1] - Ly/2;  ymax = center[1] + Ly/2;

  C[0] = par_getd("fractal", "Cx");
  C[1] = par_getd("fractal", "Cy");



  /* image properties
     -------------------------------------------------- */
  /* interior points */
  expo = par_getd("image", "exp");
  base = par_getd("image", "base");
  tol  = par_getd("image", "tolerance");

  /* exterior points */
  nstripes = par_geti("image", "nstripes");

  s = (double*) calloc_1d_array(nstripes, sizeof(double));
  for (i=0; i<nstripes; i++) {
    sprintf(label, "s%d", i+1);
    s[i] = par_getd("image", label);
  }

  color = (double**) calloc_2d_array(nstripes, 3, sizeof(double));
  for (i=0; i<nstripes; i++) {
    sprintf(label, "color_%d", i+1);
    memcpy(&(color[i][0]), SetColor(par_gets("image", label)),
	   3*sizeof(double));
  }

  /* Correct Ny, ymax and ymin to reflect the local grid */
#ifdef MPI_PARALLEL
  js = (size-rank-1) * ((double) Ny/size);
  je = (size-rank) * ((double) Ny/size);

  ymin_g = ymin;  ymax_g = ymax;

  ymin = ymin_g + (ymax_g - ymin_g)*((double) js / Ny);
  ymax = ymin_g + (ymax_g - ymin_g)*((double) je / Ny);

  center[1] = ymin + (ymax - ymin)/2;

  Ny /= size;
  Ly /= size;
#endif /* MPI_PARALLEL */

  par_close();
  return;
}
/* End process input
 * -------------------------------------------------------------------------- */



/* Rendering function
 * -------------------------------------------------------------------------- */
static void render(int s, FractalPointS **sample_data)
{
  /* coordinates */
  int i, j, ip, jp, iter, c, k;

  double dx = (xmax-xmin)/Nx, dxp = dx/3;
  double dy = (ymax-ymin)/Ny, dyp = dy/3;

  double log_M = log(M);


  /* variables for the set */
  double **Z = (double**) calloc_2d_array(7, 2, sizeof(double));
  double highest[2];
  double S[4];

  double dist, tolerance = dx/tol;

  double d, rn;
  int ind;

  int nOut, nIn;

  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {

      sample_data[j][i].inside = sample_data[j][i].outside = 0.0;
      nOut = nIn = 0;

      for (jp=-1; jp<=1; jp++) {
	for (ip=-1; ip<=1; ip++) {

	  /* initialize the orbit and iterate */
	  Z[0][1] = ymax - j*dy + jp*dyp;
	  Z[0][0] = xmin + i*dx + ip*dxp;

	  dist = jdist(Z, &iter);


	  /* "outside the set:" color with field lines */
	  if (dist > tolerance) {
	    nOut += 1;

	    rn = sqrt(Z[ins][0]*Z[ins][0] + Z[ins][1]*Z[ins][1]);
	    d = log(log_M/log(rn)) * one_log2;

	    ins = (ins+1) % 7;	/* ins now points to first */
	    for (c=0; c<4; c++) {
	      S[3-c] = 0.0;
	      for (k=0; k<4; k++) {
		ind = (ins + (c+k)) % 7;
		S[3-c] += 1 + sin(s * atan2(Z[ind][0], Z[ind][1]));
	      }
	      S[3-c] *= 0.125;
	    }

	    /* smooth interpolation */
	    sample_data[j][i].outside += 0.5 *
	      ((-pow(d,2) + pow(d,3)) * S[0]
	       + (d + 4*pow(d,2) - 3*pow(d,3)) * S[1]
	       + (2 - 5*pow(d,2) + 3*pow(d,3)) * S[2]
	       + (-d + 2*pow(d,2) - pow(d,3)) * S[3]);
	  }
	  /* "inside the set:" color by proximity */
	  else {
	    nIn += 1;
	    sample_data[j][i].inside += pow(log(iter), expo);
	  }

	}
      }
      sample_data[j][i].outside /= 9.0;
      sample_data[j][i].inside /= 9.0;


      if (nOut > nIn)
	sample_data[j][i].inside = -1.0;
      else
	sample_data[j][i].outside = -1.0;
    }
  }


  /* scale outside and inside points independently */
  highest[0] = highest[1] = 0.0;
  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {
      if (sample_data[j][i].outside > highest[0])
	highest[0] = sample_data[j][i].outside;

      if (sample_data[j][i].inside > highest[1])
	highest[1] = sample_data[j][i].inside;
    }
  }

  highest[0] = MAX(1e-20, highest[0]);
  highest[1] = MAX(1e-20, highest[1]);

#ifdef MPI_PARALLEL
  highest[0] = sync_max(highest[0]);
  highest[1] = sync_max(highest[1]);
#endif /* MPI_PARALLEL */

  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {
      if (sample_data[j][i].outside > 0.0)
	sample_data[j][i].outside /= highest[0];

      if (sample_data[j][i].inside > 0.0)
	sample_data[j][i].inside /= highest[1];
    }
  }

  free_2d_array((void**) Z);
  return;
}
/* End render function
 * -------------------------------------------------------------------------- */



/* Rasterize values to 24-bit color
 * -------------------------------------------------------------------------- */
static void rasterize(FractalPointS ***data, char *image)
{
  int i, j, c;
  double **hsv_color, temp[3];

  double num, denom;

  /* Convert the rgb colors to hsv */
  hsv_color = (double**) calloc_2d_array(nstripes, 3, sizeof(double));
  for (i=0; i<nstripes; i++)
    rgb_to_hsv(color[i], hsv_color[i]);


  /* Loop over pixels */
  for (j=0; j<Ny; j++) {
    for (i=0; i<Nx; i++) {

      if (data[0][j][i].outside > 0) { /* outside the set */

	/* weighted average of hue */
	num = denom = 0.0;
	for (c=0; c<nstripes; c++) {
	  num   += hsv_color[c][0] * data[c][j][i].outside;
	  denom += data[c][j][i].outside;
	}
	temp[0] = num / denom;


	/* weighted average of saturation */
	num = denom = 0.0;
	for (c=0; c<nstripes; c++) {
	  num   += hsv_color[c][1] * data[c][j][i].outside;
	  denom += data[c][j][i].outside;
	}
	temp[1] = num / denom;


	/* product of the values (?) */
	temp[2] = 1.0;
	for (c=0; c<nstripes; c++) {
	  temp[2] *= hsv_color[c][2] * data[c][j][i].outside;
	}
	temp[2] = sqrt(temp[2]);


	/* convert back to rgb */
	hsv_to_rgb(temp, temp);
	for (c=0; c<3; c++)
	  image[c + 3*(i + Nx*j)] = (int) (255 * temp[c]);
      }
      else {			/* inside the set */
	for (c=0; c<3; c++)
	  image[c + 3*(i + Nx*j)] = (int) (255 * (base + (1.0-base)*data[0][j][i].inside));
      }
    }
  }

  return;
}
/* End rasterize
 * -------------------------------------------------------------------------- */



/* Distance estimator
 * -------------------------------------------------------------------------- */
double jdist(double **Z, int *iter)
{
  /* distance = 2 * |zn| * log|zn| / |znp|    */
  /*   z0p = 1 + 0 i                          */
  /*   znp = 2 * z(n-1) * z(n-1)p             */

  /* Z is a 1D array of complex numbers, storing the last 7 points in
     the iteration. */

  double xp = 1, yp = 0, nz, nzp, a;
  double x = Z[0][0], y = Z[0][1];

  ins = 0;

  for (*iter = 1; *iter <= iter_max; *iter+=1) {
    nz = 2*(x*xp - y*yp) + 1;
    yp = 2*(x*yp + y*xp);
    xp = nz;

    nz = x*x - y*y + C[0];
    y = 2*x*y + C[1];
    x = nz;

    nz  = x*x + y*y;
    nzp = xp*xp + yp*yp;
    if (/* nzp > M2 || */ nz > M2)
      break;

    /* "push" into the circular buffer */
    ins = (ins + 1) % 7;
    Z[ins][0] = x;
    Z[ins][1] = y;
  }

  a = sqrt(nz);
  return 2*a*log(a)/sqrt(nzp);
}
/* End distance estimator
 * -------------------------------------------------------------------------- */



/* Synchronize
 * -------------------------------------------------------------------------- */
#ifdef MPI_PARALLEL
static double sync_max (double my_var)
/* Maximize a variable across the domain. */
{
 double global_var;
 int err;

 err = MPI_Allreduce (&my_var, &global_var, 1, MPI_DOUBLE,
		      MPI_MAX, MPI_COMM_WORLD);
 if (err) {
   ath_error("[sync_max]: MPI_Allreduce returned error code %d\n",err);
 }

 return global_var;
}
#endif /* MPI_PARALLEL */
/* End synchronize
 * -------------------------------------------------------------------------- */



/*
Local Variables:
  compile-command:
    "gcc -W -Wall -std=c99 -pedantic -O3 \
         -o fractal_mpi.x fractal_mpi.c -lm -lmpi"
End:
*/
