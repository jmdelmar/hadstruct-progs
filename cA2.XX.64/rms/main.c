#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <qhg.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

char *
flt_str(double x)
{
  char *s;
  asprintf(&s, "%g", x);
  unsigned long int n = (unsigned long int)index(s, '.');
  if((void *)n == NULL) {
    n = (unsigned long int)index(s, '\0');
    *((char *)n) = 'p';  
    *((char *)n+1) = '0';
    *((char *)n+2) = '\0';    
  } else {
    *((char *)n) = 'p';
  }
  return s;
}

void
usage(char *argv[])
{
  fprintf(stderr, " Usage: %s CONFIG ALPHA\n", argv[0]);
  return;
}

int
main(int argc, char *argv[])
{   
  if(argc != 3) {
    usage(argv);
    exit(1);
  }

  char *fname = argv[1];

  char *e;
  double alpha_gauss = (double)strtod(argv[2], &e);
  if(*e != '\0') {
    usage(argv);
    exit(2);
  }
  
  int dims[ND] = {128, 64, 64, 64}; // t,x,y,z
  
  int n_ape = 50;
  double alpha_ape = 0.5;

  int n_gauss = 90;

  enum qhg_fermion_bc_time bc = ANTIPERIODIC; // Also set this in tmLQCD's input
  int *procs = NULL;
  qhg_comms *comms = qhg_comms_init(procs);  

  qhg_lattice *lat = qhg_lattice_init(dims, comms);
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  qhg_gauge_field gf = qhg_gauge_field_init(lat);  

  /*
    Read config 
  */
  qhg_read_gauge_field_ildg(gf, fname);
  
  /*
    Plaquette 
  */
  double p = qhg_plaquette(gf);
  if(am_io_proc)
    printf("Plaquette = %12.10f\n", p);

  /*
    APE smear in 3-dimensions
   */
  double t0 = qhg_stop_watch(0);
  qhg_gauge_field gf_ape = qhg_gauge_field_init(lat);
  qhg_ape_smear_3d(gf_ape, gf, alpha_ape, n_ape);
  if(am_io_proc)
    printf("3D APE smear in %g sec\n", qhg_stop_watch(t0));
  
  /*
    Plaquette of smeared gauge-field
  */
  double p_ape = qhg_plaquette(gf_ape);
  if(am_io_proc)
    printf("3D APE plaquette = %12.10f\n", p_ape);
  
  /*
    Spinor fields
   */
  qhg_spinor_field A;
  qhg_spinor_field B;
  A = qhg_spinor_field_init(lat, bc);
  B = qhg_spinor_field_init(lat, bc);
      
  /*
    Smear computing the source r.m.s in each iteration
  */
  t0 = qhg_stop_watch(0);
  if(am_io_proc)
    printf("Smearing the source, alpha = %g\n", alpha_gauss);  

  int source_coords[] = {0,0,0,0};
  qhg_point_spinor_field(A, source_coords, 0, 0);

  qhg_spinor_field S[2] = {A, B};
  for(int i=0; i<n_gauss; i++) {
    double *rms;
    qhg_gauss_smear_iter(S[(i+1)%2], S[i%2], gf_ape, alpha_gauss);
    rms = qhg_spinor_field_rms(S[(i+1)%2], source_coords);
    if(am_io_proc)
      printf("Smearing iter = %4d, r.m.s = %+e\n", i+1, rms[source_coords[0]]);
    free(rms);
  }
  if(am_io_proc)
    printf("Done smearing in %g sec\n", qhg_stop_watch(t0));  

  char *propname;
  asprintf(&propname, "/gpfs/work/pr74yo/di56sof3/source_gN%da%s", n_gauss, flt_str(alpha_gauss));
  qhg_write_spinors(propname, 1, &S[n_gauss%2]);
  
  /* 
     Free spinor- and gauge-fields
  */
  qhg_spinor_field_finalize(A);
  qhg_spinor_field_finalize(B);
  qhg_gauge_field_finalize(gf);  
  qhg_gauge_field_finalize(gf_ape);  
  qhg_lattice_finalize(lat);
  return 0;
}

