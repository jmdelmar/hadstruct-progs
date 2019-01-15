#include <mpi.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <qhg.h>
#include <parser.h>
#include <mg_interface.h>

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
  fprintf(stderr, " Usage: %s INPUT_FILE\n", argv[0]);
  return;
}

int
main(int argc, char *argv[])
{   
  if(argc != 2) {
    usage(argv);
    exit(1);
  }

  struct gl_params rp = parse_gl_input(argv[1]);

  int n_stout = rp.smearing.n_stout;
  double omega_stout = rp.smearing.omega_stout;
  
  qhg_comms *comms = qhg_comms_init(rp.procs);  
  qhg_lattice *lat = qhg_lattice_init(rp.dims, comms);
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  qhg_gauge_field gf = qhg_gauge_field_init(lat);  
  
  double gl_timer = qhg_stop_watch(0);

  /*
    Read config 
  */
  qhg_read_gauge_field_ildg(gf, rp.config);
  
  /*
    Plaquette 
  */
  double p = qhg_plaquette(gf);
  if(am_io_proc)
    printf("Plaquette = %12.10f\n", p);

  /*
    Stout smear in 4-dimensions
   */
  double t0 = qhg_stop_watch(0);
  qhg_gauge_field gf_stout = qhg_gauge_field_init(lat);
  qhg_stout_smear(gf_stout, gf, omega_stout, n_stout);
  if(am_io_proc)
    printf("Stout smear in %g sec\n", qhg_stop_watch(t0));
  /*
    Plaquette of smeared gauge-field
  */
  double p_stout = qhg_plaquette(gf_stout);
  if(am_io_proc)
    printf("Stout plaquette = %12.10f\n", p_stout);
  /*
    Allocate memory for gluon loops
   */
  qhg_gluon_loop gl = qhg_gluon_loop_init(lat);

  /*
    Stout string used in filenames
   */
  char *stoutstr;
  asprintf(&stoutstr, "sN%do%s", n_stout, flt_str(omega_stout));
  /* 
     Calculate gluon loops
   */
  qhg_calculate_gluon_loop(gl, gf);

  /*
    Write gluon loops
  */
  if(true) {
      t0 = qhg_stop_watch(0);      
      char *fname;
      asprintf(&fname, "%s/gluon_loops_%s.h5", rp.corr_dir, stoutstr);
      char *group;
      asprintf(&group, "/gluon_loop/");

      qhg_write_gluon_loops(fname, gl, group);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 
      free(fname);
      free(group);
  }
	
  if(am_io_proc)
    printf("Done gluon loops in %g sec\n",
	   qhg_stop_watch(gl_timer));
      
  free(stoutstr);

  /* 
     Destroy gauge-fields and -loops
  */
  qhg_gauge_field_finalize(gf);  
  qhg_gauge_field_finalize(gf_stout);  
  qhg_gluon_loop_finalize(gl);  
  qhg_lattice_finalize(lat);
  qhg_comms_finalize(comms);
  return 0;
}
  
