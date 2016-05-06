#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <qhg.h>
#include <parser.h>
#include <mg4qcd_interface.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

static int write_nuclns_pos_space = 1;
static int write_mesons_pos_space = 1;
static int write_thrp_pos_space = 1;

static int read_fwd_props = 0;
static int write_fwd_props = 0;

static int read_bwd_props = 0;
static int write_bwd_props = 0;

#define NF 2
char flav_str[NF][3] = {"up","dn"};

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

  struct run_params rp = parse_input(argv[1]);
  int n_ape = 50;
  double alpha_ape = 0.5;
  int n_gauss = 50;
  double alpha_gauss = 4.0;
  enum qhg_fermion_bc_time bc = ANTIPERIODIC; // Also set this in tmLQCD's input
  qhg_comms *comms = qhg_comms_init(rp.procs);  
  qhg_lattice *lat = qhg_lattice_init(rp.dims, comms);
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  qhg_gauge_field gf = qhg_gauge_field_init(lat);  

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

  mg4qcd_state mg_state = mg4qcd_init(rp, gf);
  double mu = mg4qcd_get_mu(mg_state);
  int n_levs = mg4qcd_get_numb_levels(mg_state);
  double coarse_mu[n_levs];
  for(int i=0; i<n_levs; i++)
    coarse_mu[i] = mg4qcd_get_coarse_mu(mg_state, i);
  
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
    Source and solution spinor field
   */
  qhg_spinor_field src[NS*NC];
  qhg_spinor_field sol_u[NS*NC], sol_d[NS*NC];
  qhg_spinor_field sol_sm_u[NS*NC], sol_sm_d[NS*NC];
  for(int sp=0; sp<NS; sp++)
    for(int co=0; co<NC; co++) {
      sol_u[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_d[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_sm_u[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_sm_d[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
    }

  /*
    APE string used in filenames
   */
  char *apestr;
  asprintf(&apestr, "aN%da%s", n_ape, flt_str(alpha_ape));

  /*
    Source smearing string used in filenames
  */
  char *smrstr;
  asprintf(&smrstr, "gN%da%s", n_gauss, flt_str(alpha_gauss));
  
  /*
    Loop over source positions
  */  
  for(int isrc=0; isrc<rp.nsp; isrc++) {
    double source_timer = qhg_stop_watch(0);
    /* 
       This source position coordinates 
    */
    int *sco = &rp.spos[isrc].coords[0];
    if(am_io_proc)
      printf("Source coords (t, x, y, z) = (%d, %d, %d, %d)\n",
	     sco[0], sco[1], sco[2], sco[3]);

    /*
      Source position string to be used in filenames
    */
    char *srcstr;
    asprintf(&srcstr, "sx%02dsy%02dsz%02dst%02d", sco[1], sco[2], sco[3], sco[0]);

    /*
      Jump to reading forward props if available
     */
    if(read_fwd_props)
      goto FREAD;
    
    /*
      Smear the source for this source position
    */
    t0 = qhg_stop_watch(0);
    if(am_io_proc)
      printf("Smearing the source\n");  
    for(int sp=0; sp<NS; sp++)
      for(int co=0; co<NC; co++) {
	if(am_io_proc)
	  printf("\tcol=%d, spin=%d\n", co, sp);  
	qhg_spinor_field aux = qhg_spinor_field_init(lat, bc);
	qhg_point_spinor_field(aux, sco, sp, co);
	src[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
	qhg_gauss_smear(src[CS(sp,co)], aux, gf_ape, alpha_gauss, n_gauss);
	qhg_spinor_field_finalize(aux);
      }
    if(am_io_proc)
      printf("Done smearing in %g sec\n", qhg_stop_watch(t0));  
    
    /*
      Invert for up- and down-quark in the same color-spin loop
    */
    t0 = qhg_stop_watch(0);
    qhg_spinor_field *sol_f[] = {sol_u, sol_d};
    for(int flav=0; flav<NF; flav++) {
      int s = flav == 0 ? 1 : -1;
      double x_mu = s*mu;
      double x_coarse_mu[n_levs];
      for(int i=0; i<n_levs; i++)
	x_coarse_mu[i] = s*coarse_mu[i];

      mg4qcd_change_mu(&mg_state, x_mu, x_coarse_mu);
      for(int i=0; i<NS*NC; i++) {
	mg4qcd_invert(sol_f[flav][i], src[i], 1e-9, mg_state);
      }
    }
    if(am_io_proc)
      printf("Done up & dn inversion in %g sec\n", qhg_stop_watch(t0));  
    
    /*
      Write the propagators if selected
    */
    if(write_fwd_props) {
      char *propname;

      t0 = qhg_stop_watch(0);
      asprintf(&propname, "%s/prop_%s.up", rp.prop_dir, srcstr);      
      qhg_write_spinors(propname, NC*NS, sol_u);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", propname, qhg_stop_watch(t0));
      free(propname);
      
      t0 = qhg_stop_watch(0);
      asprintf(&propname, "%s/prop_%s.dn", rp.prop_dir, srcstr);      
      qhg_write_spinors(propname, NC*NS, sol_d);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", propname, qhg_stop_watch(t0));
      free(propname);
    }
        
    /*
      We don't need the source any more
    */
    for(int i=0; i<NS*NC; i++)
      qhg_spinor_field_finalize(src[i]);

  FREAD: if(read_fwd_props) {
      char *propname;

      t0 = qhg_stop_watch(0);      
      asprintf(&propname, "%s/prop_%s.up", rp.prop_dir, srcstr);      
      qhg_read_spinors(sol_u, NC*NS, propname);
      if(am_io_proc)
	printf("Read %s in %g sec\n", propname, qhg_stop_watch(t0));
      free(propname);
      
      t0 = qhg_stop_watch(0);      
      asprintf(&propname, "%s/prop_%s.dn", rp.prop_dir, srcstr);      
      qhg_read_spinors(sol_d, NC*NS, propname);
      if(am_io_proc)
	printf("Read %s in %g sec\n", propname, qhg_stop_watch(t0));
      free(propname);            
    }
        
    /*
      Smear the propagators sink-side
    */
    t0 = qhg_stop_watch(0);
    if(am_io_proc)
      printf("Smearing the propagator sink\n");  
    for(int sp=0; sp<NS; sp++)
      for(int co=0; co<NC; co++) {
	if(am_io_proc)
	  printf("\tcol=%d, spin=%d\n", co, sp);  
	qhg_gauss_smear(sol_sm_u[CS(sp,co)], sol_u[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
	qhg_gauss_smear(sol_sm_d[CS(sp,co)], sol_d[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
      }
    if(am_io_proc)
      printf("Done smearing in %g sec\n", qhg_stop_watch(t0));  
    
    /*
      Smeared nucleon and meson correlators and fourier transform
    */
    t0 = qhg_stop_watch(0);
    qhg_correlator mesons = qhg_mesons(sol_sm_u, sol_sm_d, sco);
    if(am_io_proc)
      printf("Done meson correlator in %g sec\n", qhg_stop_watch(t0)); 

    t0 = qhg_stop_watch(0);    
    qhg_correlator nucleons = qhg_nucleons(sol_sm_u, sol_sm_d, sco);
    if(am_io_proc)
      printf("Done nucleon correlator in %g sec\n", qhg_stop_watch(t0)); 
    
    /*
      Write the correlators
    */
    if(write_mesons_pos_space) {
      t0 = qhg_stop_watch(0);      
      char *fname;
      asprintf(&fname, "%s/mesons_%s_%s_%s.h5", rp.corr_dir, srcstr, smrstr, apestr);
      qhg_write_mesons(fname, mesons);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 
      free(fname);
    }
    qhg_correlator_finalize(mesons);
    
    if(write_nuclns_pos_space) {
      t0 = qhg_stop_watch(0);
      char *fname;
      asprintf(&fname, "%s/nucleons_%s_%s_%s.h5", rp.corr_dir, srcstr, smrstr, apestr);
      qhg_write_nucleons(fname, nucleons);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 
      free(fname);
    }
    qhg_correlator_finalize(nucleons);
    free(srcstr);
    if(am_io_proc)
      printf("Done source, coords (t, x, y, z) = (%d, %d, %d, %d), in %g sec\n",
	     sco[0], sco[1], sco[2], sco[3], qhg_stop_watch(source_timer));
  }
  free(smrstr);

  mg4qcd_finalize();
  
  /* 
     Destroy spinor- and gauge-fields
  */
  for(int i=0; i<NS*NC; i++) {
    qhg_spinor_field_finalize(sol_u[i]);
    qhg_spinor_field_finalize(sol_d[i]);
    qhg_spinor_field_finalize(sol_sm_u[i]);
    qhg_spinor_field_finalize(sol_sm_d[i]);
  }
  
  qhg_gauge_field_finalize(gf);  
  qhg_gauge_field_finalize(gf_ape);  
  qhg_lattice_finalize(lat);
  qhg_comms_finalize(comms);
  return 0;
}
  
