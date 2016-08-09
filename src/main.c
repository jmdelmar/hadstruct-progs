#include <mpi.h>
#include <stdio.h>
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

  mg_state mg_state = mg_init(rp, gf);
  
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
      for(int i=0; i<NS*NC; i++) {
	mg_invert(sol_f[flav][i], src[i], 1e-9, flav == 0 ? plus : minus, &mg_state);
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
    
    /*
      Allocate sequential source and solution
     */
    qhg_spinor_field seq_src[NC*NS], seq_sol[NC*NS];
    for(int i=0; i<NC*NS; i++) {
      seq_src[i] = qhg_spinor_field_init(lat, bc);
      seq_sol[i] = qhg_spinor_field_init(lat, bc);
    }
    
    /*
      Loop over sequential sink parameters
    */
    for(int isnk=0; isnk<rp.spos[isrc].nsinks; isnk++) {
      double sink_timer = qhg_stop_watch(0);
      qhg_thrp_nn_sink_params thrp_snk = rp.spos[isrc].sinks[isnk];
      for(int flav=0; flav<NF; flav++) {
    	/*
    	   Jump to reading backward props if available
    	*/
    	if(read_bwd_props)
    	  goto BREAD;

    	t0 = qhg_stop_watch(0);
    	switch(flav) {
    	case up:
    	  qhg_nn_sequential_sink_u(seq_src, sol_sm_u, sol_sm_d, sco[0], thrp_snk);
    	  break;
    	case dn:
    	  qhg_nn_sequential_sink_d(seq_src, sol_sm_u, sco[0], thrp_snk);
    	  break;
    	}
    	if(am_io_proc)
    	  printf("Done sequential source in %g sec\n", qhg_stop_watch(t0));

    	/*
    	  Smear the sequential source
    	*/
    	t0 = qhg_stop_watch(0);
    	if(am_io_proc)
    	  printf("Smearing the sequential sink, isnk = %2d, Proj = %s, sink-source = %2d, flav = %s\n",
    		 isnk, proj_to_str(thrp_snk.proj), thrp_snk.dt, flav_str[flav]);
    	for(int sp=0; sp<NS; sp++)
    	  for(int co=0; co<NC; co++) {
    	    if(am_io_proc)
    	      printf("\tcol=%d, spin=%d\n", co, sp);
    	    qhg_gauss_smear(seq_src[CS(sp,co)], seq_src[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
    	  }
    	if(am_io_proc)
    	  printf("Done smearing in %g sec\n", qhg_stop_watch(t0));
	
    	t0 = qhg_stop_watch(0);
	for(int i=0; i<NS*NC; i++) {
	  mg_invert(seq_sol[i], seq_src[i], 1e-9, flav == 0 ? minus : plus, &mg_state);
	}
	
    	if(am_io_proc)
    	  printf("Done %s sequential inversion in %g sec\n", flav_str[flav], qhg_stop_watch(t0));
	
    	/*
    	  Write the propagators if selected
    	*/
    	if(write_bwd_props) {
    	  t0 = qhg_stop_watch(0);
    	  char *propname;
    	  asprintf(&propname, "%s/backprop_%s_%s_dt%02d.%s", rp.prop_dir, srcstr, proj_to_str(thrp_snk.proj),
    		   thrp_snk.dt, flav_str[flav]);
    	  qhg_write_spinors(propname, NC*NS, seq_sol);
    	  if(am_io_proc)
    	    printf("Wrote %s in %g sec\n", propname, qhg_stop_watch(t0));
    	  free(propname);
    	}

      BREAD: if(read_bwd_props) {
    	  t0 = qhg_stop_watch(0);
    	  char *propname;
    	  asprintf(&propname, "%s/backprop_%s_%s_dt%02d.%s", rp.prop_dir, srcstr, proj_to_str(thrp_snk.proj),
    		   thrp_snk.dt, flav_str[flav]);
    	  qhg_read_spinors(seq_sol, NC*NS, propname);
    	  if(am_io_proc)
    	    printf("Read %s in %g sec\n", propname, qhg_stop_watch(t0));
    	  free(propname);
    	}
	
    	/*
    	  backprop = (\gamma_5 backprop)^\dagger
    	 */
    	qhg_prop_field_g5_G(seq_sol);
    	qhg_prop_field_Gdag(seq_sol);
        
    	qhg_spinor_field *fwd[2] = {sol_u, sol_d};
	
    	/*
    	  Three-point function. Needs gauge-field for derivative
    	  operators.
    	 */
    	t0 = qhg_stop_watch(0);
    	qhg_thrp_correlator thrp;
    	thrp.corr = qhg_nn_thrp(fwd[flav], seq_sol, gf, sco, thrp_snk);
    	thrp.flav = flav;
    	thrp.dt = thrp_snk.dt;
    	thrp.proj = thrp_snk.proj;

    	if(am_io_proc)
    	  printf("Done three-point correlator in %g sec\n", qhg_stop_watch(t0));
	
    	/*
    	  Write three-point function
    	 */
    	if(write_thrp_pos_space) {
    	  t0 = qhg_stop_watch(0);
    	  char *fname;
    	  asprintf(&fname, "%s/thrp_%s_%s_%s_%s_dt%02d.%s.h5",
    		   rp.corr_dir, srcstr, smrstr, apestr, proj_to_str(thrp_snk.proj),
    		   thrp_snk.dt, flav_str[flav]);
    	  qhg_write_nn_thrp(fname, thrp);
    	  if(am_io_proc)
    	    printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0));
    	  free(fname);
    	}
    	qhg_correlator_finalize(thrp.corr);
	
      }
      if(am_io_proc)
    	printf("Done sink: proj = %s, sink-source = %2d, in %g sec\n",
    	       proj_to_str(thrp_snk.proj), thrp_snk.dt, qhg_stop_watch(sink_timer));
    }
      
    /*
      Free sequential source and solution
    */
    for(int i=0; i<NC*NS; i++) {
      qhg_spinor_field_finalize(seq_src[i]);
      qhg_spinor_field_finalize(seq_sol[i]);
    }
    free(srcstr);
    if(am_io_proc)
      printf("Done source, coords (t, x, y, z) = (%d, %d, %d, %d), in %g sec\n",
	     sco[0], sco[1], sco[2], sco[3], qhg_stop_watch(source_timer));
  }
  free(smrstr);

  mg_finalize();
  
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
  
