#include <mpi.h>
#include <tmLQCD.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <qhg.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif


static int write_nuclns_pos_space = 1;
static int write_mesons_pos_space = 1;
static int write_thrp_pos_space = 1;

static int read_fwd_props = 1;
static int write_fwd_props = 0;

static int read_bwd_props = 1;
static int write_bwd_props = 0;

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

#define NF 2
#define NSRC (16)
#define NSNK (12)

char flav_str[NF][3] = {"up","dn"};

static void
invert(qhg_spinor_field psi, qhg_spinor_field eta, int op_id, enum qhg_fermion_bc_time bc)
{
  int write_prop = 0;
  qhg_spinor_field eta_theta = qhg_spinor_field_init(eta.lat, bc);
  qhg_spinor_field psi_theta = qhg_spinor_field_init(psi.lat, bc);
  double bc_angle_t;
  
  switch(bc) {
  case PERIODIC:
    bc_angle_t = 0.0;
    break;
  case ANTIPERIODIC:
    bc_angle_t = 1.0;
    break;
  }
  
  qhg_spinor_twist_t_bc(eta_theta, eta, -bc_angle_t);
  
  tmLQCD_invert((double *) psi_theta.field,
		(double *) eta_theta.field,
		op_id, write_prop);
  
  qhg_spinor_twist_t_bc(psi, psi_theta, bc_angle_t);
  
  qhg_spinor_field_finalize(eta_theta);
  qhg_spinor_field_finalize(psi_theta);
  return;
}

  
void
usage(char *argv[])
{
  fprintf(stderr, " Usage: %s CONF_NUMB\n", argv[0]);
  return;
}

int
main(int argc, char *argv[])
{   
  if(argc != 2) {
    usage(argv);
    exit(1);
  }
  char *e;
  int config = (int)strtoul(argv[1], &e, 10);
  if(*e != '\0') {
    usage(argv);
    exit(2);
  }
  
  int dims[ND] = {96, 48, 48, 48}; // t,x,y,z
  double bc_angle_t = 1.0; // Match tmLQCD's parameter. 
  int n_ape = 50;
  double alpha_ape = 0.5;
  int n_gauss = 50;
  double alpha_gauss = 4.0;
  int source_coords[NSRC][ND];
  enum qhg_fermion_bc_time bc = ANTIPERIODIC; // Also set this in tmLQCD's input
  qhg_thrp_nn_sink_params thrp_snk[NSNK];
  
  {
    char *fname_sinks;
    asprintf(&fname_sinks, "sinks.%04d.list", config);  
    FILE *fp = qhg_fopen(fname_sinks, "r");
    for(int i=0; i<NSNK; i++) {
      int dt;
      char proj[256];
      fscanf(fp, "%d %s", &dt, proj);
      thrp_snk[i].dt = dt;
      thrp_snk[i].proj = str_to_proj(proj);
    }
    fclose(fp);
  }

  {
    char *fname_sourcepos;
    asprintf(&fname_sourcepos, "sources.%04d.list", config);  
    FILE *fp = qhg_fopen(fname_sourcepos, "r");
    for(int i=0; i<NSRC; i++) {
      int *s = &source_coords[i][0];
      fscanf(fp, "%d %d %d %d", &s[0], &s[1], &s[2], &s[3]);
    }
    fclose(fp);
  }

  char prop_dir[256];
  {
    char *fname_prop_dir;
    asprintf(&fname_prop_dir, "prop.%04d.dir", config);  
    FILE *fp = qhg_fopen(fname_prop_dir, "r");
    fscanf(fp, "%s", prop_dir);
    fclose(fp);
  }
  
  char corr_dir[256];
  {
    char *fname_corr_dir;
    asprintf(&fname_corr_dir, "corr.%04d.dir", config);  
    FILE *fp = qhg_fopen(fname_corr_dir, "r");
    fscanf(fp, "%s", corr_dir);
    fclose(fp);
  }

  qhg_lattice *lat = qhg_lattice_init(dims);
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  qhg_gauge_field gf = qhg_gauge_field_init(lat);  
  
  /*
    Use tmLQCD to read config 
  */
  tmLQCD_read_gauge(config);
  _Complex double *gptr;
  tmLQCD_get_gauge_field_pointer((double **)&gptr);  

  /*
    import to qhg
  */
  qhg_import_gauge_field(gf, gptr);
  
  /*
    plaquette 
  */
  double p = qhg_plaquette(gf);
  if(am_io_proc)
    printf("Plaquette = %10.8f\n", p);

  /*
    APE smear in 3-dimensions
   */
  double t0 = qhg_stop_watch(0);
  qhg_gauge_field gf_ape = qhg_gauge_field_init(lat);
  qhg_ape_smear_3d(gf_ape, gf, alpha_ape, n_ape);
  if(am_io_proc)
    printf("3D APE smear in %g sec\n", qhg_stop_watch(t0));
  
  /*
    plaquette of smeared gauge-field
  */
  double p_ape = qhg_plaquette(gf_ape);
  if(am_io_proc)
    printf("3D APE plaquette = %10.8f\n", p_ape);
  
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
  for(int isrc=0; isrc<NSRC; isrc++) {
    double source_timer = qhg_stop_watch(0);
    /* 
       This source position coordinates 
    */
    int *sco = &source_coords[isrc][0];
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
    int op[] = {0, 1};
    qhg_spinor_field *sol_f[] = {sol_u, sol_d};
    for(int i=0; i<NS*NC; i++)
      for(int flav=0; flav<NF; flav++)
	invert(sol_f[flav][i], src[i], op[flav], bc_angle_t);
    if(am_io_proc)
      printf("Done up & dn inversion in %g sec\n", qhg_stop_watch(t0));  
    
    /*
      Write the propagators if selected
    */
    if(write_fwd_props) {
      char *propname;

      t0 = qhg_stop_watch(0);
      asprintf(&propname, "%s/prop_%s.up", prop_dir, srcstr);      
      qhg_write_spinors(propname, NC*NS, sol_u);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", propname, qhg_stop_watch(t0));
      free(propname);
      
      t0 = qhg_stop_watch(0);
      asprintf(&propname, "%s/prop_%s.dn", prop_dir, srcstr);      
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
      asprintf(&propname, "%s/prop_%s.up", prop_dir, srcstr);      
      qhg_read_spinors(sol_u, NC*NS, propname);
      if(am_io_proc)
	printf("Read %s in %g sec\n", propname, qhg_stop_watch(t0));
      free(propname);
      
      t0 = qhg_stop_watch(0);      
      asprintf(&propname, "%s/prop_%s.dn", prop_dir, srcstr);      
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
      asprintf(&fname, "%s/mesons_%s_%s_%s.h5", corr_dir, srcstr, smrstr, apestr);
      qhg_write_mesons(fname, mesons);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 
      free(fname);
    }
    qhg_correlator_finalize(mesons);
    
    if(write_nuclns_pos_space) {
      t0 = qhg_stop_watch(0);
      char *fname;
      asprintf(&fname, "%s/nucleons_%s_%s_%s.h5", corr_dir, srcstr, smrstr, apestr);
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
    for(int isnk=0; isnk<NSNK; isnk++) {
      double sink_timer = qhg_stop_watch(0);
      for(int flav=0; flav<NF; flav++) {
    	/*
    	   Jump to reading backward props if available
    	*/
    	if(read_bwd_props)
    	  goto BREAD;

    	t0 = qhg_stop_watch(0);
    	switch(flav) {
    	case up:
    	  qhg_nn_sequential_sink_u(seq_src, sol_sm_u, sol_sm_d, sco[0], thrp_snk[isnk]);
    	  break;
    	case dn:
    	  qhg_nn_sequential_sink_d(seq_src, sol_sm_u, sco[0], thrp_snk[isnk]);
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
    		 isnk, proj_to_str(thrp_snk[isnk].proj), thrp_snk[isnk].dt, flav_str[flav]);
    	for(int sp=0; sp<NS; sp++)
    	  for(int co=0; co<NC; co++) {
    	    if(am_io_proc)
    	      printf("\tcol=%d, spin=%d\n", co, sp);
    	    qhg_gauss_smear(seq_src[CS(sp,co)], seq_src[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
    	  }
    	if(am_io_proc)
    	  printf("Done smearing in %g sec\n", qhg_stop_watch(t0));
	
    	t0 = qhg_stop_watch(0);
    	int op_id[2] = {1, 0};
    	for(int i=0; i<NS*NC; i++)
	  invert(seq_sol[i], seq_src[i], op_id[flav], bc_angle_t);
	
    	if(am_io_proc)
    	  printf("Done %s sequential inversion in %g sec\n", flav_str[flav], qhg_stop_watch(t0));
	
    	/*
    	  Write the propagators if selected
    	*/
    	if(write_bwd_props) {
    	  t0 = qhg_stop_watch(0);
    	  char *propname;
    	  asprintf(&propname, "%s/backprop_%s_%s_dt%02d.%s", prop_dir, srcstr, proj_to_str(thrp_snk[isnk].proj),
    		   thrp_snk[isnk].dt, flav_str[flav]);
    	  qhg_write_spinors(propname, NC*NS, seq_sol);
    	  if(am_io_proc)
    	    printf("Wrote %s in %g sec\n", propname, qhg_stop_watch(t0));
    	  free(propname);
    	}

      BREAD: if(read_bwd_props) {
    	  t0 = qhg_stop_watch(0);
    	  char *propname;
    	  asprintf(&propname, "%s/backprop_%s_%s_dt%02d.%s", prop_dir, srcstr, proj_to_str(thrp_snk[isnk].proj),
    		   thrp_snk[isnk].dt, flav_str[flav]);
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
    	thrp.corr = qhg_nn_thrp(fwd[flav], seq_sol, gf, sco, thrp_snk[isnk]);
    	thrp.flav = flav;
    	thrp.dt = thrp_snk[isnk].dt;
    	thrp.proj = thrp_snk[isnk].proj;

    	if(am_io_proc)
    	  printf("Done three-point correlator in %g sec\n", qhg_stop_watch(t0));
	
    	/*
    	  Write three-point function
    	 */
    	if(write_thrp_pos_space) {
    	  t0 = qhg_stop_watch(0);
    	  char *fname;
    	  asprintf(&fname, "%s/thrp_%s_%s_%s_%s_dt%02d.%s.h5",
    		   corr_dir, srcstr, smrstr, apestr, proj_to_str(thrp_snk[isnk].proj),
    		   thrp_snk[isnk].dt, flav_str[flav]);
    	  qhg_write_nn_thrp(fname, thrp);
    	  if(am_io_proc)
    	    printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0));
    	  free(fname);
    	}
    	qhg_correlator_finalize(thrp.corr);
	
      }
      if(am_io_proc)
    	printf("Done sink: proj = %s, sink-source = %2d, in %g sec\n",
    	       proj_to_str(thrp_snk[isnk].proj), thrp_snk[isnk].dt, qhg_stop_watch(sink_timer));
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
  return 0;
}
  
