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

#define NF 3
char part_str[NF][7] = {"pion","kaon","kaon"};
char flav_str[NF][7] = {"up","up","strange"};

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

  int n_ape = rp.smearing.n_ape;
  double alpha_ape = rp.smearing.alpha_ape;
  int n_gauss_l = rp.smearing.n_gauss_l;
  double alpha_gauss_l = rp.smearing.alpha_gauss_l;
  int n_gauss_s = rp.smearing.n_gauss_s;
  double alpha_gauss_s = rp.smearing.alpha_gauss_s;

  enum qhg_fermion_bc_time bc = ANTIPERIODIC; // Also set this in tmLQCD's input
  qhg_comms *comms = qhg_comms_init(rp.procs);  
  qhg_lattice *lat = qhg_lattice_init(rp.dims, comms);
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  qhg_gauge_field gf = qhg_gauge_field_init(lat);  

  if(am_io_proc){
    printf("List of momenta:\n");
    for(int imom=0; imom<rp.spos[0].nmoms; imom++) {
      int *mom_vec = rp.spos[0].mom_vecs[imom];
      printf("mom = %d %d %d\n", mom_vec[0], mom_vec[1], mom_vec[2]);
    }
  }

  // Read config 

  qhg_read_gauge_field_ildg(gf, rp.config);
  
  // Plaquette 

  double p = qhg_plaquette(gf);
  if(am_io_proc)
    printf("Plaquette = %12.10f\n", p);

  mg_state mg_state = mg_init(rp, gf);

  // APE smear in 3-dimensions

  double t0 = qhg_stop_watch(0);
  qhg_gauge_field gf_ape = qhg_gauge_field_init(lat);
  qhg_ape_smear_3d(gf_ape, gf, alpha_ape, n_ape);
  if(am_io_proc)
    printf("3D APE smear in %g sec\n", qhg_stop_watch(t0));
  
  // Plaquette of smeared gauge-field

  double p_ape = qhg_plaquette(gf_ape);
  if(am_io_proc)
    printf("3D APE plaquette = %12.10f\n", p_ape);

  // Source and solution spinor field

  qhg_spinor_field src_l[NS*NC], src_s[NS*NC];
  qhg_spinor_field sol_u[NS*NC], sol_s[NS*NC];
  qhg_spinor_field sol_sm_u[NS*NC], sol_sm_s[NS*NC];
  for(int sp=0; sp<NS; sp++)
    for(int co=0; co<NC; co++) {
      sol_u[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_s[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_sm_u[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_sm_s[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
    }

  // APE string used in filenames

  char *apestr;
  asprintf(&apestr, "aN%da%s", n_ape, flt_str(alpha_ape));

  // Source smearing string used in filenames

  char *smrstr_l;
  asprintf(&smrstr_l, "gN%da%s", n_gauss_l, flt_str(alpha_gauss_l));
  char *smrstr_s;
  asprintf(&smrstr_s, "gN%da%s", n_gauss_s, flt_str(alpha_gauss_s));
  char *smrstr_f[] = {smrstr_l, smrstr_s, smrstr_s};  

  // Loop over source positions

  for(int isrc=0; isrc<rp.nsp; isrc++) {
    double source_timer = qhg_stop_watch(0);

    // This source position coordinates 

    int *sco = &rp.spos[isrc].coords[0];
    if(am_io_proc)
      printf("Source coords (t, x, y, z) = (%d, %d, %d, %d)\n",
	     sco[0], sco[1], sco[2], sco[3]);

    // Source position string to be used in filenames

    char *srcstr;
    char *fmt;
    asprintf(&fmt, "sx%%0%ddsy%%0%ddsz%%0%ddst%%0%dd",
	     (int)log10(rp.dims[1])+1,
	     (int)log10(rp.dims[2])+1,
	     (int)log10(rp.dims[3])+1,
	     (int)log10(rp.dims[0])+1);
    asprintf(&srcstr, fmt, sco[1], sco[2], sco[3], sco[0]);

    // Smear the source for this source position

    t0 = qhg_stop_watch(0);
    qhg_spinor_field *src_f[] = {src_l, src_s};
    int n_gauss_f[] = {n_gauss_l, n_gauss_s, n_gauss_s};
    double alpha_gauss_f[] = {alpha_gauss_l, alpha_gauss_s, alpha_gauss_s};

    if(am_io_proc)
      printf("Smearing the source\n");  
    for(int sp=0; sp<NS; sp++)
      for(int co=0; co<NC; co++) {
	if(am_io_proc)
	  printf("\tcol=%d, spin=%d\n", co, sp);  
	qhg_spinor_field aux = qhg_spinor_field_init(lat, bc);
	qhg_point_spinor_field(aux, sco, sp, co);
	for(int flav=0; flav<NF-1; flav++) {
	  src_f[flav][CS(sp,co)] = qhg_spinor_field_init(lat, bc);
	  qhg_gauss_smear(src_f[flav][CS(sp,co)], aux, gf_ape, alpha_gauss_f[flav], n_gauss_f[flav]);
	}
	qhg_spinor_field_finalize(aux);
      }
    if(am_io_proc)
      printf("Done smearing in %g sec\n", qhg_stop_watch(t0));  
    
    // Invert for up- and strange-quark in the same color-spin loop

    t0 = qhg_stop_watch(0);
    qhg_spinor_field *sol_f[] = {sol_u, sol_s};

    for(int flav=0; flav<NF-1; flav++) {
      mg_state.params.mu = flav == 0 ? rp.act.mu_l : rp.act.mu_s;
      DDalphaAMG_update_parameters(&mg_state.params, &mg_state.status);
      for(int i=0; i<NS*NC; i++) {
	mg_invert(sol_f[flav][i], src_f[flav][i], 1e-9, flav == 0 ? minus : plus, &mg_state);
      }
    }
    if(am_io_proc)
      printf("Done up & strange inversion in %g sec\n", qhg_stop_watch(t0));  
    

    // We don't need the sources any more

    for(int i=0; i<NS*NC; i++) {
      qhg_spinor_field_finalize(src_l[i]);
      qhg_spinor_field_finalize(src_s[i]);
    }

    // Smear the propagators sink-side

    t0 = qhg_stop_watch(0);
    if(am_io_proc)
      printf("Smearing the propagator sink\n");  
    for(int sp=0; sp<NS; sp++)
      for(int co=0; co<NC; co++) {
	if(am_io_proc)
	  printf("\tcol=%d, spin=%d\n", co, sp);  
	qhg_gauss_smear(sol_sm_u[CS(sp,co)], sol_u[CS(sp,co)], gf_ape, alpha_gauss_l, n_gauss_l);
	qhg_gauss_smear(sol_sm_s[CS(sp,co)], sol_s[CS(sp,co)], gf_ape, alpha_gauss_s, n_gauss_s);
      }
    if(am_io_proc)
      printf("Done smearing in %g sec\n", qhg_stop_watch(t0));  
    
    // Smeared meson correlators and fourier transform

    t0 = qhg_stop_watch(0);
    qhg_correlator corr_pion = qhg_mesons_pseudoscalar(sol_sm_u, sol_sm_u, sco);
    qhg_correlator_shift(corr_pion, corr_pion.origin);
    qhg_correlator corr_kaon = qhg_mesons_pseudoscalar(sol_sm_u, sol_sm_s, sco);
    qhg_correlator_shift(corr_kaon, corr_kaon.origin);
    if(am_io_proc)
      printf("Done meson correlator in %g sec\n", qhg_stop_watch(t0)); 

    // Write the correlators

    if(true) {
      t0 = qhg_stop_watch(0);      
      char *fname;
      asprintf(&fname, "%s/twop_pion_%s_%s_%s.h5", rp.corr_dir, srcstr, smrstr_l, apestr);
      char *group;
      asprintf(&group, "twop_pion/%s/", srcstr);      
      qhg_write_single_meson(fname, corr_pion, group);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 

      asprintf(&fname, "%s/twop_kaon_%s_%s_%s.h5", rp.corr_dir, srcstr, smrstr_s, apestr);
      asprintf(&group, "twop_kaon/%s/", srcstr);      
      qhg_write_single_meson(fname, corr_kaon, group);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 
      free(fname);
      free(group);
    }
    qhg_correlator_finalize(corr_pion);
    qhg_correlator_finalize(corr_kaon);
        
    // Allocate sequential source and solution

    qhg_spinor_field seq_src[NC*NS], seq_sol[NC*NS];
    for(int i=0; i<NC*NS; i++) {
      seq_src[i] = qhg_spinor_field_init(lat, bc);
      seq_sol[i] = qhg_spinor_field_init(lat, bc);
    }

    // Loop over sequential sink parameters

    for(int isnk=0; isnk<rp.spos[isrc].nsinks; isnk++) {
      
      // NEW: loop over nmoms
      for(int imom=0; imom<rp.spos[isrc].nmoms; imom++) {
	int *mom_vec = rp.spos[isrc].mom_vecs[imom];

	double sink_timer = qhg_stop_watch(0);
	qhg_thrp_nn_sink_params thrp_snk = rp.spos[isrc].sinks[isnk];
	for(int flav=0; flav<NF; flav++) {
	  t0 = qhg_stop_watch(0);
	  switch(flav) {
	  case 0: // pion
	  case 2: // kaon, strange part
	    qhg_mesons_sequential_sink(seq_src, sol_sm_u, sco[0], thrp_snk);
	    break;
	  case 1: // kaon, up part
	    qhg_mesons_sequential_sink(seq_src, sol_sm_s, sco[0], thrp_snk);
	    break;
	  }
	  if(am_io_proc)
	    printf("Done sequential source in %g sec\n", qhg_stop_watch(t0));

	  // NEW: this function add the phase to the vector
	  qhg_phase_sequential_sink(seq_src, seq_src, mom_vec, sco, -1);

	
	  // Smear the sequential source
	  t0 = qhg_stop_watch(0);
	  if(am_io_proc)
	    printf("Smearing the sequential sink, isnk = %2d, sink-source = %2d, flav = %s\n",
		   isnk, thrp_snk.dt, flav_str[flav]);
	  for(int sp=0; sp<NS; sp++)
	    for(int co=0; co<NC; co++) {
	      if(am_io_proc)
		printf("\tcol=%d, spin=%d\n", co, sp);
	      qhg_gauss_smear(seq_src[CS(sp,co)], seq_src[CS(sp,co)], gf_ape, alpha_gauss_f[flav], n_gauss_f[flav]);
	    }
	  if(am_io_proc)
	    printf("Done smearing in %g sec\n", qhg_stop_watch(t0));

	  // Set mu to correct flavor

	  mg_state.params.mu = flav == 2 ? rp.act.mu_s : rp.act.mu_l;
	  DDalphaAMG_update_parameters(&mg_state.params, &mg_state.status);
	
	  t0 = qhg_stop_watch(0);
	  for(int i=0; i<NS*NC; i++) {
	    mg_invert(seq_sol[i], seq_src[i], 1e-9, plus, &mg_state);
	    //mg_invert(seq_sol[i], seq_src[i], 1e-9, flav == 2 ? minus : plus, &mg_state);
	  }
	
	  if(am_io_proc)
	    printf("Done %s sequential inversion in %g sec\n", flav_str[flav], qhg_stop_watch(t0));
		
	  // backprop = (\gamma_5 backprop)^\dagger

	  qhg_prop_field_g5_G(seq_sol);
	  qhg_prop_field_Gdag(seq_sol);
        
	  qhg_spinor_field *fwd[3] = {sol_u, sol_u, sol_s};
	
	  // Three-point function. Needs gauge-field for derivative
	  // operators.

	  t0 = qhg_stop_watch(0);
	  qhg_thrp_correlator thrp;
	  thrp.corr = qhg_nn_thrp(fwd[flav], seq_sol, gf, sco, thrp_snk);
	  switch( flav ) {
	  case 0: // pion
	    thrp.flav = up;
	    break;
	  case 1: // kaon, up part
	    thrp.flav = up;
	    break;
	  case 2: // kaon, strange part
	    thrp.flav = strange;
	    break;
	  }

	  thrp.dt = thrp_snk.dt;

	  // If we are on the kaon strange part, we have calculated the 
	  // three-point functions for the down part of K^0 so we need 
	  // to take the complex conjugate to get the three-point functions
	  // for the strange part of K^+

	  if( flav == 2 )
	    qhg_conjugate_thrp( thrp.corr, 1 );

	  if(am_io_proc)
	    printf("Done three-point correlator in %g sec\n", qhg_stop_watch(t0));
	
	  // Write three-point function

	  if(true) {
	    t0 = qhg_stop_watch(0);
	    char *fname;
	    // NEW: changed the name for having the momentum
	    asprintf(&fname, "%s/thrp_%s_%s_%s_%s_dt%02d_mom_%+d_%+d_%+d.%s.h5",
		     rp.corr_dir, part_str[flav], srcstr, smrstr_f[flav], apestr, 
		     thrp_snk.dt, mom_vec[0], mom_vec[1], mom_vec[2], flav_str[flav]);
	    char *group;
	    asprintf(&group, "thrp/%s/dt%02d/%s",
		     srcstr, thrp_snk.dt, flav_str[flav]);
	    qhg_correlator_shift(thrp.corr, thrp.corr.origin);
	    qhg_write_mesons_thrp(fname, thrp, group);
	    if(am_io_proc)
	      printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0));
	    free(fname);
	    free(group);
	  }
	  qhg_correlator_finalize(thrp.corr);

	  /*
	   * NEW: this function does the second derivative. 
	   * It's just for crosscheck,
	   * At the moment let's keep it but after can be removed
	   *
	   * Three-point double derivative function. 
	   */
	  if(false) {
	    t0 = qhg_stop_watch(0);
	    thrp.corr = qhg_nn_thrp_der2(fwd[flav], seq_sol, gf, sco, thrp_snk);
	    thrp.flav = flav;
	    thrp.dt = thrp_snk.dt;
	    thrp.proj = thrp_snk.proj;
          
	    if(am_io_proc)
	      printf("Done three-point correlator in %g sec\n", qhg_stop_watch(t0));

	    t0 = qhg_stop_watch(0);
	    char *fname;
	    // NEW: changed the name for having the momentum
	    asprintf(&fname, "%s/thrp_der2_%s_%s_%s_%s_dt%02d_mom_%+d_%+d_%+d.%s.h5",
		     rp.corr_dir, part_str[flav], srcstr, smrstr_f[flav], apestr, 
		     thrp_snk.dt, mom_vec[0], mom_vec[1], mom_vec[2], flav_str[flav]);
	    char *group;
	    asprintf(&group, "thrp/%s/dt%02d/%s",
		     srcstr, thrp_snk.dt, flav_str[flav]);
	    qhg_correlator_shift(thrp.corr, thrp.corr.origin);
	    qhg_write_nn_thrp_der2(fname, thrp, group);
	    if(am_io_proc)
	      printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0));
	    free(fname);
	    free(group);

	    qhg_correlator_finalize(thrp.corr);
	  }

	  /*
	   * NEW: this function can be used to skip some combination of derivative.
	   * At the moment returns false by default,
	   * which means all the directions will be done
	   */
	  bool to_skip(int dirs[8]) {
	    if(am_io_proc) {
	      char dir[4] = {'t', 'x', 'y', 'z'};
	      printf("Running over dirs: ");
	      for(int i=0; i<8; i++) {
		if(dirs[i] > 0) 
		  printf("%+d%c ",(1-2*(i/4))*dirs[i], dir[i%4]);
	      }
	      printf("\n");
	    }
	    return false;
	  }
	
	  /*
	   * NEW: Here is the loop over all the derivative. 
	   * We do all of them so we can crosscheck, 
	   * but after only 2nd and 3rd will be done here.
	   */
	  for(int der_order = 1; der_order <=3; der_order++) {
	    t0 = qhg_stop_watch(0);
	    qhg_der_correlator corr = qhg_nn_thrp_der(fwd[flav], seq_sol, gf, sco, thrp_snk, der_order, to_skip, true, true);

	    switch( flav ) {
	    case 0: // pion
	      corr.flav = up;
	      break;
	    case 1: // kaon, up part
	      corr.flav = up;
	      break;
	    case 2: // kaon, strange part
	      corr.flav = strange;
	      break;
	    }
	    corr.dt = thrp_snk.dt;

	    if( der_order == 2 || der_order == 3) {
	      qhg_der_correlator corr_avg = qhg_avg_der_combos(corr, mom_vec); 
	      corr = qhg_averaged_der_correlator_copy(corr_avg);
	      qhg_der_correlator_finalize(corr_avg);
	    }

	    if(am_io_proc)
	      printf("Done three-point %d derivative correlator in %g sec\n", der_order, qhg_stop_watch(t0));
	
	    /*
	      Write three-point function
	    */
	    if(true) {
	      t0 = qhg_stop_watch(0);
	      char *fname;
	      // NEW: changed the name for having the momentum
	      asprintf(&fname, "%s/thrp_der%d_%s_%s_%s_%s_dt%02d_mom_%+d_%+d_%+d.%s.h5",
		       rp.corr_dir, der_order, part_str[flav], srcstr, smrstr_f[flav], apestr, 
		       thrp_snk.dt, mom_vec[0], mom_vec[1], mom_vec[2], flav_str[flav]);
	      char *group;
	      asprintf(&group, "thrp/%s/%s/dt%02d/%s",
		       srcstr, proj_to_str(thrp_snk.proj),
		       thrp_snk.dt, flav_str[flav]);
	      thrp.corr.lat = corr.lat;
	      thrp.corr.site_size = corr.site_size;
	      int origin[4] = {corr.origin[0], corr.origin[1], corr.origin[2], corr.origin[3]};
	      thrp.corr.origin = origin;
	      if(am_io_proc)
		printf("Preparing to write file %s\n", fname);
	      for(int id=0; id < corr.ncorr; id++) {
		thrp.corr.C = corr.C[id];
		if(thrp.corr.C != NULL)
		  qhg_correlator_shift(thrp.corr, thrp.corr.origin);
	      }
	      for(int i=0; i < 4; i++)
		corr.origin[i] = 0;

	      if( der_order == 2 || der_order == 3) {
		qhg_write_mesons_averaged_thrp_der(fname, corr, group);
	      } else {
		qhg_write_mesons_thrp_der(fname, corr, group);
	      }
	      
	      for(int id=0; id < corr.ncorr; id++)
		if(corr.C[id] != NULL)
		  free(corr.C[id]);

	      if(am_io_proc)
		printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0));
	      free(fname);
	      free(group);
	    }
	  }
	}
	if(am_io_proc)
	  printf("Done sink: sink-source = %2d, in %g sec\n",
		 thrp_snk.dt, qhg_stop_watch(sink_timer));
      }
    }
      
    // Free sequential source and solution
    
    for(int i=0; i<NC*NS; i++) {
      qhg_spinor_field_finalize(seq_src[i]);
      qhg_spinor_field_finalize(seq_sol[i]);
    }
    free(srcstr);
    if(am_io_proc)
      printf("Done source, coords (t, x, y, z) = (%d, %d, %d, %d), in %g sec\n",
	     sco[0], sco[1], sco[2], sco[3], qhg_stop_watch(source_timer));
    
  }
  free(smrstr_l);
  free(smrstr_s);
  
  mg_finalize();

  // Destroy spinor- and gauge-fields

  for(int i=0; i<NS*NC; i++) {
    qhg_spinor_field_finalize(sol_u[i]);
    qhg_spinor_field_finalize(sol_s[i]);
    qhg_spinor_field_finalize(sol_sm_u[i]);
    qhg_spinor_field_finalize(sol_sm_s[i]);
  }

  qhg_gauge_field_finalize(gf);  
  qhg_gauge_field_finalize(gf_ape);  
  qhg_lattice_finalize(lat);
  qhg_comms_finalize(comms);

  return 0;
}
  
