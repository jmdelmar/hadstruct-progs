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

#define NF 2
char part_str[NF][7] = {"kaon","kaon"};
char flav_str[NF][7] = {"up","strange"};

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
  // Up propagator smeared with n_gauss_s to be
  // consistent with up part kaon 
  // three-point functions
  qhg_spinor_field sol_u_sm_s[NS*NC];
  // Strange propagator smeared with n_gauss_u to be
  // consistent with strange part of kaon 
  // three-point functions
  qhg_spinor_field sol_s_sm_u[NS*NC];
  for(int sp=0; sp<NS; sp++)
    for(int co=0; co<NC; co++) {
      sol_u[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_s[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_sm_u[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_sm_s[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_u_sm_s[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
      sol_s_sm_u[CS(sp,co)] = qhg_spinor_field_init(lat, bc);
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
    int n_gauss_f[] = {n_gauss_l, n_gauss_s};
    double alpha_gauss_f[] = {alpha_gauss_l, alpha_gauss_s};

    if(am_io_proc)
      printf("Smearing the source\n");  
    for(int sp=0; sp<NS; sp++)
      for(int co=0; co<NC; co++) {
	if(am_io_proc)
	  printf("\tcol=%d, spin=%d\n", co, sp);  
	qhg_spinor_field aux = qhg_spinor_field_init(lat, bc);
	qhg_point_spinor_field(aux, sco, sp, co);
	for(int flav=0; flav<NF; flav++) {
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

    for(int flav=0; flav<NF; flav++) {
      mg_state.params.mu = flav == 0 ? rp.act.mu_l : rp.act.mu_s;
      DDalphaAMG_update_parameters(&mg_state.params, &mg_state.status);
      mg_state.current_mu_sign = plus;

      for(int i=0; i<NS*NC; i++) {
	mg_invert(sol_f[flav][i], src_f[flav][i], 1e-9, minus, &mg_state);
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
	qhg_gauss_smear(sol_u_sm_s[CS(sp,co)], sol_u[CS(sp,co)], gf_ape, alpha_gauss_s, n_gauss_s);
	qhg_gauss_smear(sol_s_sm_u[CS(sp,co)], sol_s[CS(sp,co)], gf_ape, alpha_gauss_l, n_gauss_l);
      }
    if(am_io_proc)
      printf("Done smearing in %g sec\n", qhg_stop_watch(t0));  
    
    // Smeared meson correlators and fourier transform

    t0 = qhg_stop_watch(0);
    qhg_correlator corr_kaon_u = qhg_mesons_pseudoscalar(sol_u_sm_s, sol_sm_s, sco);
    qhg_correlator_shift(corr_kaon_u, corr_kaon_u.origin);
    qhg_correlator corr_kaon_s = qhg_mesons_pseudoscalar(sol_sm_u, sol_s_sm_u, sco);
    qhg_correlator_shift(corr_kaon_s, corr_kaon_s.origin);
    if(am_io_proc)
      printf("Done meson correlator in %g sec\n", qhg_stop_watch(t0)); 

    // Write the correlators

    if(true) {
      t0 = qhg_stop_watch(0);      
      char *fname;
      asprintf(&fname, "%s/twop_kaon_%s_%s_%s_%s.h5", rp.corr_dir, srcstr, smrstr_s, smrstr_l, apestr);
      char *group;
      asprintf(&group, "twop_kaon/%s/", srcstr);
      qhg_write_single_meson(fname, corr_kaon_u, group);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 
      t0 = qhg_stop_watch(0);      
      asprintf(&fname, "%s/twop_kaon_%s_%s_%s_%s.h5", rp.corr_dir, srcstr, smrstr_l, smrstr_s, apestr);
      asprintf(&group, "twop_kaon/%s/", srcstr);
      qhg_write_single_meson(fname, corr_kaon_s, group);
      if(am_io_proc)
	printf("Wrote %s in %g sec\n", fname, qhg_stop_watch(t0)); 
      free(fname);
      free(group);
    }
    qhg_correlator_finalize(corr_kaon_u);
    qhg_correlator_finalize(corr_kaon_s);
        
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
    qhg_spinor_field_finalize(sol_u_sm_s[i]);
    qhg_spinor_field_finalize(sol_s_sm_u[i]);
  }

  qhg_gauge_field_finalize(gf);  
  qhg_gauge_field_finalize(gf_ape);  
  qhg_lattice_finalize(lat);
  qhg_comms_finalize(comms);

  return 0;
}
  
