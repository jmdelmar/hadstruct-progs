#include <qhg.h>
#include <DDalphaAMG.h>
#include <mg_interface_types.h>

static int dims[ND];
static int procs[ND];

static int
conf_index_fct(int t, int z, int y, int x, int mu) {
  int d[] = {
    dims[0]/procs[0],
    dims[1]/procs[1],
    dims[2]/procs[2],
    dims[3]/procs[3],
  };
  int co[] = {t,x,y,z};
  int pos = IDX(co, d);
  static int size_per_pos = ND*NC*NC*2;
  int dir = (mu==1)?3:((mu==3)?1:mu);
  int size_per_dir = NC*NC*2;     
  return pos*size_per_pos + dir*size_per_dir;
}

static int
vector_index_fct(int t, int z, int y, int x) {
  int d[] = {
    dims[0]/procs[0],
    dims[1]/procs[1],
    dims[2]/procs[2],
    dims[3]/procs[3],
  };
  int co[] = {t,x,y,z};
  int pos = IDX(co, d);
  int size_per_pos = NS*NC*2;
  return pos*size_per_pos;
}

static int
cart_rank(MPI_Comm comm, const int *c, int *rank) {
  int coords[] = {c[0], c[3], c[2], c[1]};
  int ierr = MPI_Cart_rank(comm, coords, rank);
  return ierr;
}

static int
cart_coords(MPI_Comm comm, int rank, int maxrank, int coords[ND]) {
  int c[ND];
  int ierr = MPI_Cart_coords(comm, rank, maxrank, c);
  coords[0] = c[0];
  coords[3] = c[1];
  coords[2] = c[2];
  coords[1] = c[3];
  return ierr;
}

mg_state
mg_init(struct multigrid_params mg, struct action_params act, qhg_gauge_field gf)
{
  for(int i=0; i<ND; i++)
    dims[i] = gf.lat->dims[i];

  for(int i=0; i<ND; i++)
    procs[i] = gf.lat->comms->proc_dims[i];
  
  mg_state state;
  DDalphaAMG_init init;
  DDalphaAMG_parameters params;
  DDalphaAMG_status status;
  init.comm_cart = gf.lat->comms->comm;
  init.Cart_rank = cart_rank;
  init.Cart_coords = cart_coords;
  init.global_lattice[0] = gf.lat->dims[0];
  init.global_lattice[1] = gf.lat->dims[3];
  init.global_lattice[2] = gf.lat->dims[2];
  init.global_lattice[3] = gf.lat->dims[1];
  init.procs[0] = gf.lat->comms->proc_dims[0];
  init.procs[1] = gf.lat->comms->proc_dims[3];
  init.procs[2] = gf.lat->comms->proc_dims[2];
  init.procs[3] = gf.lat->comms->proc_dims[1];
  init.bc = 1;
  init.block_lattice[0] = mg.block[0];
  init.block_lattice[3] = mg.block[1];
  init.block_lattice[2] = mg.block[2];
  init.block_lattice[1] = mg.block[3];
  init.number_openmp_threads = gf.lat->comms->nthreads;
  init.number_of_levels = mg.n_levels;
  init.kappa = act.kappa;
  init.mu = act.mu;
  init.csw = act.csw;
  init.init_file = NULL;
  init.rnd_seeds = NULL;
  DDalphaAMG_initialize(&init, &params, &status);  


  for(int i=0; i<mg.n_levels; i++) {
    params.setup_iterations[i] = mg.setup_iterations[i];
    params.mg_basis_vectors[i] = mg.n_basis_vectors[i];
    params.mu_factor[i] = mg.coarse_mu[i]/act.mu;
  }

  params.mu_odd_shift = 0;
  params.mu_even_shift = 0;
  params.mixed_precision = 1;
  params.print = mg.verbosity;
  params.conf_index_fct = conf_index_fct;  
  params.vector_index_fct = vector_index_fct;  
  DDalphaAMG_update_parameters(&params, &status);

  DDalphaAMG_set_configuration( (double*) &(gf.field[0]), &status);

  
  DDalphaAMG_setup(&status);
  DDalphaAMG_update_parameters(&params, &status);
  
  state.init = init;
  state.params = params;
  state.status = status;
  state.current_mu_sign = (act.mu >= 0) - (act.mu < 0) == 1 ? plus : minus;
  return state;
}

void
mg_invert(qhg_spinor_field x, qhg_spinor_field b, double eps, enum mu_sign s, mg_state *state)
{
  if(s != state->current_mu_sign) {
    state->params.mu = -state->params.mu;
    /* for(int i=0; i<state->init.number_of_levels; i++) { */
    /*   state->params.mu_factor[i] = -state->params.mu_factor[i]; */
    /* } */
    DDalphaAMG_update_parameters(&state->params, &state->status);
    state->current_mu_sign = s;
  }

  DDalphaAMG_solve((double *)x.field, (double *)b.field, eps, &state->status);
  int n_coarse = state->status.coarse_iter_count;
  double t_coarse = state->status.coarse_time / n_coarse;
  int n_lev = state->init.number_of_levels;
  int coarse_dims[ND];
  for(int d=0; d<ND; d++)
    coarse_dims[d] = state->init.global_lattice[d];

  for(int i=0; i<n_lev; i++)
    for(int d=0; d<ND; d++)
      coarse_dims[d] /= state->params.block_lattice[i][d];

  int coarse_vol = 1;
  for(int d=0; d<ND; d++)
    coarse_vol *= coarse_dims[d];

  int nv = state->params.mg_basis_vectors[n_lev-2];
  double perf = (double)(coarse_vol*(4.0*nv-1.0)*nv*108.0)/t_coarse;
  if(x.lat->comms->proc_id == 0)
    printf(" Coarse grid peformance: %12.2f Gflop/s\n\n\n", perf/1e9);
  return;
}

void
mg_finalize(void)
{
  DDalphaAMG_finalize();
  return;
}
