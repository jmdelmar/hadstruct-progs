#include <qhg.h>
#include <parser_types.h>
#include <mg4qcd.h>
#include <mg4qcd_interface_types.h>

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

mg4qcd_state
mg4qcd_init(struct run_params rp, qhg_gauge_field gf)
{
  for(int i=0; i<ND; i++)
    dims[i] = gf.lat->dims[i];

  for(int i=0; i<ND; i++)
    procs[i] = gf.lat->comms->proc_dims[i];
  
  mg4qcd_state state;
  MG4QCD_Init mg_init;
  MG4QCD_Parameters mg_params;
  MG4QCD_Status mg_status;
  mg_init.comm_cart = gf.lat->comms->comm;
  mg_init.Cart_rank = cart_rank;
  mg_init.Cart_coords = cart_coords;
  mg_init.global_lattice[0] = gf.lat->dims[0];
  mg_init.global_lattice[1] = gf.lat->dims[3];
  mg_init.global_lattice[2] = gf.lat->dims[2];
  mg_init.global_lattice[3] = gf.lat->dims[1];
  mg_init.procs[0] = gf.lat->comms->proc_dims[0];
  mg_init.procs[1] = gf.lat->comms->proc_dims[3];
  mg_init.procs[2] = gf.lat->comms->proc_dims[2];
  mg_init.procs[3] = gf.lat->comms->proc_dims[1];
  mg_init.bc = 2;
  mg_init.block_lattice[0] = rp.mg.block[0];
  mg_init.block_lattice[3] = rp.mg.block[1];
  mg_init.block_lattice[2] = rp.mg.block[2];
  mg_init.block_lattice[1] = rp.mg.block[3];
  mg_init.number_openmp_threads = gf.lat->comms->nthreads;
  mg_init.number_of_levels = rp.mg.n_levels;
  mg_init.kappa = rp.act.kappa;
  mg_init.mu = rp.act.mu;
  mg_init.csw = rp.act.csw;
  MG4QCD_init(&mg_init, &mg_params, &mg_status);  

  for(int i=0; i<rp.mg.n_levels; i++) {
    mg_params.setup_iterations[i] = rp.mg.setup_iterations[i];
    mg_params.mg_basis_vectors[i] = rp.mg.n_basis_vectors[i];
    mg_params.coarse_mu[i] = rp.mg.coarse_mu[i];
  }
  
  mg_params.print = rp.mg.verbosity;
  mg_params.conf_index_fct = conf_index_fct;  
  mg_params.vector_index_fct = vector_index_fct;  
  MG4QCD_update_parameters(&mg_params, &mg_status);
  MG4QCD_set_configuration( (double*) &(gf.field[0]), &mg_status);
  if(gf.lat->comms->proc_id == 0)
    printf("Plaquette according to MG4QCD = %12.10f\n", mg_status.info);
  MG4QCD_setup(&mg_status);
  MG4QCD_update_parameters(&mg_params, &mg_status);
  
  state.mg_init = mg_init;
  state.mg_params = mg_params;
  state.mg_status = mg_status;
  state.current_mu_sign = (rp.act.mu >= 0) - (rp.act.mu < 0) == 1 ? plus : minus;
  return state;
}

void
mg4qcd_invert(qhg_spinor_field x, qhg_spinor_field b, double eps, enum mu_sign s, mg4qcd_state *state)
{
  if(s != state->current_mu_sign) {
    state->mg_params.mu = -state->mg_params.mu;
    for(int i=0; i<state->mg_params.number_of_levels; i++) {
      state->mg_params.coarse_mu[i] = -state->mg_params.coarse_mu[i];
    }
    MG4QCD_update_parameters(&state->mg_params, &state->mg_status);
    state->current_mu_sign = s;
  }

  MG4QCD_solve((double *)x.field, (double *)b.field, eps, &state->mg_status);  
  return;
}

void
mg4qcd_finalize(void)
{
  MG4QCD_finalize();
  return;
}
