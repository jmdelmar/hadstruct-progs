#include <stdio.h>
#include <stdlib.h>
#include <mxml.h>
#include <parser_types.h>

/*
 * Trims leading and trailing white-space
 */
static char *
trim(const char *in_str)
{
  char *str = strdup(in_str);  
  char *end;

  /* Trim leading space */
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  /* Trim trailing space */
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  /* Write new null terminator */
  *(end+1) = 0;

  return str;
}

/*
 * Function to use to exit if tag not found
 */
static void
exit_tag_not_found(char tag[], char fname[])
{
  fprintf(stderr, "</%s>: not found in %s\n", tag, fname);
  exit(1);
}

/*
 * Function to use to exit if value conversion fails
 */
static void
exit_convertions_error(char tag[], char fname[], const char val[])
{
  fprintf(stderr, "</%s>: got %s which seems to have incorrect type in %s\n", tag, val, fname);
  exit(1);
}

/*
 * Main external function. Parses input file. Returns parameters in
 * single struct
 */
struct run_params
parse_input(char fname[])
{
  FILE *fp = fopen(fname, "r");
  if(fp == NULL) {
    fprintf(stderr, "%s: error opening file\n", fname);
    exit(1);
  }
  
  /* Find the top-level "hadstruct-input" element */
  mxml_node_t *tree = mxmlLoadFile(NULL, fp, MXML_OPAQUE_CALLBACK);
  mxml_node_t *node = mxmlFindElement(tree, tree, "hadstruct-input", NULL, NULL, MXML_DESCEND);
  const char *config = mxmlGetOpaque(mxmlFindPath(node, "config"));
  const char *corr_dir = mxmlGetOpaque(mxmlFindPath(node, "corrs-dir"));
  const char *prop_dir = mxmlGetOpaque(mxmlFindPath(node, "props-dir"));
  const char *d = mxmlGetOpaque(mxmlFindPath(node, "dims"));  
  const char *p = mxmlGetOpaque(mxmlFindPath(node, "procs"));  
  if(config == NULL) {
    exit_tag_not_found("config", fname);
  }
  if(corr_dir == NULL) {
    exit_tag_not_found("corrs-dir", fname);
  }
  if(prop_dir == NULL) {
    exit_tag_not_found("props-dir", fname);
  }
  if(d == NULL) {
    exit_tag_not_found("dims", fname);
  }
  if(p == NULL) {
    exit_tag_not_found("procs", fname);
  }

  struct action_params a;
  {
    /* 
       Get the parameters of the action 
    */
    mxml_node_t *anode = mxmlFindElement(node, node, "action", NULL, NULL, MXML_DESCEND);
    /* First as strings */
    const char *mu = mxmlGetOpaque(mxmlFindPath(anode, "mu"));
    const char *mu_l = mxmlGetOpaque(mxmlFindPath(anode, "mu_l"));
    const char *mu_s = mxmlGetOpaque(mxmlFindPath(anode, "mu_s"));
    const char *kappa = mxmlGetOpaque(mxmlFindPath(anode, "kappa"));
    const char *csw = mxmlGetOpaque(mxmlFindPath(anode, "csw"));
    if(mu == NULL) {
      mu = mu_l;
      if(mu_l == NULL || mu_s == NULL) {
	exit_tag_not_found("action/mu", fname);
      }
    }
    if(kappa == NULL) {
      exit_tag_not_found("action/kappa", fname);
    }
    if(csw == NULL) {
      exit_tag_not_found("action/csw", fname);
    }

    /* Convert to float. Check for errors */
    char *e;
    a.mu = strtod(mu, &e);
    if(e == mu) {
      exit_convertions_error("action/mu", fname, mu);
    } 
    a.mu_l = strtod(mu_l, &e);
    if(e == mu_l) {
      exit_convertions_error("action/mu_l", fname, mu_l);
    } 
    a.mu_s = strtod(mu_s, &e);
    if(e == mu_s) {
      exit_convertions_error("action/mu_s", fname, mu_s);
    } 
    a.kappa = strtod(kappa, &e);
    if(e == kappa) {
      exit_convertions_error("action/kappa", fname, kappa);
    }
    a.csw = strtod(csw, &e);
    if(e == csw) {
      exit_convertions_error("action/csw", fname, csw);
    }
  }
 
  struct multigrid_params mg;
  {
    /* 
       Get the multigrid parameters
    */
    mxml_node_t *mnode = mxmlFindElement(node, node, "multi-grid", NULL, NULL, MXML_DESCEND);
    /* First as strings */
    const char *verb = mxmlGetOpaque(mxmlFindPath(mnode, "verbosity"));
    const char *nl = mxmlGetOpaque(mxmlFindPath(mnode, "n_levels"));
    const char *bl = mxmlGetOpaque(mxmlFindPath(mnode, "block"));
    if(verb == NULL) {
      exit_tag_not_found("multi-grid/verbosity", fname);
    }
    if(nl == NULL) {
      exit_tag_not_found("multi-grid/n_levels", fname);
    }
    if(bl == NULL) {
      exit_tag_not_found("multi-grid/block", fname);
    }


    sscanf(bl, "%d %d %d %d",
	   &mg.block[0],
	   &mg.block[1],
	   &mg.block[2],
	   &mg.block[3]);
    
    /* Convert. Check for errors */
    char *e;
    mg.n_levels = strtoul(nl, &e, 10);
    if(e == nl) {
      exit_convertions_error("multi-grid/n_levels", fname, nl);
    } 

    mg.verbosity = strtoul(verb, &e, 10);
    if(e == nl) {
      exit_convertions_error("multi-grid/verbosity", fname, verb);
    }

    for(int il=0; il<mg.n_levels; il++) {
      char *ltag;
      asprintf(&ltag, "level_%d", il+1);
      mxml_node_t *lnode = mxmlFindElement(mnode, mnode, ltag, NULL, NULL, MXML_DESCEND);
      const char *n_iter = mxmlGetOpaque(mxmlFindPath(lnode, "setup_iters"));
      const char *n_base = mxmlGetOpaque(mxmlFindPath(lnode, "n_basis_vectors"));
      const char *mu = mxmlGetOpaque(mxmlFindPath(lnode, "mu"));
      
      if(n_base == NULL) {
	char *e;
	asprintf(&e, "multi-grid/%s/n_basis_vectors", ltag);	
	exit_tag_not_found(e, fname);
      }
      if(n_iter == NULL) {
	char *e;
	asprintf(&e, "multi-grid/%s/setup_iters", ltag);	
	exit_tag_not_found(e, fname);
      }
      if(mu == NULL) {
	char *e;
	asprintf(&e, "multi-grid/%s/mu", ltag);	
	exit_tag_not_found(e, fname);
      }
      /* Convert. Check for errors */
      char *e;
      mg.setup_iterations[il] = strtoul(n_iter, &e, 10);
      if(e == n_iter) {
	char *e;
	asprintf(&e, "multi-grid/%s/setup_iters", ltag);	
	exit_convertions_error(e, fname, n_iter);
      } 
      mg.n_basis_vectors[il] = strtoul(n_base, &e, 10);
      if(e == n_base) {
	char *e;
	asprintf(&e, "multi-grid/%s/n_basis_vectors", ltag);	
	exit_convertions_error(e, fname, n_base);
      } 
      mg.coarse_mu[il] = strtod(mu, &e);
      if(e == mu) {
	char *e;
	asprintf(&e, "multi-grid/%s/mu", ltag);	
	exit_convertions_error(e, fname, mu);
      } 
      
    }
  }
 
  /* Enumerate "sp" (source position) elements */
  mxml_index_t *ind = mxmlIndexNew(node, "sp", NULL);
  int nsp = mxmlIndexGetCount(ind);
  struct source_position *r = malloc(sizeof(struct source_position)*(nsp));

  {
    /* Iterate over "sp" elements */
    int isp = 0;
    for(mxml_node_t *n=mxmlIndexEnum(ind); n!=NULL; n=mxmlIndexEnum(ind)) {
      const char *c = mxmlGetOpaque(mxmlFindPath(n, "coords"));
      if(c == NULL) {
	exit_tag_not_found("coords", fname);
      }
      sscanf(c, "%d %d %d %d",
	     &r[isp].coords[0],
	     &r[isp].coords[1],
	     &r[isp].coords[2],
	     &r[isp].coords[3]);
      /* Enumerate "vec" elements (if any) */
      mxml_index_t *mom_index = mxmlIndexNew(n, "mom", NULL);
      if(mom_index == NULL) {
	exit_tag_not_found("mom", fname);
	}
      int nmoms = mxmlIndexGetCount(mom_index);
      r[isp].nmoms = nmoms;
      int imom = 0;
      /* Iterate over "mom" elements */
      for(mxml_node_t *m=mxmlIndexEnum(mom_index); m!=NULL; m=mxmlIndexEnum(mom_index)) {
	const char *vec = mxmlGetOpaque(mxmlFindPath(m, "vec"));
	sscanf(vec, "%d %d %d",
	       &r[isp].mom_vecs[imom][0],
	       &r[isp].mom_vecs[imom][1],
	       &r[isp].mom_vecs[imom][2]);
	imom++;
      }
      /* Enumerate "sink" elements (if any) */
      mxml_index_t *it = mxmlIndexNew(n, "sink", NULL);
      int nsnk = mxmlIndexGetCount(it);
      r[isp].nsinks = nsnk;
      r[isp].sinks = malloc(sizeof(qhg_thrp_nn_sink_params)*(nsnk));
      int isnk = 0;
      /* Iterate over "sink" elements */
      for(mxml_node_t *s=mxmlIndexEnum(it); s!=NULL; s=mxmlIndexEnum(it)) {
	r[isp].sinks[isnk].proj = str_to_proj((char *)trim(mxmlGetOpaque(mxmlFindPath(s, "proj"))));
	r[isp].sinks[isnk].dt = atoi(mxmlGetOpaque(mxmlFindPath(s, "dt")));
	r[isp].sinks[isnk].dt = atoi(mxmlGetOpaque(mxmlFindPath(s, "dt")));
	isnk++;
      }
      isp++;
    }
  }

  struct smearing_params s;
  {
    /* 
       Get the smearing parameters
    */
    mxml_node_t *snode = mxmlFindElement(node, node, "smearing", NULL, NULL, MXML_DESCEND);
    /* First as strings */
    const char *alpha_ape = mxmlGetOpaque(mxmlFindPath(snode, "alpha_ape"));
    const char *alpha_gauss = mxmlGetOpaque(mxmlFindPath(snode, "alpha_gauss"));
    const char *alpha_gauss_l = mxmlGetOpaque(mxmlFindPath(snode, "alpha_gauss_l"));
    const char *alpha_gauss_s = mxmlGetOpaque(mxmlFindPath(snode, "alpha_gauss_s"));
    const char *n_ape = mxmlGetOpaque(mxmlFindPath(snode, "n_ape"));
    const char *n_gauss = mxmlGetOpaque(mxmlFindPath(snode, "n_gauss"));
    const char *n_gauss_l = mxmlGetOpaque(mxmlFindPath(snode, "n_gauss_l"));
    const char *n_gauss_s = mxmlGetOpaque(mxmlFindPath(snode, "n_gauss_s"));

    if(alpha_ape == NULL) {
      exit_tag_not_found("smearing/alpha_ape", fname);
    }
    if(alpha_gauss == NULL) {
      alpha_gauss = alpha_gauss_l;
      if( alpha_gauss_l == NULL || alpha_gauss_s == NULL ) {
	exit_tag_not_found("smearing/alpha_gauss", fname);
      }
    }
    if(n_ape == NULL) {
      exit_tag_not_found("smearing/n_ape", fname);
    }
    if(n_gauss == NULL) {
      n_gauss = n_gauss_l;
      if( n_gauss_l == NULL || n_gauss_s == NULL ) {
	exit_tag_not_found("smearing/n_gauss", fname);
      }
    }

    /* Convert to float. Check for errors */
    char *e;
    s.alpha_ape = strtod(alpha_ape, &e);
    if(e == alpha_ape) {
      exit_convertions_error("smearing/alpha_ape", fname, alpha_ape);
    } 

    s.alpha_gauss = strtod(alpha_gauss, &e);
    if(e == alpha_gauss) {
      exit_convertions_error("smearing/alpha_gauss", fname, alpha_gauss);
    }

    s.alpha_gauss_l = strtod(alpha_gauss_l, &e);
    if(e == alpha_gauss_l) {
      exit_convertions_error("smearing/alpha_gauss_l", fname, alpha_gauss_l);
    }
    
    s.alpha_gauss_s = strtod(alpha_gauss_s, &e);
    if(e == alpha_gauss_s) {
      exit_convertions_error("smearing/alpha_gauss_s", fname, alpha_gauss_s);
    }

    s.n_ape = strtoul(n_ape, &e, 10);
    if(e == n_ape) {
      exit_convertions_error("smearing/n_ape", fname, n_ape);
    } 

    s.n_gauss = strtoul(n_gauss, &e, 10);
    if(e == n_gauss) {
      exit_convertions_error("smearing/n_gauss", fname, n_gauss);
    } 

    s.n_gauss_l = strtoul(n_gauss_l, &e, 10);
    if(e == n_gauss_l) {
      exit_convertions_error("smearing/n_gauss_l", fname, n_gauss_l);
    } 

    s.n_gauss_s = strtoul(n_gauss_s, &e, 10);
    if(e == n_gauss_s) {
      exit_convertions_error("smearing/n_gauss_s", fname, n_gauss_s);
      }
  }

  struct run_params rp;
  sscanf(d, "%d %d %d %d",
	 &rp.dims[0],
	 &rp.dims[1],
	 &rp.dims[2],
	 &rp.dims[3]);
  sscanf(p, "%d %d %d %d",
	 &rp.procs[0],
	 &rp.procs[1],
	 &rp.procs[2],
	 &rp.procs[3]);
  strcpy(rp.prop_dir, trim(prop_dir));
  strcpy(rp.corr_dir, trim(corr_dir));
  strcpy(rp.config, trim(config));
  rp.act = a;
  rp.mg = mg;
  rp.spos = r;
  rp.nsp = nsp;
  rp.smearing = s;
  
  fclose(fp);
  mxmlDelete(tree);
  return rp;
}
