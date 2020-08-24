import sys
import numpy as np

r = np.random.randint
L,T = 32, 64
rpos = lambda : tuple((r(L),r(L),r(L),r(T)))
nsrc = 32

#
# seed for different replicas
#
seeds = {"0": 290159766,
         "1": 127301172,
         "2": 651320871,}

with open( sys.argv[2], "r" ) as conf_file:

    conf_list = conf_file.readlines()

with open( sys.argv[3], 'r' ) as old_src_file:

    old_src_list = old_src_file.readlines()

with open( sys.argv[1], 'w' ) as src_file:

    for conf, old_src in zip( conf_list, old_src_list ):

        [ repl, traj ] = conf.strip().split( '-' )

        np.random.seed(int(traj)+seeds[repl])
        buf = traj  + ":  "

        i = 0
        while i < nsrc:

            src = "sx%02dsy%02dsz%02dst%03d" % (rpos())

            if src not in old_src:

                buf = buf + " " + src
                i+=1

            else:

                print( "DUPLICATE", conf )

        buf = buf + "\n"
                
        src_file.write( buf )
