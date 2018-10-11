import sys
import numpy as np

r = np.random.randint
L,T = 64,128
rpos = lambda : tuple((r(L),r(L),r(L),r(T)))
nsrc = 32

#
# seed for different replicas
#
seeds = {"0": 321778,
         "1": 908912,
         "2": 743627,}

repl = sys.argv[1]
traj = int(sys.argv[2])
np.random.seed(traj+seeds[repl])
buf = list()
for i in range(nsrc):
    buf.append("sx%02dsy%02dsz%02dst%03d" % (rpos()))
print(" ".join(buf))
