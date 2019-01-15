import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import re

L = 32
T = 64
procs = 4,4,4,2

mu = 0.003
mu_l = 0.003
mu_s = 0.022
csw = 1.74
kappa = 0.1400645

n_stout = 10
omega_stout = 0.1315

repl = sys.argv[1]
traj = sys.argv[2]

conf_dir = "/home/tuf47161/scratch/cA211.30.32/Confs/"
corr_dir = "/home/tuf47161/scratch/cA211.30.32/Corr/%s-%s" % (repl, traj)

tree = ET.Element("hadstruct-input")
ET.SubElement(tree, "dims").text = "%d %d %d %d" % (T, L, L, L)
ET.SubElement(tree, "procs").text = "%d %d %d %d" % tuple(procs)
ET.SubElement(tree, "config").text = conf_dir + ("conf.%s-%04d" % (repl, int(traj)))
ET.SubElement(tree, "corrs-dir").text = corr_dir
s = ET.SubElement(tree, "smearing")
ET.SubElement(s, "n_stout").text = str(n_stout)
ET.SubElement(s, "omega_stout").text = str(omega_stout)
a = ET.SubElement(tree, "action")
ET.SubElement(a, "mu").text = str(mu)
ET.SubElement(a, "kappa").text = str(kappa)
ET.SubElement(a, "csw").text = str(csw)

flat = ET.tostring(tree, "utf-8")
frmt = minidom.parseString(flat)
print(frmt.toprettyxml(indent = "   "))
