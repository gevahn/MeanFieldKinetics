import subprocess as sp
import sys

for s in range(50):
	sigma = s * 0.2
	sp.call(["./KMC","-n 40", "--ka " + sys.argv[1], "--kd " + sys.argv[2], "--kdiff 1", "--kr " + sys.argv[3], "--Keq 50", "--sigma %d"%(sigma),"--mean -6"])
