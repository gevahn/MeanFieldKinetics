import subprocess as sp


for sigma in range(100):
	sp.call(["./KMC","-n 40", "--ka 0.1", "--kd 0.1", "--kdiff 0", "--kr 0", "--Keq 50", "--sigma %d"%(sigma),"--mean -6"])
