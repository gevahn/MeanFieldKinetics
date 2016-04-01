import numpy as np
import csv
import subprocess as sp

with open('rates.csv','rb') as csvfile:
	rateReader = csv.reader(csvfile, delimiter='\t')
	rateReader.next()
	for rowS in rateReader:
		row = [float(x) for x in rowS]
		sp.call(["./KMC","-n %d"%(int(row[0])),"--ka %f"%(row[1]),"--kd %f"%(row[2]),"--kdiff %f"%(row[3]),"--kr %f"%(row[4]),"--Keq %f"%(row[5]),"--mean %f"%(row[6]),"--sigma %f"%(row[7])]) 







