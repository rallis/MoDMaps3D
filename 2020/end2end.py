#######################################################
## MoDMaps3D                                         ##
## MultiDimensional Scaling in Python                ##
## Coded by Rallis Karamichalis                      ##
## Github Repo: https://github.com/rallis/MoDMaps3D  ##
#######################################################


import math, sys, os, time
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

import numpy as np

import requests

groups = [
	{
		"ids": "NC_017870,NC_017871,NC_015788,NC_015790,NC_015791,NC_015792,NC_015794,NC_015795,NC_015796,NC_014568,NC_014571,NC_014572,NC_013825,NC_013762,NC_012430,NC_010224,NC_009335,NC_008090,NC_008085,NC_008088,NC_008089,NC_008091,NC_008076,NC_008077,NC_008078,NC_008079,NC_008080,NC_008081,NC_008082,NC_008083,NC_008084,NC_007446,NC_006887,NC_006888,NC_006889,NC_006890,NC_006407,NC_006325,NC_006326,NC_006327,NC_006328,NC_006329,NC_006330,NC_006331,NC_006332,NC_006333,NC_006334,NC_006335,NC_006336,NC_006337,NC_006338,NC_006339,NC_006340,NC_006341,NC_006342,NC_006343,NC_006344,NC_006345,NC_006346,NC_005797,NC_004926,NC_004021,NC_002756",
		"name": "Caudata",
		"color": "green"
	},
	{
		"ids": "NC_007911,NC_006404,NC_006301,NC_006302,NC_006303,NC_006304,NC_006305,NC_002471",
		"name": "Gymnophiona",
		"color": "cyan"
	},
	{
		"ids": "NC_016119,NC_016059,NC_015615,NC_015617,NC_015618,NC_015620,NC_015305,NC_014685,NC_014691,NC_014581,NC_014584,NC_013270,NC_012837,NC_012647,NC_011049,NC_010232,NC_010233,NC_009886,NC_009422,NC_009423,NC_009258,NC_009264,NC_008410,NC_008144,NC_007888,NC_007440,NC_007178,NC_006839,NC_006688,NC_006689,NC_006690,NC_006402,NC_006403,NC_006405,NC_006406,NC_006408,NC_008975,NC_005794,NC_005055,NC_002805,NC_001573",
		"name": "Anura",
		"color": "blue"
	},
]

groups = [
	{
		"ids": "AY521630,AM000053,AM000054,AY521629,AY521631,AM000055,DQ396400,AF004394,AF042101,AF256204,AF086817,AF042103,AB428555,AB287363,EU786678,U46016,AY713415,AY713416,AY255826,EF514713,DQ207941,DQ369994,EU786673,AY773338,AY773341,EF633445,AY322189,AF484516,AJ488927,AY371157,AY795907,U88826,AF061642,AY772535,AY586549,AF423760,AY371121,AB231893,EU786670",
		"name": "test",
		"color": "testcolor"
	}
]


def request_with_retry(retries=5, backoff_factor=0.5, status_forcelist=(500, 502, 504), session=None):
	session = session or requests.Session()
	retry = Retry(
		total=retries,
		read=retries,
		connect=retries,
		backoff_factor=backoff_factor,
		status_forcelist=status_forcelist,
	)
	adapter = HTTPAdapter(max_retries=retry)
	session.mount('http://', adapter)
	session.mount('https://', adapter)
	return session

def parse_seq_from_fast(str):
	lines = str.split("\n")
	print("lines=", len(lines))
	output = ""
	for ln in lines:
		if not ln.startswith(">"):
			output += ln
	# print(output)
	return output

def replaceMultiple(mainString, toBeReplaces, newString):
    for elem in toBeReplaces :
        if elem in mainString :
            mainString = mainString.replace(elem, newString)
    
    return  mainString

def build_fcgr(seq, k):
	num_kmers = 0
	kmers = {}
	mapACGT = {"A": 0, "C": 1, "G": 2, "T": 3}

	newseq = replaceMultiple(seq, ["B", "D", "E", "F", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "U", "V", "W", "X", "Y", "Z"], "")
	print("len seq=", len(seq), "new seq=", len(newseq))

	fcgr_dim = int(math.pow(4,k))
	kmer_inds = np.zeros(fcgr_dim, dtype=int)
	# print(kmer_inds)

	for i in range(len(newseq) - k + 1):
		curmer = newseq[i:i+k]
		curmerInd = 0
		for j in range(len(curmer)):
			curmerInd += int(math.pow(4, k-1-j)*mapACGT[curmer[j]]) 
		
		if kmer_inds[curmerInd]==0:
			num_kmers += 1
			kmer_inds[curmerInd] = 1    # for AID

	return {"size":num_kmers, "flatFCGR":kmer_inds}

###############################################
root_dir = os.getcwd()
print("RootDir={}".format(root_dir))

fcgr_resolution = 9
all_fcgrs = []
for gr in groups:
	print(gr["name"], gr["color"])
	accession_numbers = gr["ids"].split(",")
	print("Contains {} accession_numbers".format(len(accession_numbers)))

	for ind, acc_num in enumerate(accession_numbers):
		print("{}/{} {}".format(ind, len(accession_numbers), acc_num))
		r = request_with_retry().get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={}".format(acc_num))
		print(acc_num, r.status_code)
		if r.status_code != 200:
			raise Exception("error")
					
		with open(os.path.join(root_dir, "data", "{}.fasta".format(acc_num)), "w") as fout:
			fout.write(r.text)
		
		seq = parse_seq_from_fast(r.text)
		all_fcgrs.append(build_fcgr(seq, fcgr_resolution))


print(len(all_fcgrs))
print(all_fcgrs[0])

dist_matrix = np.zeros((len(all_fcgrs),len(all_fcgrs)), dtype=np.float64)
for i in range(len(all_fcgrs)):
	for j in range(i+1, len(all_fcgrs)):
		numerator = all_fcgrs[i]['size'] + all_fcgrs[j]['size']
		denominator = 0
		for union_ind in range(int(math.pow(4,fcgr_resolution))):
			if (all_fcgrs[i]['flatFCGR'][union_ind] == 1) or  (all_fcgrs[j]['flatFCGR'][union_ind] == 1):
				denominator += 1

		val = 2- numerator*1.0/denominator
		print(i,j, numerator, denominator, val)
		dist_matrix[i][j] = val 
		dist_matrix[j][i] = val 
		


print(len(dist_matrix))
print(dist_matrix)




print("Computing eigenvalues..\n")
eigValues, eigVectors = np.linalg.eig(dist_matrix)
idx = eigValues.argsort()[::-1][0:5]  
selEigValues = eigValues[idx]
selEigVectors = eigVectors[:,idx]
print(eigValues)
print(idx) 
print(selEigValues)
print(selEigValues.shape)
print(selEigVectors.shape)

if False in (selEigValues > 0):
    print("First 5 largest eigenvalues are not all positive. Exiting..")
    sys.exit(-1)


selEigVectors = np.array(selEigVectors)

diagValues = []
for i in range(len(selEigValues)):
    diagValues.append(math.sqrt(eigValues[i]))
print(diagValues)

diag = np.diag(diagValues)
points = np.dot(selEigVectors,diag)
print("pointsSize=", points.shape)

minmaxScaling = []
for i in range(5):
	minmaxScaling.append([ min(points[:,i]), max(points[:,i]) ])
print(minmaxScaling)

scaledPoints = []
for i in range(len(accessions)):
	scaledPoints.append([0, 0, 0, 0, 0])
	for j in range(5):
		scaledPoints[i][j] = 2.0 *(points[i][j] - minmaxScaling[j][0]) / ( minmaxScaling[j][1] - minmaxScaling[j][0]) - 1

print("Writing output [map_output.txt]\n")
for i in range(len(accessions)):
	map_output += str(scaledPoints[i][0]) + "\n"
	map_output += str(scaledPoints[i][1]) + "\n"
	map_output += str(scaledPoints[i][2]) + "\n"
	map_output += str(scaledPoints[i][3]) + "\n"
	map_output += str(scaledPoints[i][4]) + "\n"
	map_output += str(i) + "\n"
	map_output += allSequences[i][0] + "\n"
	map_output += allSequences[i][3].split("|")[-1] + "\n"
	map_output += allSequences[i][1] + "\n"
	map_output += allSequences[i][2] + "\n"


fout = open("map_output.txt", "w")
fout.write(map_output)
fout.close()





# map_output = ""


# print ("Reading [mapheader.txt]\n"
# header = open ("mapheader.txt", "r")
# for line in header:
# 	map_output += line
# header.close()


# print ("Reading [accIDs.txt]\n"
# accIDs = open ("accIDs.txt", "r")
# accessions = accIDs.readlines()
# accIDs.close()	
# tmp = accessions[0].split('"')
# accessions = [tmp[it] for it in range(1,len(tmp),2)]


# print ("Reading [allSequences.txt]\n"
# allSeq = open ("allSequences.txt", "r")
# allSequences = allSeq.readlines()
# allSeq.close()	
# tmp = allSequences[0].split('"')
# allSequences = [[tmp[it], tmp[it+2], tmp[it+4], tmp[it+6]] for it in range(1,len(tmp),8)]


# print ("Reading [dMatrix.txt]\n"
# dMat = open ("dMatrix.txt", "r")
# distMat = dMat.readlines()
# dMat.close()
# tmp = distMat[0].split('{')
# dMat=[]
# for row in range(2,len(tmp)):
# 	dMat.append(tmp[row].split("}")[0].split(","))
# dMat = np.array(dMat, dtype="float64")


# print ("Computing eigenvalues..\n"
# eigValues, eigVectors = np.linalg.eig(dMat)
# idx = eigValues.argsort()[::-1][0:5]  
# selEigValues = eigValues[idx]
# selEigVectors = eigVectors[:,idx]
# #print (idx 
# #print (selEigValues.shape
# #print (selEigVectors.shape

# if False in (selEigValues > 0):
#     print ("First 5 largest eigenvalues are not all positive. Exiting.."
#     sys.exit(-1)


# selEigVectors = np.array(selEigVectors)

# diagValues = []
# for i in range(len(selEigValues)):
#     diagValues.append(math.sqrt(eigValues[i]))
# #print (diagValues
    
# diag = np.diag(diagValues)
# points = np.dot(selEigVectors,diag)
# #print ("pointsSize=", points.shape

# minmaxScaling = []
# for i in range(5):
# 	minmaxScaling.append([ min(points[:,i]), max(points[:,i]) ])
# #print (minmaxScaling

# scaledPoints = []
# for i in range(len(accessions)):
# 	scaledPoints.append([0, 0, 0, 0, 0])
# 	for j in range(5):
# 		scaledPoints[i][j] = 2.0 *(points[i][j] - minmaxScaling[j][0]) / ( minmaxScaling[j][1] - minmaxScaling[j][0]) - 1

# print ("Writing output [map_output.txt]\n"
# for i in range(len(accessions)):
# 	map_output += str(scaledPoints[i][0]) + "\n"
# 	map_output += str(scaledPoints[i][1]) + "\n"
# 	map_output += str(scaledPoints[i][2]) + "\n"
# 	map_output += str(scaledPoints[i][3]) + "\n"
# 	map_output += str(scaledPoints[i][4]) + "\n"
# 	map_output += str(i) + "\n"
# 	map_output += allSequences[i][0] + "\n"
# 	map_output += allSequences[i][3].split("|")[-1] + "\n"
# 	map_output += allSequences[i][1] + "\n"
# 	map_output += allSequences[i][2] + "\n"


# fout = open("map_output.txt", "w")
# fout.write(map_output)
# fout.close()

