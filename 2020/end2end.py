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

# groups = [
# 	{
# 		"ids": "NC_017870,NC_017871,NC_015788,NC_015790,NC_015791,NC_015792,NC_015794,NC_015795,NC_015796,NC_014568,NC_014571,NC_014572,NC_013825,NC_013762,NC_012430,NC_010224,NC_009335,NC_008090,NC_008085,NC_008088,NC_008089,NC_008091,NC_008076,NC_008077,NC_008078,NC_008079,NC_008080,NC_008081,NC_008082,NC_008083,NC_008084,NC_007446,NC_006887,NC_006888,NC_006889,NC_006890,NC_006407,NC_006325,NC_006326,NC_006327,NC_006328,NC_006329,NC_006330,NC_006331,NC_006332,NC_006333,NC_006334,NC_006335,NC_006336,NC_006337,NC_006338,NC_006339,NC_006340,NC_006341,NC_006342,NC_006343,NC_006344,NC_006345,NC_006346,NC_005797,NC_004926,NC_004021,NC_002756",
# 		"name": "Caudata",
# 		"color": "green"
# 	},
# 	{
# 		"ids": "NC_007911,NC_006404,NC_006301,NC_006302,NC_006303,NC_006304,NC_006305,NC_002471",
# 		"name": "Gymnophiona",
# 		"color": "cyan"
# 	},
# 	{
# 		"ids": "NC_016119,NC_016059,NC_015615,NC_015617,NC_015618,NC_015620,NC_015305,NC_014685,NC_014691,NC_014581,NC_014584,NC_013270,NC_012837,NC_012647,NC_011049,NC_010232,NC_010233,NC_009886,NC_009422,NC_009423,NC_009258,NC_009264,NC_008410,NC_008144,NC_007888,NC_007440,NC_007178,NC_006839,NC_006688,NC_006689,NC_006690,NC_006402,NC_006403,NC_006405,NC_006406,NC_006408,NC_008975,NC_005794,NC_005055,NC_002805,NC_001573",
# 		"name": "Anura",
# 		"color": "blue"
# 	},
# ]
#
# groups = [
# 	{
# 		"ids": "AY521630,AM000053,AM000054,AY521629,AY521631,AM000055,DQ396400,AF004394,AF042101,AF256204,AF086817,AF042103,AB428555,AB287363,EU786678,U46016,AY713415,AY713416,AY255826,EF514713,DQ207941,DQ369994,EU786673,AY773338,AY773341,EF633445,AY322189,AF484516,AJ488927,AY371157,AY795907,U88826,AF061642,AY772535,AY586549,AF423760,AY371121,AB231893,EU786670",
# 		"name": "test",
# 		"color": "testcolor"
# 	}
# ]
#
#
# groups = [
# 	{
# 		"ids": "MN908947.3",
# 		"name": "test",
# 		"color": "testcolor"
# 	}
# ]
#
#
# groups = [
# 	{
# 		"ids": "A07867.1, AB289590.1, AB428555.1, AB565479.1, AB604946.1, AB731666.1, AF069140.1, AF070521.1, AF538302.1, AY037268.1, AY314049.1, AY314052.1, AY314061.1, AY331285.1, AY331297.1, AY423383.1, AY560108.1, AY779554.1, AY779561.1, AY779563.1, AY795905.1, AY835755.1, AY835772.1, AY970948.1, DQ127549.1, DQ295192.1, DQ322227.1, DQ676887.1, DQ853439.1, DQ853454.1, DQ854716.1, DQ886036.1, DQ990880.1, EF363122.1, EF514705.1, EF637046.1, EF637049.1, EU547186.1, EU786678.1, FJ460501.1, FJ469688.1, FJ469692.1, FJ469696.1, FJ469711.1, FJ469727.1, FJ469734.1, FJ469747.1, FJ469756.1, FJ469761.1, FJ496075.1, FJ496146.1, FJ496158.1, GQ386788.1, GQ386792.1, GU362883.1, JF320006.1, JF320010.1, JF320020.1, JF320043.1, JF320069.1, JF320082.1, JF320129.1, JF320146.1, JF320148.1, JF320158.1, JF320161.1, JF320191.1, JF320197.1, JF320204.1, JF320267.1, JF320420.1, JF320488.1, JF320534.1, JF320561.1, JF320565.1, JF320567.1, JF320628.1, JF683750.1, JF683784.1, JF683807.1, JF683809.1, JF689867.1, JF689887.1, JF932474.1, JF932486.1, JF932489.1, JF932497.1, JN024100.1, JN024135.1, JN024182.1, JN024208.1, JN024209.1, JN024234.1, JN024245.1, JN024253.1, JN024263.1, JN024266.1, JN024306.1, JN024344.1, JN024363.1, JN024375.1, JN024379.1, JN024385.1, JN024389.1, JN024397.1, JN024399.1, JN024402.1, JN024403.1, JN024477.1, JN024481.1, JN024487.1, JN024529.1, JN024531.1, JN024538.1, JN024547.1, JN024552.1, JN692446.1, JN692466.1, JN692468.1, JN692473.1, JN692475.1, JN944901.1, JN944931.1, JQ403024.1, JQ403030.1, JQ403054.1, JQ403055.1, JQ416159.1, JX446796.1, JX960598.1, JX960599.1, KC473828.1, KC473833.1, KC473842.1, KC596069.1, KC797175.1, KC797177.1, KC797229.1, KF384807.1, KF384812.1, KF526120.1, KF526122.1, KF526174.1, KF526298.1, KF526312.1, KJ140252.1, KJ140256.1, KJ140259.2, KJ704790.1, KJ704794.1, KJ849815.1, KJ849818.1, KJ948660.1, KP411824.1, KT124744.1, KT124748.1, KT124772.1, KT124796.1, KT124805.1, KT124808.1, KT200349.1, KT284378.1, KT427696.1, KT427699.1, KT427726.1, KT427730.1, KT427731.1, KT427751.1, KT427770.1, KT427797.1, KT427805.1, KU685589.1, KX505450.1, KX505536.1, U69593.1",
# 		"name": "test",
# 		"color": "testcolor"
# 	},
#
# 	{
# 		"ids": "AB254150.1, AF290030.1, AF443079.1, AF443095.1, AF443098.1, AX455929.1, AY043174.1, AY162223.1, AY253321.1, AY463226.1, AY463230.1, AY585268.1, AY713415.1, AY734556.1, DQ011166.1, DQ011179.1, DQ093585.1, DQ164111.1, DQ275655.1, DQ275664.1, DQ351221.1, DQ369977.1, DQ369978.1, DQ369995.1, DQ396366.1, DQ396374.1, DQ445637.1, EU293445.1, FJ496191.1, FJ496197.1, FJ496198.1, GQ999973.1, JX976690.1, JX976693.1, KC156117.1, KF250403.1, KF250404.1, KP109482.1, KP109483.1, KP109494.1, KP411831.1, KP411836.1, KR820294.1, KR820336.1, KR820344.1, KR820372.1, KR820396.1, KR820407.1, KR820410.1, KR820429.1, KR820431.1, KR820444.1, KT183052.1, KT183068.1, KT183123.1, KT183131.1, KT183142.1, KT183155.1, KT183172.1, KT183177.1, KT183188.1, KT183202.1, KT183212.1, KT183217.1, KT183221.1, KT183235.1, KT183243.1, KT183259.1, KT183266.1, KT183300.1, KT276258.1, KU319536.1, KU319541.1, KU319543.1, KU319546.1",
# 		"name": "test",
# 		"color": "testcolor"
# 	},
#
# 	{
# 		"ids": "AF193275.1, AF361872.1, AF457063.1, AF484508.1, AY829206.1, DQ207944.1, FJ388906.1, JQ292895.1, KF716472.1, KF716478.1, KP223766.1, KP223775.1, KP223780.1, KP223841.1, KP718928.1, KT022363.1, KT983615.1",
# 		"name": "test",
# 		"color": "testcolor"
# 	},
#
# 	{
# 		"ids": "A14116.1, AB485649.1, AF484480.1, AF484494.1, AY371155.1, KF716476.1",
# 		"name": "test",
# 		"color": "testcolor"
# 	},
#
# 	{
# 		"ids": "AF084936.1, FJ670530.1, KF716477.1, KX389645.1",
# 		"name": "test",
# 		"color": "testcolor"
# 	},
#
# 	{
# 		"ids": "AB485656.1, EU446022.1, KJ883138.1, KJ883139.1, KJ883142.1, KJ883145.1, KJ883148.1",
# 		"name": "test",
# 		"color": "testcolor"
# 	},
#
# 	{
# 		"ids": "EF029067.1",
# 		"name": "test",
# 		"color": "testcolor"
# 	},
# ]



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
	accession_numbers = [x.strip() for x in accession_numbers]
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

