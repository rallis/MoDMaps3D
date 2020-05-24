#######################################################
## MoDMaps3D                                         ##
## MultiDimensional Scaling in Python                ##
## Coded by Rallis Karamichalis                      ##
## Github Repo: https://github.com/rallis/MoDMaps3D  ##
#######################################################

import json
import math, sys
import numpy as np

map_output = ""

print("Reading [mapheader.txt]\n")
header = open("mapheader.txt", "r")
for line in header:
	map_output += line
header.close()

print("Reading [accIDs.txt]\n")
accIDs = open("accIDs.txt", "r")
accessions = accIDs.readlines()
accIDs.close()
tmp = accessions[0].split('"')
accessions = [tmp[it] for it in range(1, len(tmp), 2)]

print("Reading [allSequences.txt]\n")
allSeq = open("allSequences.txt", "r")
allSequences = allSeq.readlines()
allSeq.close()
tmp = allSequences[0].split('"')
allSequences = [[tmp[it], tmp[it + 2], tmp[it + 4], tmp[it + 6]] for it in range(1, len(tmp), 8)]

# print("Reading [dMatrix.txt]\n")
# dMat = open ("dMatrix.txt", "r")
# distMat = dMat.readlines()
# dMat.close()
# tmp = distMat[0].split('{')
# dMat=[]
# for row in range(2,len(tmp)):
# 	dMat.append(tmp[row].split("}")[0].split(","))
# dMat = np.array(dMat, dtype="float64")

with open("dMatrix.txt", "r") as fin:
	dMat = json.loads(fin.read())

print(dMat)

# print("Computing eigenvalues..\n")
# eigValues, eigVectors = np.linalg.eig(dMat)
# idx = eigValues.argsort()[::-1][0:5]
# selEigValues = eigValues[idx]
# selEigVectors = eigVectors[:, idx]
# # print(idx
# # print(selEigValues.shape
# # print(selEigVectors.shape
# if False in (selEigValues > 0):
# 	print("First 5 largest eigenvalues are not all positive. Exiting..")
# 	sys.exit(-1)
#
# selEigVectors = np.array(selEigVectors)
#
# diagValues = []
# for i in range(len(selEigValues)):
# 	diagValues.append(math.sqrt(eigValues[i]))
# # print(diagValues
#
# diag = np.diag(diagValues)
# points = np.dot(selEigVectors, diag)
# # print("pointsSize=", points.shape
#
# minmaxScaling = []
# for i in range(5):
# 	minmaxScaling.append([min(points[:, i]), max(points[:, i])])
# # print(minmaxScaling
#
# scaledPoints = []
# for i in range(len(accessions)):
# 	scaledPoints.append([0, 0, 0, 0, 0])
# 	for j in range(5):
# 		scaledPoints[i][j] = 2.0 * (points[i][j] - minmaxScaling[j][0]) / (
# 				minmaxScaling[j][1] - minmaxScaling[j][0]) - 1
#
# print("Writing output [map_output.txt]\n")
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
#
# fout = open("map_output.txt", "w")
# fout.write(map_output)
# fout.close()

# from __future__ import absolute_import, division, unicode_literals
#
# import json
# import logging
# import numpy as np
# import scipy.sparse.linalg as linalg
#
#
# def mds(delta, dim):
#     delta = delta.astype(float)
#     (n, n) = delta.shape
#
#     deltasq = delta**2
#     deltatotals = np.sum(deltasq, axis=0)/n
#     sumOfDelta = np.sum(deltatotals)/n
#
#     # this way avoids temporaries (and is *much* faster than a loop)
#     bMatr = deltasq
#     bMatr -= deltatotals
#     bMatr = bMatr.transpose()
#     bMatr -= deltatotals
#     bMatr = bMatr.transpose()
#     bMatr += sumOfDelta
#     bMatr *= -0.5
#
#     (eigenvals, eigenvecs) = linalg.eigsh(bMatr, k=dim)
#     if (eigenvals < 0).any():
#         logging.getLogger('kameris.mds') \
#                .warning('some eigenvalues were negative')
#
#     # not sure why eigensystem is sorted in reverse order for eigsh...
#     points = np.fliplr(np.dot(eigenvecs, np.diag(np.sqrt(eigenvals))))
#     # TODO: see why the NaNs sometimes happen (maybe run dim+1?)
#     return points[:, ~np.isnan(points).any(axis=0)]
#
#
# def run_mds_step(options, exp_options):
#     dists = kameris_formats.dist_reader.read_matrix(options['dists_file'])
#     points = mds(dists, options['dimensions']).tolist()
#
#     with open(options['output_file'], 'w') as outfile:
#         json.dump(points, outfile)
