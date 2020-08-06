#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import math
import sys
import numpy as np
import scipy.linalg

with open("dMatrix.json", "r") as fin:
    dMatrix = json.loads(fin.read())

with open("dataset.json", "r") as fin:
    dataset = json.loads(fin.read())

with open("all_sequences_info.json", "r") as fin:
    all_sequences_info = json.loads(fin.read())


print("Computing eigenvalues..")
# eigValues, eigVectors = np.linalg.eig(dMatrix)
eigValues, eigVectors = scipy.linalg.eigh(dMatrix)
idx = eigValues.argsort()[::-1][0:5]
selEigValues = eigValues[idx]
selEigVectors = eigVectors[:,idx]

if False in (selEigValues > 0):
    print("First 5 largest eigenvalues are not all positive. Exiting..")
    sys.exit(-1)

selEigVectors = np.array(selEigVectors)

diagValues = []
for i in range(len(selEigValues)):
    diagValues.append(math.sqrt(selEigValues[i]))
# print(diagValues)

diag = np.diag(diagValues)
points = np.dot(selEigVectors, diag)
# print("pointsSize=", points.shape)

minmaxScaling = []
for i in range(5):
    minmaxScaling.append([min(points[:, i]), max(points[:, i])])
# print(minmaxScaling)

scaledPoints = []
for i in range(len(all_sequences_info)):
    scaledPoints.append([0, 0, 0, 0, 0])
    for j in range(5):
        scaledPoints[i][j] = 2.0 * (points[i][j] - minmaxScaling[j][0]) / (minmaxScaling[j][1] - minmaxScaling[j][0]) - 1
# print(scaledPoints)

json_mapfile = {
    "description": dataset["description"],
    "kmer": dataset["kmer"],
    "enable_taxa_info": dataset["enable_taxa_info"],
    "groups": dataset["groups"],
    "points": scaledPoints,
    "all_sequences_info": all_sequences_info,
}

print("Exporting MoDMap..")
with open("json_mapfile.json", "w", encoding="utf8") as fout:
    fout.write(json.dumps(json_mapfile, indent=4))
