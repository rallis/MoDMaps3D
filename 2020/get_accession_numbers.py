#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
import shutil
import requests
import time


def download_file(url, local_filename):
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    return local_filename


root_dir = os.getcwd()
accession_numbers = []
data_dir = "fasta"


# EXTRACT ACCESSION NUMBERS

with open("nuccore_result.txt", "r", encoding="utf8") as fin:
    lines = fin.readlines()
    print(len(lines))

for idx in range(3, len(lines), 4):
    parts = lines[idx].split(" ")
    accession_numbers.append(parts[0])
print(len(accession_numbers))

with open("accession_numbers.json", "w", encoding="utf8") as fout:
    fout.write(json.dumps(accession_numbers))

#
# # DOWNLOAD FASTA FILES
#
# if not os.path.exists(data_dir):
#     os.makedirs(data_dir)
# # else:
# #     shutil.rmtree(data_dir)
# #     os.makedirs(data_dir)
# os.chdir(data_dir)
#
# for idx, acc_num in enumerate(accession_numbers[:5]):
#     # if idx < 5230:
#     #     continue
#
#     print(idx, acc_num, len(accession_numbers))
#     # get_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={}".format(acc_num)
#     get_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=xml&id={}".format(acc_num)
#     local_filename = "{}.fasta".format(acc_num)
#
#     # download_file(get_url, local_filename)
#
#     r = requests.get(get_url)
#     resp = r.text
#     with open(local_filename, "w", encoding="utf8") as fout:
#         fout.write(resp)
#
#     try:
#         json_response = r.json()
#         json_error = (json_response.get("error") == "API rate limit exceeded")
#     except Exception as ex:
#         json_error = False
#
#     if json_error:
#         raise Exception("API limit, over 3TPS")
#
#     if idx % 3 == 0:
#         time.sleep(1.3)


#
# import glob
# os.chdir(data_dir)
# files = glob.glob("*.fasta")
# for fl in files:
#     print(fl.split(".")[0])
#     os.rename(fl, fl.split(".")[0])