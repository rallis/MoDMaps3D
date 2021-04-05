

import os, json, requests, time

root_dir = os.getcwd()

# url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="
url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=xml&id="

# extension_to_use = ".fasta"
extension_to_use = ".gb"


os.chdir(os.path.join(root_dir,"2020"))

with open("sequence.txt", "r", encoding="utf-8") as fin:
    data = fin.readlines()

accession_numbers = [x.strip() for x in data]
print(len(accession_numbers))


c = 0
for acc in accession_numbers:
    acc_without_version = acc.split(".")[0]
    print(c, acc)
    c += 1

    if c < 22670:
        continue

    while True:
        response = requests.get(url + acc)
        # print(response.text)
        has_error = "error" in response.text
        if has_error:
            print("To sleep because: ", response.text)
            time.sleep(1)
        else:
            break
    with open(os.path.join(root_dir, "2020", "data_mtdna_genbank", acc_without_version + extension_to_use), 'w') as fout:
        fout.write(response.text)

# NC_008964