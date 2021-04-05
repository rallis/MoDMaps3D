

import os, json, requests, time

root_dir = os.getcwd()

url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="


os.chdir(os.path.join(root_dir,"maps","samplemaps"))

with open("sample3.json", "r", encoding="utf-8") as fin:
    data = json.load(fin)

# print(data)
c = 0
groups = data["groups"]
for gr in groups:

    accession_numbers = gr["accession_numbers"]

    for acc in accession_numbers:
        print(c, acc)
        c += 1

        # if c < 2522:
        #     continue

        while True:
            response = requests.get(url + acc)
            # print(response.text)
            has_error = "error" in response.text
            if has_error:
                print("To sleep because: ", response.text)
                time.sleep(1)
            else:
                break
        with open(os.path.join(root_dir, "2020", "data_repo", acc + ".fasta"), 'w') as fout:
            fout.write(response.text)

# NC_008964