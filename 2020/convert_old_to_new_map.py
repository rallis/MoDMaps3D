import os, glob, json

root_dir = os.getcwd()
os.chdir("..")
os.chdir("maps")
os.chdir("refmaps")
print(os.getcwd())

files = glob.glob("*.txt")
print(len(files))

for fl in files:
    print(fl)
    with open(fl, "r", encoding="utf8") as fin:
        data = json.loads(fin.read())

    print(data)

    new_map = {
        "description": data["mapDescription"],
        "kmer": 9,
        "enable_taxa_info": True,
        "groups": []
    }
    names = data["namesets"]
    colors = data["colorsets"]
    sets = data["sets"]

    for idx in range(int(data["numofsets"])):
        new_map["groups"].append({
            "name": names[idx],
            "color": colors[idx],
            "accession_numbers": sets[idx].split(", "),
        })

    print(new_map)

    with open(fl[:-3] + "json", "w") as fout:
        fout.write(json.dumps(new_map, indent=4))
