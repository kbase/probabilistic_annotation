#!/usr/bin/python

import json

def GenomeJsonToFasta(genome_json_file, out_file):
    # Read the existing annotation_file.
    resp = json.load(open(genome_json_file, "r"))
    fout = open(out_file, "w")
    features = resp["features"]
    for feature in features:
        # Not a protein-encoding gene
        if "protein_translation" not in feature:
            continue
        myid = feature["id"]
        if "function" in feature:
            function = feature["function"]
        else:
            function = ""
        seq = feature["protein_translation"]
        fout.write(">%s %s\n%s\n" %(myid, function, seq))
    return

GenomeJsonToFasta("TEST_GENOME_JSON", "kb|g.22438.faa")
