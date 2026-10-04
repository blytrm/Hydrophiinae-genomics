#!/usr/bin/env python3
"""ene stats via gffutils. Usage: gffutils_script.py <gff3>"""
import gffutils
import sys

gff = sys.argv[1]
db = gffutils.create_db(gff, dbfn=":memory:", force=True, keep_order=True,
                        merge_strategy="create_unique", sort_attribute_values=True)

genes = list(db.features_of_type("gene"))
cds = list(db.features_of_type("CDS"))

# genes lacking any descriptive attribute
named_keys = ("Name", "product", "description")
unknown = [g for g in genes if not any(k in g.attributes for k in named_keys)]

lengths = [g.end - g.start + 1 for g in genes] or [0]

print(f"gene_count\t{len(genes)}")
print(f"cds_count\t{len(cds)}")
print(f"gene_unknown_count\t{len(unknown)}")
print(f"gene_length_min\t{min(lengths)}")
print(f"gene_length_max\t{max(lengths)}")
print(f"gene_length_mean\t{sum(lengths) / len(lengths):.1f}")
