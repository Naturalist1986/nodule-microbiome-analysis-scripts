from Bio import Phylo
from io import StringIO
import csv

# Load tree
with open("/mnt/c/Temp/Brady_tree.txt") as f:
    tree_str = f.read()

# Bio.Phylo struggles with bracket annotations like [0.981]; strip them
import re
tree_str_clean = re.sub(r'\[[\d\.]+\]', '', tree_str)

tree = Phylo.read(StringIO(tree_str_clean), "newick")

# Define clade groups from image
clade_groups = {
    "B. canariense": ["RH_binette_bin1", "hok_binette_bin4", "hok_binette_bin1"],
    "B. diazoefficiens 1": ["hok_binette_bin3"],
    "B. diazoefficiens 2": ["mtz_binette_bin3"],
    "B. retamae": ["ces_binette_bin1"],
    "B. algeriense": [
        "mtz_binette_bin10", "mtz_binette_bin12", "carR_binette_bin3",
        "carK_binette_bin6", "ces_binette_bin5", "carR_binette_bin2",
        "carK_binette_bin2", "RH_binette_bin2", "hok_binette_bin2"
    ],
}

all_mags = {mag for mags in clade_groups.values() for mag in mags}

# Verify all MAGs are in the tree
all_leaves = {c.name for c in tree.get_terminals()}
for mag in all_mags:
    if mag not in all_leaves:
        print(f"WARNING: {mag} not found in tree!")

# Build parent map
parents = {}
for clade in tree.find_clades(order="level"):
    for child in clade.clades:
        parents[child] = clade

# Keyword to search for in leaf names when climbing the tree for single-MAG clades
clade_keywords = {
    "B. diazoefficiens 1": "diazoefficiens",
    "B. diazoefficiens 2": "diazoefficiens",
    "B. retamae": "retamae",
}

def climb_until_keyword(tree, leaf_name, keyword):
    """Climb the tree from a leaf until the clade contains a leaf matching keyword."""
    node = next(tree.find_clades({"name": leaf_name}))
    while node in parents:
        node = parents[node]
        leaves = [c.name for c in node.get_terminals()]
        if any(keyword in (l or "") for l in leaves):
            return node
    return node  # return root if keyword never found

rows = []
for clade_name, mags in clade_groups.items():
    if len(mags) == 1 and clade_name in clade_keywords:
        keyword = clade_keywords[clade_name]
        mrca = climb_until_keyword(tree, mags[0], keyword)
    else:
        mrca = tree.common_ancestor(mags)

    leaves_in_clade = [c.name for c in mrca.get_terminals()]
    for leaf in leaves_in_clade:
        is_mag = "Yes" if leaf in all_mags else "No"
        rows.append({"Clade": clade_name, "Genome": leaf, "Is_MAG": is_mag})

with open("/mnt/c/Temp/clade_members.tsv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Clade", "Genome", "Is_MAG"], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

# Print summary
for clade_name in clade_groups:
    clade_rows = [r for r in rows if r["Clade"] == clade_name]
    mags_in = [r["Genome"] for r in clade_rows if r["Is_MAG"] == "Yes"]
    others = [r["Genome"] for r in clade_rows if r["Is_MAG"] == "No"]
    print(f"\n{clade_name}: {len(clade_rows)} total ({len(mags_in)} MAGs, {len(others)} references)")
    print(f"  References: {', '.join(others)}")

print("\nDone! Output: /mnt/c/Temp/clade_members.tsv")
