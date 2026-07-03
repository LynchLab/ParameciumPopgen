from pathlib import Path
import shutil
import pandas as pd

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
SRC = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/dnds/etc_ribo/M8_vs_M8a_icode5")

DEST = Path("/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/final_etc_ribo/dnds/etc/M8_i3")

ETC_LIST = SRC / "etc_genes.txt"
M8_CSV = SRC / "M8_beta_pq_omega.csv"

# ------------------------------------------------------------
# Make destination folder
# ------------------------------------------------------------
DEST.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Read ETC gene list
# ------------------------------------------------------------
with open(ETC_LIST) as f:
    etc_genes = [line.strip() for line in f if line.strip()]

etc_set = set(etc_genes)

print(f"ETC genes in list: {len(etc_genes)}")

# ------------------------------------------------------------
# Copy only folders listed in etc_genes.txt
# ------------------------------------------------------------
copied = []
missing_folders = []

for gene in etc_genes:
    src_folder = SRC / gene
    dest_folder = DEST / gene

    if src_folder.is_dir():
        if dest_folder.exists():
            shutil.rmtree(dest_folder)
        shutil.copytree(src_folder, dest_folder)
        copied.append(gene)
        print(f"copied: {gene}")
    else:
        missing_folders.append(gene)
        print(f"missing folder: {gene}")

# ------------------------------------------------------------
# Copy etc_genes.txt
# ------------------------------------------------------------
shutil.copy2(ETC_LIST, DEST / "etc_genes.txt")

# ------------------------------------------------------------
# Filter M8_beta_pq_omega.csv to only ETC genes
# ------------------------------------------------------------
df = pd.read_csv(M8_CSV)

df_etc = df[df["gene"].isin(etc_set)].copy()

df_etc.to_csv(DEST / "M8_beta_pq_omega.csv", index=False)

# ------------------------------------------------------------
# Report removed genes from CSV
# ------------------------------------------------------------
csv_genes = set(df["gene"].astype(str))
kept_genes = set(df_etc["gene"].astype(str))

removed_from_csv = sorted(csv_genes - etc_set)
missing_from_csv = sorted(etc_set - csv_genes)

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------
print("\nSummary")
print("-------")
print(f"Destination: {DEST}")
print(f"Folders copied: {len(copied)}")
print(f"Missing folders: {len(missing_folders)}")
print(f"Rows in original M8 csv: {len(df)}")
print(f"Rows kept in filtered M8 csv: {len(df_etc)}")
print(f"Rows removed from M8 csv: {len(removed_from_csv)}")

if missing_folders:
    print("\nGenes in etc_genes.txt but folder missing:")
    for g in missing_folders:
        print(g)

if missing_from_csv:
    print("\nGenes in etc_genes.txt but missing from M8_beta_pq_omega.csv:")
    for g in missing_from_csv:
        print(g)

if removed_from_csv:
    print("\nGenes removed from M8_beta_pq_omega.csv because not in etc_genes.txt:")
    for g in removed_from_csv:
        print(g)

print("\nDone.")
