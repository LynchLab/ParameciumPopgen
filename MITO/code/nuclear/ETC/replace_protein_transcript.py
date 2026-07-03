#!/usr/bin/env python3

import re



# --------------------------------------------------

# 1) Standard .P → .T replacements

# --------------------------------------------------

PREFIXES = [

    "PBIA.V1_4.1",

    "PBUR.Dd1.1",

    "PCAU.43c3d.1",

    "PDEC.223.1",

    "PDODEC.274.1",

    "PJENN.M.1",

    "PNOV.TE.1",

    "POCT.K8.1",

    "PPENT.87.1",

    "PQUADEC.NiA.1",

    "PSEX.AZ8_4.1",

    "PSON.ATCC_30995.1",

    "PTET.51.1",

    "PTRED.209.2",

]



re_dotP = re.compile(

    r'(' + '|'.join(map(re.escape, PREFIXES)) + r')\.P(\d+)'

)



# --------------------------------------------------

# 2) Compact protein → transcript prefixes

# --------------------------------------------------

re_pmmnp = re.compile(r'\bPMMNP(\d+)\b')    # multimicronucleatum

re_pprimp = re.compile(r'\bPPRIMP(\d+)\b') # primaurelia



def normalize_ids(text: str) -> str:

    """

    Normalize all Paramecium protein IDs to transcript IDs.

    """

    # .P → .T

    text = re_dotP.sub(r'\1.T\2', text)



    # compact prefixes

    text = re_pmmnp.sub(r'PMMNT\1', text)

    text = re_pprimp.sub(r'PPRIMT\1', text)



    return text





# --------------------------------------------------

# File runner

# --------------------------------------------------

if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:

        print("Usage: fix_ids.py input.tsv output.tsv")

        sys.exit(1)



    infile, outfile = sys.argv[1], sys.argv[2]



    with open(infile) as fin, open(outfile, "w") as fout:

        for line in fin:

            fout.write(normalize_ids(line))


