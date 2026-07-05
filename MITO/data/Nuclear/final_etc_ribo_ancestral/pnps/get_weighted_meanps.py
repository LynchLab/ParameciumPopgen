import csv

species = [
    "biaurelia",
    "caudatum",
    "multimicronucleatum",
    "primaurelia",
    "sexaurelia",
    "tetraurelia",
    "tribursaria"
]

# ETC mean_piS
etc_ps = [
    0.0062534960234979725,
    0.03923678125343537,
    0.05410180763945808,
    0.04094593660577718,
    0.021016389182603505,
    0.007336328670063226,
    0.04467702961362838
]

etc_n = [52, 66, 28, 60, 61, 52, 65]

# RIBO mean_piS
ribo_ps = [
    0.007587898411690726,
    0.052935514859021426,
    0.036867541901551525,
    0.03514380488256197,
    0.02276194725597661,
    0.007334159171250127,
    0.051968934185386534
]

ribo_n = [66, 51, 40, 73, 77, 64, 45]

weighted_ps = []

for e_ps, r_ps, e_n, r_n in zip(etc_ps, ribo_ps, etc_n, ribo_n):
    w = (e_ps*e_n + r_ps*r_n) / (e_n + r_n)
    weighted_ps.append(w)

for s, w in zip(species, weighted_ps):
    print(s, w)

# save CSV
with open("weighted_mean_pS_nuclear_ETC_RIBO_7species.csv","w",newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["species","mean_pS"])
    for s, w in zip(species, weighted_ps):
        writer.writerow([s, w])
