// HFE.cpp: haploid frequency estimator.
// Site-by-site ML estimation of the major-allele frequency (p) and the error rate (e)
// from individual nucleotide read quartets of haploid / fully homozygous individuals,
// assuming at most two alleles per site in the population. 
// Uses the GFE input format (Maruki & Lynch 2015, Package-GFE) and a reduced GFE output layout.
//
// Model, individual k with n_M,k major, n_m,k minor and n_e,k other reads:
//   L_k = p (1-e)^n_M (e/3)^(n_m+n_e) + (1-p) (1-e)^n_m (e/3)^(n_M+n_e)
// 
// Maximised by EM:
//   E-step  g_k = 1 / (1 + exp(-d_k)),  d_k = logit(p) + (n_M,k - n_m,k) ln[(1-e)/(e/3)]
//   M-step  p = mean_k g_k
//           e = sum_k [g_k (n_m,k + n_e,k) + (1-g_k)(n_M,k + n_e,k)] / N
// 
// Start: pooled-read closed form p0 = (N_M - N_e/2)/(N_M + N_m - N_e), e0 = 1.5 N_e / N.
// 
// Null model (monomorphic, all individuals carry M): e_null = (N - N_M)/N.
// 
// pol_llstat = 2 (lnL_full - lnL_null). p = 1 is on the boundary, so under the null the
// statistic is asymptotically a 50:50 mixture of 0 and chi-square(1); 5% critical value 2.706.
//
// Build:  g++ -O3 -std=c++17 -o HFE HFE.cpp
//
// Usage:  ./HFE [-in In_GFE.txt] [-out Out_HFE.txt] [-min_cov 1] [-max_cov 2000000000]
//               [-min_ind 1] [-max_e 1.0] [-mode f|F] [-pop pop] [-ref_info 1|0]
//   -min_cov / -max_cov  individuals outside this coverage range are ignored at the site
//   -min_ind             minimum number of individuals left after the coverage filter (Nr)
//   -max_e               maximum allowed error-rate estimate
//   -mode                f (default): full HFE output described below;
//                        F: input layout of Maruki's PSA (as GFE v3.0 -mode F), see below
//   -pop                 population label prefixed to the column names in the F mode
//   -ref_info            F mode only: 1 (default) prints scaffold, site, ref_nuc first; 0 omits them
//
// F mode, for site-wise W&C (1984) fixation indices with PSA. Run HFE on each population with
// the same site list, the first with -ref_info 1 and the rest with -ref_info 0, then
//   paste Out_F_pop1.txt Out_F_pop2.txt ... > In_PSA.txt
//   ./PSA -in In_PSA.txt -out Out_PSA.txt -cv 2.706 -min_Ni <min haploid individuals>
// Nine columns per population, in PSA's order: major_allele, minor_allele, pop_coverage,
//   Ni, best_Maf (p), best_maf (q), best_error, best_H, pol_llstat, where
//   Ni = Nr (a haploid individual with >= 1 read contributes one fully sampled genome; not
//        GFE's diploid sum of 1 - 0.5^cov), and best_H = 0 (no heterozygotes). With H = 0,
//   PSA's theta is W&C's random-union-of-gametes estimator with n_i = Nr; its F and f
//   columns are identically 1 and meaningless for haploids.
// 
//
// Sites with no coverage, with Nr < min_ind, or with best_error > max_e are printed with NA
// estimates, as in GFE. 
// Sites with a single observed nucleotide, or where the null model fits at least as well,
// are reported with p = 1, minor_allele NA and pol_llstat 0.
//
// DEFAULT output columns: scaffold, site, ref_nuc, major_allele, minor_allele, pop_coverage,
//   Nr (individuals passing the coverage filter), best_Mac, best_mac (major/minor allele
//   counts, Nr*p rounded), best_p, best_q, best_error, null_best_error, pol_llstat, n_iter.


#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static const char NUC[4] = {'A', 'C', 'G', 'T'};
static const double EMIN = 1e-12;  // floor on e inside logs

struct Ind { int nM, nm, ne; };

// log-likelihood of all individuals at (p, e); p in [0,1]
static double loglik(const vector<Ind> &I, double p, double e) {
  double ee = max(e, EMIN);
  double l1e = log1p(-ee), le3 = log(ee / 3);
  double lp = p > 0 ? log(p) : -INFINITY, lq = p < 1 ? log1p(-p) : -INFINITY;
  double ll = 0;
  for (const Ind &x : I) {
    double a = lp + x.nM * l1e + (x.nm + x.ne) * le3;
    double b = lq + x.nm * l1e + (x.nM + x.ne) * le3;
    double mx = max(a, b);
    ll += mx + log1p(exp(min(a, b) - mx));
  }
  return ll;
}

// EM; returns iterations used; p, e updated in place, ll set to final log-likelihood
static int em(const vector<Ind> &I, long N, double &p, double &e, double &ll,
              double tol = 1e-10, int maxit = 10000) {
  ll = loglik(I, p, e);
  int it = 0;
  while (it < maxit) {
    it++;
    double ee = max(e, EMIN);
    double lr = log1p(-ee) - log(ee / 3);
    double lgt = log(p) - log1p(-p);
    double sg = 0, R = 0;
    for (const Ind &x : I) {
      double d = lgt + (x.nM - x.nm) * lr;
      double g = 1 / (1 + exp(-d));
      sg += g;
      R += g * (x.nm + x.ne) + (1 - g) * (x.nM + x.ne);
    }
    double pn = sg / I.size(), en = R / N;
    // keep p strictly inside (0,1) so logit stays finite; the null model covers p = 1
    pn = min(max(pn, 1e-12), 1 - 1e-12);
    double lln = loglik(I, pn, en);
    double dl = lln - ll;
    p = pn; e = en; ll = lln;
    if (fabs(dl) <= tol * (fabs(ll) + tol)) break;
  }
  return it;
}

static int parse_quartet(const string &s, array<int, 4> &c) {
  int j = 0; long v = 0; bool any = false;
  for (char ch : s) {
    if (ch == '/') { if (j > 3) return 0; c[j++] = (int)v; v = 0; any = false; }
    else if (ch >= '0' && ch <= '9') { v = v * 10 + (ch - '0'); any = true; }
    else return 0;
  }
  if (j != 3 || !any) return 0;
  c[3] = (int)v;
  return 1;
}

int main(int argc, char *argv[]) {
  const char *in_name = "In_GFE.txt", *out_name = "Out_HFE.txt";
  int min_cov = 1, max_cov = 2000000000, min_ind = 1;
  double max_e = 1.0;
  string mode = "f", pop_label = "pop";
  int ref_info = 1;
  bool help = false;
  for (int a = 1; a < argc; a++) {
    string s = argv[a];
    bool hasval = a + 1 < argc;
    if (s == "-h") help = true;
    else if (s == "-in" && hasval) in_name = argv[++a];
    else if (s == "-out" && hasval) out_name = argv[++a];
    else if (s == "-min_cov" && hasval) min_cov = atoi(argv[++a]);
    else if (s == "-max_cov" && hasval) max_cov = atoi(argv[++a]);
    else if (s == "-min_ind" && hasval) min_ind = atoi(argv[++a]);
    else if (s == "-max_e" && hasval) max_e = atof(argv[++a]);
    else if (s == "-mode" && hasval) mode = argv[++a];
    else if (s == "-pop" && hasval) pop_label = argv[++a];
    else if (s == "-ref_info" && hasval) ref_info = atoi(argv[++a]);
    else { fprintf(stderr, "unknown option or missing value: %s\n", s.c_str()); help = true; break; }
  }
  if (!help && mode != "f" && mode != "F") { fprintf(stderr, "unknown mode: %s\n", mode.c_str()); help = true; }
  if (!help && ref_info != 0 && ref_info != 1) { fprintf(stderr, "-ref_info must be 0 or 1\n"); help = true; }
  if (help) {
    fprintf(stderr, "USAGE: %s {<options>}\n", argv[0]);
    fprintf(stderr, "\t-h: print this message\n");
    fprintf(stderr, "\t-in <s>: input file (default In_GFE.txt)\n");
    fprintf(stderr, "\t-out <s>: output file (default Out_HFE.txt)\n");
    fprintf(stderr, "\t-min_cov <d>: minimum coverage required for each individual (default 1)\n");
    fprintf(stderr, "\t-max_cov <d>: maximum coverage allowed for each individual (default 2000000000)\n");
    fprintf(stderr, "\t-min_ind <d>: minimum number of individuals passing the coverage filter (default 1)\n");
    fprintf(stderr, "\t-max_e <f>: maximum allowed error-rate estimate (default 1.0)\n");
    fprintf(stderr, "\t-mode <s>: f (default, full output) or F (PSA input layout)\n");
    fprintf(stderr, "\t-pop <s>: population label for the column names in the F mode (default pop)\n");
    fprintf(stderr, "\t-ref_info <d>: F mode: 1 prints scaffold, site, ref_nuc; 0 omits them (default 1)\n");
    return 1;
  }

  ifstream in(in_name);
  if (!in.is_open()) { fprintf(stderr, "Cannot open %s for reading.\n", in_name); return 1; }
  FILE *out = fopen(out_name, "w");
  if (!out) { fprintf(stderr, "Cannot open %s for writing.\n", out_name); return 1; }

  string line;
  getline(in, line);
  int nsample = 0;
  {
    istringstream ss(line);
    string tok;
    int i = 0;
    while (ss >> tok) if (i++ >= 3) nsample++;
  }
  printf("%d individuals to be analyzed\n", nsample);

  const bool Fm = mode == "F";
  if (!Fm) {
    fprintf(out, "scaffold\tsite\tref_nuc\tmajor_allele\tminor_allele\tpop_coverage\tNr\tbest_Mac\tbest_mac\t"
                 "best_p\tbest_q\tbest_error\tnull_best_error\tpol_llstat\tn_iter\n");
  } else {
    const char *P = pop_label.c_str();
    if (ref_info) fprintf(out, "scaffold\tsite\tref_nuc\t");
    fprintf(out, "%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Maf\t%s_best_maf\t"
                 "%s_best_error\t%s_best_H\t%s_pol_llstat\n", P, P, P, P, P, P, P, P, P);
  }

  vector<array<int, 4>> rc(nsample);
  vector<Ind> I;
  long lineno = 1;
  while (getline(in, line)) {
    lineno++;
    if (line.empty()) continue;
    istringstream ss(line);
    string scaffold, site, ref, q;
    ss >> scaffold >> site >> ref;
    array<long, 4> pop{0, 0, 0, 0};
    int Nr = 0;
    bool bad = false;
    for (int k = 0; k < nsample; k++) {
      if (!(ss >> q) || !parse_quartet(q, rc[k])) { bad = true; break; }
      int cov = rc[k][0] + rc[k][1] + rc[k][2] + rc[k][3];
      if (cov < min_cov || cov > max_cov) cov = 0, rc[k] = {0, 0, 0, 0};
      if (cov > 0) {
        Nr++;
        for (int j = 0; j < 4; j++) pop[j] += rc[k][j];
      }
    }
    // F mode: leading reference columns only with -ref_info 1
    auto print_ref = [&]() {
      if (ref_info) fprintf(out, "%s\t%s\t%s\t", scaffold.c_str(), site.c_str(), ref.c_str());
    };
    if (bad) {
      if (!Fm) { fprintf(stderr, "Skipping line %ld: bad or missing read quartet\n", lineno); continue; }
      // keep pasted files aligned; pop_coverage must stay numeric for PSA
      fprintf(stderr, "Line %ld: bad or missing read quartet, printed as NA\n", lineno);
      print_ref();
      fprintf(out, "NA\tNA\t0\t0\tNA\tNA\tNA\tNA\tNA\n");
      continue;
    }
    long N = pop[0] + pop[1] + pop[2] + pop[3];

    auto print_na = [&]() {
      if (Fm) {
        print_ref();
        fprintf(out, "NA\tNA\t%ld\t%d\tNA\tNA\tNA\tNA\tNA\n", N, Nr);
      } else {
        fprintf(out, "%s\t%s\t%s\tNA\tNA\t%ld\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n",
                scaffold.c_str(), site.c_str(), ref.c_str(), N, Nr);
      }
    };
    if (N == 0 || Nr < min_ind) { print_na(); continue; }

    // candidate alleles: two most abundant nucleotides in the pooled sample (ties: A<C<G<T order)
    array<int, 4> ord{0, 1, 2, 3};
    stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return pop[a] > pop[b]; });
    int n1 = ord[0], n2 = ord[1];

    if (pop[n1] == N) {  // only one nucleotide observed: monomorphic, no errors
      if (Fm) {
        print_ref();
        fprintf(out, "%c\tNA\t%ld\t%d\t1.000000\t0.000000\t0.000000\t0.000000\tNA\n", NUC[n1], N, Nr);
      } else {
        fprintf(out, "%s\t%s\t%s\t%c\tNA\t%ld\t%d\t%d\t0\t1.000000\t0.000000\t0\t0\t0.000000\t0\n",
                scaffold.c_str(), site.c_str(), ref.c_str(), NUC[n1], N, Nr, Nr);
      }
      continue;
    }
    if (pop[n2] > 0 && pop[ord[2]] == pop[n2])
      fprintf(stderr, "The identity of the minor allele is ambiguous at site %s on %s\n", site.c_str(), scaffold.c_str());

    I.clear();
    long NM = 0, Nm = 0, Ne = 0;
    for (int k = 0; k < nsample; k++) {
      int cov = rc[k][0] + rc[k][1] + rc[k][2] + rc[k][3];
      if (cov == 0) continue;
      Ind x{rc[k][n1], rc[k][n2], cov - rc[k][n1] - rc[k][n2]};
      NM += x.nM; Nm += x.nm; Ne += x.ne;
      I.push_back(x);
    }

    // null model: all individuals carry the major allele
    double null_e = (double)(N - NM) / N;
    double null_ll = loglik(I, 1.0, null_e);

    // full model by EM from the pooled-read closed form
    double den = (double)(NM + Nm - Ne);
    double p = den > 0 ? (NM - 0.5 * Ne) / den : 0.5;
    p = min(max(p, 0.01), 0.99);
    double e = min(max(1.5 * Ne / N, 1e-4), 0.25);
    double ll;
    int iters = em(I, N, p, e, ll);

    if (null_ll >= ll) { p = 1.0; e = null_e; ll = null_ll; }
    double pol = 2 * (ll - null_ll);

    if (e > max_e) { print_na(); continue; }

    // report the more frequent of the two as major; no minor allele when p = 1
    int maj = n1, mnr = n2;
    if (p < 0.5) { p = 1 - p; swap(maj, mnr); }
    int Mac = (int)lround(Nr * p), mac = Nr - Mac;
    string minor = p < 1.0 ? string(1, NUC[mnr]) : "NA";
    if (Fm) {
      // Ni = Nr; best_H = 0 under the haploid model
      print_ref();
      fprintf(out, "%c\t%s\t%ld\t%d\t%f\t%f\t%.6g\t0.000000\t%f\n",
              NUC[maj], minor.c_str(), N, Nr, p, 1 - p, e, pol);
    } else {
      fprintf(out, "%s\t%s\t%s\t%c\t%s\t%ld\t%d\t%d\t%d\t%f\t%f\t%.6g\t%.6g\t%f\t%d\n",
              scaffold.c_str(), site.c_str(), ref.c_str(), NUC[maj], minor.c_str(), N, Nr, Mac, mac,
              p, 1 - p, e, null_e, pol, iters);
    }
  }
  fclose(out);
  return 0;
}
