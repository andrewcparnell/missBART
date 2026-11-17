// missBART: C++ samplers for multivariate BART with (possibly non-ignorable)
// missingness in the response.
//
// One sampler, three models:
//   model 0: multivariate BART on complete responses (mvBART)
//   model 1: multivariate BART + multivariate probit regression (missBART1)
//   model 2: multivariate BART + multivariate probit BART (missBART2)
//
// All random numbers come from R's RNG, so set.seed() gives reproducible fits.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <map>
#include <cmath>

using namespace Rcpp;

// ----------------------------------------------------------------------------
// Small sampling helpers
// ----------------------------------------------------------------------------

// Standard normal truncated to (a, Inf). Robert (1995).
static inline double rtn_tail(double a) {
  double x;
  if (a <= 0.45) {
    do { x = norm_rand(); } while (x <= a);
    return x;
  }
  const double lam = 0.5 * (a + std::sqrt(a * a + 4.0));
  do {
    x = a + exp_rand() / lam;
  } while (unif_rand() > std::exp(-0.5 * (x - lam) * (x - lam)));
  return x;
}

// N(mu, sd^2) truncated to (0, Inf) if positive, (-Inf, 0] otherwise.
static inline double rtn_sign(double mu, double sd, bool positive) {
  return positive ? mu + sd * rtn_tail(-mu / sd) : mu - sd * rtn_tail(mu / sd);
}

static bool safe_chol(arma::mat& out, const arma::mat& A, const char* layout) {
  arma::mat S = 0.5 * (A + A.t());
  if (arma::chol(out, S, layout)) return true;
  double jitter = 1e-10 * std::max(1.0, arma::trace(S) / S.n_rows);
  for (int k = 0; k < 12; ++k) {
    arma::mat Sj = S;
    Sj.diag() += jitter;
    if (arma::chol(out, Sj, layout)) return true;
    jitter *= 10.0;
  }
  return false;
}

static arma::mat safe_inv(const arma::mat& A) {
  arma::mat U;
  if (!safe_chol(U, A, "upper")) stop("matrix is not positive definite");
  arma::mat Ui = arma::inv(arma::trimatu(U));
  return Ui * Ui.t();
}

// Wishart(df, S) by the Bartlett decomposition.
static arma::mat rwish(double df, const arma::mat& S) {
  const int p = S.n_rows;
  arma::mat L;
  if (!safe_chol(L, S, "lower")) stop("Wishart scale is not positive definite");
  arma::mat A(p, p, arma::fill::zeros);
  for (int i = 0; i < p; ++i) {
    A(i, i) = std::sqrt(R::rchisq(df - i));
    for (int j = 0; j < i; ++j) A(i, j) = norm_rand();
  }
  arma::mat LA = L * A;
  return LA * LA.t();
}

// Inverse-Wishart(df, S): the inverse of a Wishart(df, S^{-1}) draw.
static arma::mat riwish(double df, const arma::mat& S) {
  return safe_inv(rwish(df, safe_inv(S)));
}

// ----------------------------------------------------------------------------
// Trees
// ----------------------------------------------------------------------------

struct Node {
  int parent;    // -1 for the root
  int left;      // -1 for a terminal node; right child is always left + 1
  int var;       // splitting variable (internal nodes)
  double val;    // go left if x <= val
  bool na_left;  // direction taken by missing covariate values
  int depth;
};

struct Tree {
  std::vector<Node> nodes;
  arma::mat mu;            // p x nodes.size(); only terminal columns are used
  std::vector<int> leaf;   // terminal node of each observation
};

// Terminal node reached by one observation. x points at the first covariate of
// the observation and consecutive covariates are `stride` apart.
static inline int descend(const std::vector<Node>& nd, const double* x,
                          std::size_t stride) {
  int k = 0;
  while (nd[k].left >= 0) {
    const double v = x[stride * nd[k].var];
    const bool go_left = std::isnan(v) ? nd[k].na_left : (v <= nd[k].val);
    k = go_left ? nd[k].left : nd[k].left + 1;
  }
  return k;
}

// Tree-dependent part of the log marginal likelihood of the partial residuals
// in one terminal node holding `cnt` observations with column sums `s`:
// (mu_r' Sigma_r^{-1} mu_r + log|Sigma_r|) / 2 in the notation of the paper,
// using the eigendecomposition U diag(lam) U' of the residual precision.
static inline double leaf_loglik(int cnt, const double* s, int p, double kappa,
                                 const arma::mat& U, const arma::vec& lam,
                                 const arma::vec& kUmu0) {
  double val = 0.0;
  for (int k = 0; k < p; ++k) {
    const double d = cnt * lam[k] + kappa;
    const double* u = U.colptr(k);
    double us = 0.0;
    for (int j = 0; j < p; ++j) us += u[j] * s[j];
    const double c = lam[k] * us + kUmu0[k];
    val += c * c / d - std::log(d);
  }
  return 0.5 * val;
}

class Forest {
public:
  int n, p, q, K;
  const arma::mat& X;              // n x q, may change between sweeps
  std::vector<Tree> trees;
  arma::mat fit;                   // p x n, sum over trees
  double kappa;
  arma::vec mu0;
  double prior_alpha, prior_beta;
  int min_node, max_attempt;
  bool check_rows;                 // probit trees: no node made only of
  std::vector<char> row_ok;        //   rows flagged as not ok
  arma::vec xmin, xmax;
  double n_prop, n_acc;

  Forest(const arma::mat& X_, int p_, int K_, double kappa_, const arma::vec& mu0_,
         double a_, double b_, int min_node_, int max_attempt_)
    : n(X_.n_rows), p(p_), q(X_.n_cols), K(K_), X(X_), kappa(kappa_), mu0(mu0_),
      prior_alpha(a_), prior_beta(b_), min_node(min_node_),
      max_attempt(max_attempt_), check_rows(false), n_prop(0), n_acc(0) {
    trees.resize(K);
    fit.zeros(p, n);
    r.set_size(p, n);
    const double sd = 1.0 / std::sqrt(kappa);
    for (int t = 0; t < K; ++t) {
      Node root = {-1, -1, -1, 0.0, false, 0};
      trees[t].nodes.push_back(root);
      trees[t].mu.set_size(p, 1);
      for (int j = 0; j < p; ++j) trees[t].mu(j, 0) = mu0[j] + sd * norm_rand();
      trees[t].leaf.assign(n, 0);
      fit.each_col() += trees[t].mu.col(0);
    }
    update_ranges();
  }

  void update_ranges() {
    xmin.set_size(q);
    xmax.set_size(q);
    for (int v = 0; v < q; ++v) {
      double lo = R_PosInf, hi = R_NegInf;
      const double* col = X.colptr(v);
      for (int i = 0; i < n; ++i) {
        if (std::isnan(col[i])) continue;
        if (col[i] < lo) lo = col[i];
        if (col[i] > hi) hi = col[i];
      }
      if (lo > hi) { lo = 0.0; hi = 0.0; }
      xmin[v] = lo;
      xmax[v] = hi;
    }
  }

  // One Metropolis-Hastings + Gibbs pass over all trees. `target` is p x n and
  // the residual precision is U diag(lam) U'.
  void sweep(const arma::mat& target, const arma::mat& U, const arma::vec& lam) {
    const arma::vec kUmu0 = kappa * (U.t() * mu0);
    const double leaf_const = 0.5 * (p * std::log(kappa) - kappa * arma::dot(mu0, mu0));
    arma::vec e(p);

    for (int t = 0; t < K; ++t) {
      Tree& tr = trees[t];

      // partial residuals
      for (int i = 0; i < n; ++i) {
        const double* m = tr.mu.colptr(tr.leaf[i]);
        double* f = fit.colptr(i);
        double* ri = r.colptr(i);
        const double* y = target.colptr(i);
        for (int j = 0; j < p; ++j) {
          f[j] -= m[j];
          ri[j] = y[j] - f[j];
        }
      }
      leaf_stats(tr.nodes, tr.leaf, S0, cnt0, ok0);

      // propose
      int move = 0;  // 0 grow, 1 prune, 2 change
      if (tr.nodes.size() > 1) {
        const double u = unif_rand() * 9.0;
        move = (u < 2.5) ? 0 : ((u < 5.0) ? 1 : 2);
      }
      bool valid = false;
      for (int attempt = 0; attempt < max_attempt && !valid; ++attempt) {
        propose(tr, move);
        leaf_stats(nodes1, leaf1, S1, cnt1, ok1);
        valid = is_valid(nodes1, cnt1, ok1);
      }
      n_prop += 1.0;
      if (valid) {
        const double ratio =
          (tree_prior(nodes1) + tree_loglik(nodes1, S1, cnt1, U, lam, kUmu0, leaf_const)) -
          (tree_prior(tr.nodes) + tree_loglik(tr.nodes, S0, cnt0, U, lam, kUmu0, leaf_const)) +
          log_hastings(tr.nodes, nodes1, move);
        if (ratio >= 0.0 || -exp_rand() <= ratio) {
          tr.nodes.swap(nodes1);
          tr.leaf.swap(leaf1);
          S0.swap(S1);
          cnt0.swap(cnt1);
          n_acc += 1.0;
        }
      }

      // terminal node parameters
      const int nn = tr.nodes.size();
      tr.mu.zeros(p, nn);
      for (int l = 0; l < nn; ++l) {
        if (tr.nodes[l].left >= 0) continue;
        const double* s = S0.colptr(l);
        for (int k = 0; k < p; ++k) {
          const double d = cnt0[l] * lam[k] + kappa;
          const double* u = U.colptr(k);
          double us = 0.0;
          for (int j = 0; j < p; ++j) us += u[j] * s[j];
          e[k] = (lam[k] * us + kUmu0[k]) / d + norm_rand() / std::sqrt(d);
        }
        double* m = tr.mu.colptr(l);
        for (int j = 0; j < p; ++j) {
          double acc = 0.0;
          for (int k = 0; k < p; ++k) acc += U(j, k) * e[k];
          m[j] = acc;
        }
      }
      for (int i = 0; i < n; ++i) {
        const double* m = tr.mu.colptr(tr.leaf[i]);
        double* f = fit.colptr(i);
        for (int j = 0; j < p; ++j) f[j] += m[j];
      }
    }
  }

  // Sum-of-trees prediction for new covariates (n_new x q); out is p x n_new.
  void predict(const arma::mat& Xnew, arma::mat& out) const {
    const int nnew = Xnew.n_rows;
    out.zeros(p, nnew);
    for (int t = 0; t < K; ++t) {
      const Tree& tr = trees[t];
      for (int i = 0; i < nnew; ++i) {
        const int l = descend(tr.nodes, Xnew.memptr() + i, nnew);
        const double* m = tr.mu.colptr(l);
        double* o = out.colptr(i);
        for (int j = 0; j < p; ++j) o[j] += m[j];
      }
    }
  }

  // Prediction and terminal nodes for a single covariate vector.
  void row_fit(const double* xrow, arma::vec& out, std::vector<int>& leaves) const {
    out.zeros(p);
    leaves.resize(K);
    for (int t = 0; t < K; ++t) {
      const int l = descend(trees[t].nodes, xrow, 1);
      leaves[t] = l;
      const double* m = trees[t].mu.colptr(l);
      for (int j = 0; j < p; ++j) out[j] += m[j];
    }
  }

  void set_row(int i, const std::vector<int>& leaves, const arma::vec& f) {
    for (int t = 0; t < K; ++t) trees[t].leaf[i] = leaves[t];
    fit.col(i) = f;
  }

  // All trees as one matrix with a row per node: tree number, row of the left
  // child (0 for a terminal node; the right child is the next row), splitting
  // variable, split value, direction of missing values, and the p terminal
  // node parameters. Indices are 1-based.
  NumericMatrix export_trees() const {
    int total = 0;
    for (int t = 0; t < K; ++t) total += trees[t].nodes.size();
    NumericMatrix out(total, 5 + p);
    int row = 0;
    for (int t = 0; t < K; ++t) {
      const Tree& tr = trees[t];
      const int base = row;
      for (std::size_t k = 0; k < tr.nodes.size(); ++k, ++row) {
        const Node& nd = tr.nodes[k];
        out(row, 0) = t + 1;
        out(row, 1) = nd.left >= 0 ? base + nd.left + 1 : 0;
        out(row, 2) = nd.left >= 0 ? nd.var + 1 : 0;
        out(row, 3) = nd.left >= 0 ? nd.val : 0.0;
        out(row, 4) = (nd.left >= 0 && nd.na_left) ? 1.0 : 0.0;
        if (nd.left < 0) for (int j = 0; j < p; ++j) out(row, 5 + j) = tr.mu(j, k);
      }
    }
    return out;
  }

  // Number of splitting rules using each variable.
  void var_counts(double* out) const {
    for (int v = 0; v < q; ++v) out[v] = 0.0;
    for (int t = 0; t < K; ++t)
      for (std::size_t k = 0; k < trees[t].nodes.size(); ++k)
        if (trees[t].nodes[k].left >= 0) out[trees[t].nodes[k].var] += 1.0;
  }

private:
  arma::mat r, S0, S1;
  std::vector<int> cnt0, cnt1, ok0, ok1, leaf1, pool;
  std::vector<Node> nodes1;

  void leaf_stats(const std::vector<Node>& nd, const std::vector<int>& leaf,
                  arma::mat& S, std::vector<int>& cnt, std::vector<int>& ok) const {
    const int nn = nd.size();
    S.zeros(p, nn);
    cnt.assign(nn, 0);
    if (check_rows) ok.assign(nn, 0);
    for (int i = 0; i < n; ++i) {
      const int l = leaf[i];
      ++cnt[l];
      if (check_rows && row_ok[i]) ++ok[l];
      double* s = S.colptr(l);
      const double* ri = r.colptr(i);
      for (int j = 0; j < p; ++j) s[j] += ri[j];
    }
  }

  bool is_valid(const std::vector<Node>& nd, const std::vector<int>& cnt,
                const std::vector<int>& ok) const {
    for (std::size_t l = 0; l < nd.size(); ++l) {
      if (nd[l].left >= 0) continue;
      if (cnt[l] < min_node) return false;
      if (check_rows && ok[l] == 0) return false;
    }
    return true;
  }

  double tree_prior(const std::vector<Node>& nd) const {
    double lp = 0.0;
    for (std::size_t k = 0; k < nd.size(); ++k) {
      const double ps = prior_alpha * std::pow(1.0 + nd[k].depth, -prior_beta);
      lp += std::log(nd[k].left >= 0 ? ps : 1.0 - ps);
    }
    return lp;
  }

  // Log marginal likelihood of the partial residuals given the tree, up to
  // terms that do not depend on the tree structure.
  double tree_loglik(const std::vector<Node>& nd, const arma::mat& S,
                     const std::vector<int>& cnt, const arma::mat& U,
                     const arma::vec& lam, const arma::vec& kUmu0,
                     double leaf_const) const {
    double ll = 0.0;
    for (std::size_t l = 0; l < nd.size(); ++l) {
      if (nd[l].left >= 0) continue;
      ll += leaf_loglik(cnt[l], S.colptr(l), p, kappa, U, lam, kUmu0) + leaf_const;
    }
    return ll;
  }

  // Probability of proposing each type of move from a given tree.
  static double move_prob(const std::vector<Node>& nd, int move) {
    if (nd.size() == 1) return move == 0 ? 1.0 : 0.0;
    return move == 2 ? 4.0 / 9.0 : 2.5 / 9.0;
  }

  static int n_leaves(const std::vector<Node>& nd) {
    int k = 0;
    for (std::size_t i = 0; i < nd.size(); ++i) if (nd[i].left < 0) ++k;
    return k;
  }

  // Internal nodes whose two children are both terminal.
  static int n_prunable(const std::vector<Node>& nd) {
    int k = 0;
    for (std::size_t i = 0; i < nd.size(); ++i) {
      const int c = nd[i].left;
      if (c >= 0 && nd[c].left < 0 && nd[c + 1].left < 0) ++k;
    }
    return k;
  }

  // Log of q(new -> old) / q(old -> new). A grow move picks one of the terminal
  // nodes and is reversed by pruning one of the prunable nodes of the new tree,
  // and vice versa. The splitting rule is proposed from its prior, so its
  // density cancels against the prior and does not appear, and a change move
  // contributes nothing.
  static double log_hastings(const std::vector<Node>& old_nd, const std::vector<Node>& new_nd, int move) {
    if (move == 0)
      return std::log(move_prob(new_nd, 1) / n_prunable(new_nd)) - std::log(move_prob(old_nd, 0) / n_leaves(old_nd));
    if (move == 1)
      return std::log(move_prob(new_nd, 0) / n_leaves(new_nd)) - std::log(move_prob(old_nd, 1) / n_prunable(old_nd));
    return 0.0;
  }

  // Interval of variable v still available at `node`, given the splitting
  // rules of its ancestors.
  void bounds(const std::vector<Node>& nd, int node, int v, double& lo, double& hi) const {
    lo = xmin[v];
    hi = xmax[v];
    int cur = node;
    while (nd[cur].parent >= 0) {
      const int par = nd[cur].parent;
      if (nd[par].var == v) {
        if (nd[par].left == cur) hi = std::min(hi, nd[par].val);
        else lo = std::max(lo, nd[par].val);
      }
      cur = par;
    }
  }

  inline int pick(int size) const {
    int k = static_cast<int>(unif_rand() * size);
    return k >= size ? size - 1 : k;
  }

  // Fills nodes1 / leaf1 with a proposal.
  void propose(const Tree& tr, int move) {
    nodes1 = tr.nodes;
    leaf1 = tr.leaf;
    const int nn = nodes1.size();
    pool.clear();

    if (move == 0) {  // grow a terminal node
      for (int k = 0; k < nn; ++k) if (nodes1[k].left < 0) pool.push_back(k);
      const int g = pool[pick(pool.size())];
      const int v = pick(q);
      double lo, hi;
      bounds(nodes1, g, v, lo, hi);
      const double val = (hi > lo) ? R::runif(lo, hi) : lo;
      const bool na_left = unif_rand() < 0.5;
      const int d = nodes1[g].depth + 1;
      Node child = {g, -1, -1, 0.0, false, d};
      nodes1.push_back(child);
      nodes1.push_back(child);
      nodes1[g].left = nn;
      nodes1[g].var = v;
      nodes1[g].val = val;
      nodes1[g].na_left = na_left;
      const double* col = X.colptr(v);
      for (int i = 0; i < n; ++i) {
        if (leaf1[i] != g) continue;
        const bool go_left = std::isnan(col[i]) ? na_left : (col[i] <= val);
        leaf1[i] = go_left ? nn : nn + 1;
      }
    } else if (move == 1) {  // prune a pair of terminal nodes
      for (int k = 0; k < nn; ++k) {
        const int c = nodes1[k].left;
        if (c >= 0 && nodes1[c].left < 0 && nodes1[c + 1].left < 0) pool.push_back(k);
      }
      const int g = pool[pick(pool.size())];
      const int c = nodes1[g].left;
      nodes1.erase(nodes1.begin() + c, nodes1.begin() + c + 2);
      nodes1[g].left = -1;
      nodes1[g].var = -1;
      for (std::size_t k = 0; k < nodes1.size(); ++k) {
        if (nodes1[k].parent > c) nodes1[k].parent -= 2;
        if (nodes1[k].left > c) nodes1[k].left -= 2;
      }
      for (int i = 0; i < n; ++i) {
        if (leaf1[i] == c || leaf1[i] == c + 1) leaf1[i] = g;
        else if (leaf1[i] > c) leaf1[i] -= 2;
      }
    } else {  // change the rule of an internal node
      for (int k = 0; k < nn; ++k) if (nodes1[k].left >= 0) pool.push_back(k);
      const int g = pool[pick(pool.size())];
      const int v = pick(q);
      double lo, hi;
      bounds(nodes1, g, v, lo, hi);
      nodes1[g].var = v;
      nodes1[g].val = (hi > lo) ? R::runif(lo, hi) : lo;
      nodes1[g].na_left = unif_rand() < 0.5;
      for (int i = 0; i < n; ++i) leaf1[i] = descend(nodes1, X.memptr() + i, n);
    }
  }
};

// ----------------------------------------------------------------------------
// Gaussian update of missing responses
// ----------------------------------------------------------------------------

struct Pattern {
  arma::uvec mis, obs;
  std::vector<int> rows;
};

// Each row's density is proportional to exp(-y'Qy/2 + b'y) with b the
// corresponding column of Bt. conditional = true draws the missing entries
// given the observed entries of the same row; conditional = false uses the
// corresponding margin of N(Q^{-1}b, Q^{-1}) instead (Appendix C.1 as printed).
static void draw_ymis(arma::mat& Yt, const arma::mat& Q, const arma::mat& Bt,
                      const std::vector<Pattern>& patterns, bool conditional) {
  arma::mat Sfull;
  if (!conditional) Sfull = safe_inv(Q);
  for (std::size_t g = 0; g < patterns.size(); ++g) {
    const Pattern& pt = patterns[g];
    const int pm = pt.mis.n_elem;
    arma::vec eps(pm), draw;
    if (conditional) {
      arma::mat Uc;
      if (!safe_chol(Uc, Q.submat(pt.mis, pt.mis), "upper")) stop("precision of missing responses is not positive definite");
      arma::mat Qmo;
      if (pt.obs.n_elem > 0) Qmo = Q.submat(pt.mis, pt.obs);
      for (std::size_t k = 0; k < pt.rows.size(); ++k) {
        const int i = pt.rows[k];
        arma::vec yi = Yt.col(i), bi = Bt.col(i);
        arma::vec rhs = bi.elem(pt.mis);
        if (pt.obs.n_elem > 0) rhs -= Qmo * yi.elem(pt.obs);
        for (int j = 0; j < pm; ++j) eps[j] = norm_rand();
        draw = arma::solve(arma::trimatu(Uc), arma::solve(arma::trimatl(Uc.t()), rhs) + eps);
        for (int j = 0; j < pm; ++j) Yt(pt.mis[j], i) = draw[j];
      }
    } else {
      arma::mat Lc;
      if (!safe_chol(Lc, Sfull.submat(pt.mis, pt.mis), "lower")) stop("covariance of missing responses is not positive definite");
      for (std::size_t k = 0; k < pt.rows.size(); ++k) {
        const int i = pt.rows[k];
        arma::vec mean = Sfull * Bt.col(i);
        for (int j = 0; j < pm; ++j) eps[j] = norm_rand();
        draw = mean.elem(pt.mis) + Lc * eps;
        for (int j = 0; j < pm; ++j) Yt(pt.mis[j], i) = draw[j];
      }
    }
  }
}

static inline double quad_form(const double* a, const double* b, const arma::mat& Q, int p) {
  double out = 0.0;
  for (int j = 0; j < p; ++j) {
    const double dj = a[j] - b[j];
    const double* qc = Q.colptr(j);
    double acc = 0.0;
    for (int k = 0; k < p; ++k) acc += qc[k] * (a[k] - b[k]);
    out += dj * acc;
  }
  return out;
}

// Transposes a p x n matrix on the BART scale and maps it back to the scale of
// the data.
static arma::mat to_original(const arma::mat& At, const arma::vec& y_min, const arma::vec& y_rng) {
  arma::mat out = At.t();
  for (arma::uword j = 0; j < out.n_cols; ++j) out.col(j) = (out.col(j) + 0.5) * y_rng[j] + y_min[j];
  return out;
}

// ----------------------------------------------------------------------------
// The sampler
// ----------------------------------------------------------------------------

// [[Rcpp::export(name = ".missbart_cpp")]]
List missbart_cpp(const arma::mat& x, const arma::mat& y_init, const IntegerMatrix& m,
                  const arma::mat& x_new, int model, List ctrl) {
  const int n = x.n_rows, q = x.n_cols, p = y_init.n_cols;
  const int n_new = x_new.n_rows;
  const bool do_predict = n_new > 0;

  const int n_reg_trees = ctrl["n_reg_trees"], n_class_trees = ctrl["n_class_trees"];
  const int burn = ctrl["burn"], iters = ctrl["iters"], thin = ctrl["thin"];
  const double kappa = ctrl["kappa"], kappa_class = ctrl["kappa_class"];
  const arma::vec mu0 = as<arma::vec>(ctrl["mu0"]);
  const double wish_df = ctrl["wish_df"];
  const arma::mat Vinv = as<arma::mat>(ctrl["Vinv"]);
  const arma::vec omega_init = as<arma::vec>(ctrl["omega_init"]);
  const double prior_alpha = ctrl["prior_alpha"], prior_beta = ctrl["prior_beta"];
  const int min_node = ctrl["min_node"], max_attempt = ctrl["max_attempt"];
  const bool include_x = ctrl["include_x"], include_y = ctrl["include_y"];
  const bool conditional = ctrl["conditional"];
  const arma::vec z_offset = as<arma::vec>(ctrl["z_offset"]);
  const double mh_sd = ctrl["MH_sd"];
  const arma::vec y_min = as<arma::vec>(ctrl["y_min"]), y_rng = as<arma::vec>(ctrl["y_rng"]);
  const double a0 = ctrl["a0"], b0 = ctrl["b0"], ax = ctrl["ax"], bx = ctrl["bx"],
               ay = ctrl["ay"], by = ctrl["by"];
  const bool save_y_pred = ctrl["save_y_pred"], save_z = ctrl["save_z"];
  const bool make_pdp = ctrl["make_pdp"];
  const bool keep_trees = ctrl["keep_trees"];
  const arma::vec pdp_grid = as<arma::vec>(ctrl["pdp_grid"]);
  const bool show_progress = ctrl["show_progress"];
  const int progress_every = std::max(1, as<int>(ctrl["progress_every"]));
  const int total = burn + thin * iters;

  // ----- data, stored one observation per column -----
  arma::mat Yt = y_init.t();
  std::vector<int> mis_i, mis_j;                 // missing cells, column-major
  std::vector<std::vector<int> > row_mis(n);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < n; ++i)
      if (m(i, j) == 0) { mis_i.push_back(i); mis_j.push_back(j); row_mis[i].push_back(j); }
  const int n_mis = mis_i.size();
  std::vector<int> rows_mis;
  std::vector<Pattern> patterns;
  {
    std::map<std::vector<int>, int> seen;
    for (int i = 0; i < n; ++i) {
      if (row_mis[i].empty()) continue;
      rows_mis.push_back(i);
      std::map<std::vector<int>, int>::iterator it = seen.find(row_mis[i]);
      if (it == seen.end()) {
        Pattern pt;
        std::vector<char> is_mis(p, 0);
        for (std::size_t k = 0; k < row_mis[i].size(); ++k) is_mis[row_mis[i][k]] = 1;
        pt.mis.set_size(row_mis[i].size());
        pt.obs.set_size(p - row_mis[i].size());
        int a = 0, b = 0;
        for (int j = 0; j < p; ++j) { if (is_mis[j]) pt.mis[a++] = j; else pt.obs[b++] = j; }
        patterns.push_back(pt);
        seen[row_mis[i]] = patterns.size() - 1;
        it = seen.find(row_mis[i]);
      }
      patterns[it->second].rows.push_back(i);
    }
  }

  // ----- regression forest and residual precision -----
  Forest reg(x, p, n_reg_trees, kappa, mu0, prior_alpha, prior_beta, min_node, max_attempt);
  arma::mat Omega = arma::diagmat(omega_init);
  arma::mat U;
  arma::vec lam;
  arma::eig_sym(lam, U, Omega);

  // ----- missingness model -----
  const int qx = include_x ? q : 0, py = include_y ? p : 0;
  const int r1 = 1 + qx + py;       // probit regression predictors
  const int r2 = qx + py;           // probit BART predictors
  arma::mat Zd, Zc;                 // design matrices of the two models
  arma::mat B, Rm, Zl;              // model 1: coefficients, correlation, latent (n x p)
  arma::vec tau;
  arma::mat Zt, Zc_target;          // model 2: latent, p x n, and latent minus offset
  const arma::mat Ip = arma::eye(p, p);
  const arma::vec ones_p = arma::ones(p);
  Forest* cls = 0;
  int ycol0 = 0;

  if (model == 1) {
    Zd.ones(n, r1);
    if (include_x) Zd.cols(1, q) = x;
    ycol0 = 1 + qx;
    if (include_y) Zd.cols(ycol0, ycol0 + p - 1) = y_init;
    B.zeros(r1, p);
    Rm = Ip;
    tau.set_size(r1);
    tau[0] = R::rgamma(0.5 * p + a0, 1.0 / b0);
    if (qx > 0) tau.subvec(1, qx).fill(R::rgamma(0.5 * p * qx + ax, 1.0 / bx));
    if (py > 0) tau.subvec(1 + qx, r1 - 1).fill(R::rgamma(0.5 * p * py + ay, 1.0 / by));
    Zl.set_size(n, p);
    for (int j = 0; j < p; ++j) for (int i = 0; i < n; ++i) Zl(i, j) = m(i, j) == 0 ? 1.0 : -1.0;
  } else if (model == 2) {
    Zc.set_size(n, r2);
    if (include_x) Zc.cols(0, q - 1) = x;
    ycol0 = qx;
    if (include_y) Zc.cols(ycol0, ycol0 + p - 1) = y_init;
    cls = new Forest(Zc, p, n_class_trees, kappa_class, arma::zeros(p), prior_alpha,
                     prior_beta, min_node, max_attempt);
    cls->check_rows = true;
    cls->row_ok.assign(n, 1);
    for (int i = 0; i < n; ++i) if (!row_mis[i].empty()) cls->row_ok[i] = 0;
    Zt.set_size(p, n);
    for (int j = 0; j < p; ++j) for (int i = 0; i < n; ++i) Zt(j, i) = m(i, j) == 0 ? 1.0 : -1.0;
  }

  // ----- storage -----
  List y_post(iters), y_pred(save_y_pred ? iters : 0), sigma_post(iters);
  List new_y_hat(do_predict ? iters : 0), new_y_post(do_predict ? iters : 0);
  NumericMatrix y_impute(iters, n_mis);
  NumericMatrix var_imp(iters, q);
  List B_post(model == 1 ? iters : 0), R_post(model == 1 ? iters : 0);
  NumericMatrix tau_post(model == 1 ? iters : 0, model == 1 ? r1 : 0);
  List z_post((save_z && model > 0) ? iters : 0);
  NumericMatrix var_imp_miss(model == 2 ? iters : 0, model == 2 ? r2 : 0);
  const bool pdp_on = make_pdp && model == 2 && p == 1 && include_y;
  List pdp_out(pdp_on ? iters : 0);
  List reg_trees(keep_trees ? iters : 0), class_trees((keep_trees && model == 2) ? iters : 0);
  NumericVector mh_accept(n_mis);
  std::vector<double> row_acc(n, 0.0);

  arma::mat Bt(p, n), Res(p, n), new_fit, tmp;
  arma::vec vbuf(q), vbuf2(std::max(r2, 1)), xrow(std::max(r2, 1)), zprop(p), yprop(p);
  std::vector<int> leaves;
  int s = 0;

  try {
  for (int it = 1; it <= total; ++it) {
    if (it % 16 == 0) checkUserInterrupt();
    if (show_progress && (it == 1 || it % progress_every == 0 || it == total)) {
      Rprintf("\r  iteration %d / %d", it, total);
      if (it == total) Rprintf("\n");
    }

    // 1. regression trees and their terminal node parameters
    reg.sweep(Yt, U, lam);

    if (model == 1) {
      // 2. missing responses
      if (n_mis > 0) {
        arma::mat Q = Omega;
        Bt = Omega * reg.fit;
        if (include_y) {
          const arma::mat P = safe_inv(Rm);
          const arma::mat By = B.rows(ycol0, r1 - 1);
          const arma::mat BP = By * P;
          Q += BP * By.t();
          Bt += BP * (Zl - Zd.cols(0, ycol0 - 1) * B.rows(0, ycol0 - 1)).t();
        }
        draw_ymis(Yt, Q, Bt, patterns, conditional);
        if (include_y)
          for (int k = 0; k < n_mis; ++k) Zd(mis_i[k], ycol0 + mis_j[k]) = Yt(mis_j[k], mis_i[k]);
      }
    }

    // 3. residual precision
    Res = Yt - reg.fit;
    Omega = rwish(wish_df + n, safe_inv(Res * Res.t() + Vinv));
    arma::eig_sym(lam, U, Omega);

    if (model == 1) {
      // 4. latent probit variables
      const arma::mat Mu = Zd * B;
      if (p == 1) {
        for (int i = 0; i < n; ++i) Zl(i, 0) = rtn_sign(Mu(i, 0), 1.0, m(i, 0) == 0);
      } else {
        const arma::mat P = safe_inv(Rm);
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < p; ++j) {
            double acc = 0.0;
            for (int k = 0; k < p; ++k) if (k != j) acc += P(j, k) * (Zl(i, k) - Mu(i, k));
            const double cm = Mu(i, j) - acc / P(j, j);
            Zl(i, j) = rtn_sign(cm, 1.0 / std::sqrt(P(j, j)), m(i, j) == 0);
          }
        }
      }

      // 5. probit coefficients and correlation matrix
      arma::mat Xi_inv = Zd.t() * Zd;
      Xi_inv.diag() += tau;
      if (p == 1) {
        arma::mat Uc;
        if (!safe_chol(Uc, Xi_inv, "upper")) stop("probit precision is not positive definite");
        arma::vec eps(r1);
        for (int k = 0; k < r1; ++k) eps[k] = norm_rand();
        B = arma::solve(arma::trimatu(Uc), arma::solve(arma::trimatl(Uc.t()), Zd.t() * Zl) + eps);
      } else {
        // parameter expansion of Talhouk, Doucet & Murphy (2012)
        const arma::mat P = safe_inv(Rm);
        arma::rowvec d(p);
        for (int j = 0; j < p; ++j) d[j] = std::sqrt(1.0 / R::rgamma(0.5 * (p + 1), 2.0 / P(j, j)));
        arma::mat W = Zl;
        W.each_row() %= d;
        const arma::mat Xi = safe_inv(Xi_inv);
        const arma::mat M = Xi * (Zd.t() * W);
        arma::mat Sc = W.t() * W + Ip - M.t() * Xi_inv * M;
        const arma::mat Sigma = riwish(n + p + 1.0, Sc);
        arma::mat Lxi, Usig;
        if (!safe_chol(Lxi, Xi, "lower") || !safe_chol(Usig, Sigma, "upper"))
          stop("probit covariance is not positive definite");
        arma::mat E(r1, p);
        for (int j = 0; j < p; ++j) for (int k = 0; k < r1; ++k) E(k, j) = norm_rand();
        const arma::mat gamma = M + Lxi * E * Usig;
        const arma::rowvec Dv = arma::sqrt(Sigma.diag().t());
        Rm = Sigma / (Dv.t() * Dv);
        Rm.diag().ones();
        B = gamma;
        B.each_row() /= Dv;
        // The latent variables move with the expansion parameters: they are kept
        // between iterations (only one Gibbs sweep is made above), so they must
        // be mapped back to the correlation scale together with B and R.
        Zl = W;
        Zl.each_row() /= Dv;
      }

      // 6. prior precisions of the probit coefficients
      {
        const arma::mat P = safe_inv(Rm);
        const arma::vec Ad = arma::sum((B * P) % B, 1);
        tau[0] = R::rgamma(0.5 * p + a0, 1.0 / (0.5 * Ad[0] + b0));
        if (qx > 0)
          tau.subvec(1, qx).fill(R::rgamma(0.5 * p * qx + ax, 1.0 / (0.5 * arma::sum(Ad.subvec(1, qx)) + bx)));
        if (py > 0)
          tau.subvec(1 + qx, r1 - 1).fill(R::rgamma(0.5 * p * py + ay, 1.0 / (0.5 * arma::sum(Ad.subvec(1 + qx, r1 - 1)) + by)));
      }
    } else if (model == 2) {
      // 4. missingness trees, then the latent probit variables
      cls->update_ranges();
      Zc_target = Zt;
      Zc_target.each_col() -= z_offset;
      cls->sweep(Zc_target, Ip, ones_p);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j) Zt(j, i) = rtn_sign(z_offset[j] + cls->fit(j, i), 1.0, m(i, j) == 0);

      // 5. missing responses
      if (n_mis > 0) {
        if (include_y) {
          for (std::size_t k = 0; k < rows_mis.size(); ++k) {
            const int i = rows_mis[k];
            for (int j = 0; j < p; ++j) yprop[j] = Yt(j, i);
            for (std::size_t c = 0; c < row_mis[i].size(); ++c) yprop[row_mis[i][c]] += mh_sd * norm_rand();
            for (int v = 0; v < r2; ++v) xrow[v] = Zc(i, v);
            for (int j = 0; j < p; ++j) xrow[ycol0 + j] = yprop[j];
            cls->row_fit(xrow.memptr(), zprop, leaves);
            double l1 = 0.0, l2 = 0.0;
            for (int j = 0; j < p; ++j) {
              const double d1 = Zt(j, i) - z_offset[j] - cls->fit(j, i), d2 = Zt(j, i) - z_offset[j] - zprop[j];
              l1 += d1 * d1;
              l2 += d2 * d2;
            }
            const double ratio =
              -0.5 * (quad_form(yprop.memptr(), reg.fit.colptr(i), Omega, p) + l2) +
               0.5 * (quad_form(Yt.colptr(i), reg.fit.colptr(i), Omega, p) + l1);
            if (ratio >= 0.0 || -exp_rand() <= ratio) {
              for (int j = 0; j < p; ++j) { Yt(j, i) = yprop[j]; Zc(i, ycol0 + j) = yprop[j]; }
              cls->set_row(i, leaves, zprop);
              row_acc[i] += 1.0;
            }
          }
        } else {
          Bt = Omega * reg.fit;
          draw_ymis(Yt, Omega, Bt, patterns, conditional);
        }
      }
    }

    // ----- store -----
    if (it > burn && (it - burn) % thin == 0) {
      const arma::mat Sigma = safe_inv(Omega);
      arma::mat Lsig;
      safe_chol(Lsig, Sigma, "lower");

      tmp = to_original(reg.fit, y_min, y_rng);
      y_post[s] = tmp;
      if (save_y_pred) {
        arma::mat E(p, n);
        for (int i = 0; i < n; ++i) for (int j = 0; j < p; ++j) E(j, i) = norm_rand();
        tmp = to_original(reg.fit + Lsig * E, y_min, y_rng);
        y_pred[s] = tmp;
      }
      sigma_post[s] = Sigma % (y_rng * y_rng.t());
      for (int k = 0; k < n_mis; ++k)
        y_impute(s, k) = (Yt(mis_j[k], mis_i[k]) + 0.5) * y_rng[mis_j[k]] + y_min[mis_j[k]];
      reg.var_counts(vbuf.memptr());
      for (int v = 0; v < q; ++v) var_imp(s, v) = vbuf[v];
      if (keep_trees) {
        reg_trees[s] = reg.export_trees();
        if (model == 2) class_trees[s] = cls->export_trees();
      }

      if (do_predict) {
        reg.predict(x_new, new_fit);
        tmp = to_original(new_fit, y_min, y_rng);
        new_y_hat[s] = tmp;
        arma::mat E(p, n_new);
        for (int i = 0; i < n_new; ++i) for (int j = 0; j < p; ++j) E(j, i) = norm_rand();
        tmp = to_original(new_fit + Lsig * E, y_min, y_rng);
        new_y_post[s] = tmp;
      }

      if (model == 1) {
        B_post[s] = B;
        R_post[s] = Rm;
        for (int k = 0; k < r1; ++k) tau_post(s, k) = tau[k];
        if (save_z) z_post[s] = Zl;
      } else if (model == 2) {
        cls->var_counts(vbuf2.memptr());
        for (int v = 0; v < r2; ++v) var_imp_miss(s, v) = vbuf2[v];
        if (save_z) { tmp = Zt.t(); z_post[s] = tmp; }
        if (pdp_on) {
          // latent missingness predictor over a grid of response values
          const int G = pdp_grid.n_elem, nc = include_x ? n : 1;
          arma::mat pd(nc, G);
          for (int i = 0; i < nc; ++i) {
            for (int v = 0; v < r2; ++v) xrow[v] = Zc(i, v);
            for (int g = 0; g < G; ++g) {
              xrow[ycol0] = pdp_grid[g];
              cls->row_fit(xrow.memptr(), zprop, leaves);
              pd(i, g) = z_offset[0] + zprop[0];
            }
          }
          pdp_out[s] = pd;
        }
      }
      ++s;
    }
  }
  } catch (...) {
    delete cls;
    throw;
  }

  for (int k = 0; k < n_mis; ++k) mh_accept[k] = row_acc[mis_i[k]] / total;
  const double reg_acc = reg.n_acc / std::max(1.0, reg.n_prop);
  const double cls_acc = cls ? cls->n_acc / std::max(1.0, cls->n_prop) : NA_REAL;
  delete cls;

  return List::create(
    _["y_post"] = y_post, _["y_pred"] = y_pred, _["sigma_post"] = sigma_post,
    _["y_impute"] = y_impute, _["new_y_hat"] = new_y_hat, _["new_y_post"] = new_y_post,
    _["var_imp"] = var_imp, _["B_post"] = B_post, _["R_post"] = R_post,
    _["tau_b"] = tau_post, _["z_post"] = z_post, _["var_imp_miss"] = var_imp_miss,
    _["pdp_out"] = pdp_out, _["y_miss_accept"] = mh_accept,
    _["reg_trees"] = reg_trees, _["class_trees"] = class_trees,
    _["tree_accept"] = NumericVector::create(_["regression"] = reg_acc, _["missingness"] = cls_acc));
}

// ----------------------------------------------------------------------------
// Hooks used by the unit tests
// ----------------------------------------------------------------------------

// [[Rcpp::export(name = ".rtnorm_cpp")]]
NumericVector rtnorm_cpp(NumericVector mu, double sd, bool positive) {
  NumericVector out(mu.size());
  for (R_xlen_t i = 0; i < mu.size(); ++i) out[i] = rtn_sign(mu[i], sd, positive);
  return out;
}

// [[Rcpp::export(name = ".rwish_cpp")]]
arma::mat rwish_cpp(double df, const arma::mat& S, bool inverse) {
  return inverse ? riwish(df, S) : rwish(df, S);
}

// Terminal node of every row of x for a tree given as parallel vectors
// (0-based indices, left = -1 for terminal nodes).
// [[Rcpp::export(name = ".descend_cpp")]]
IntegerVector descend_cpp(const arma::mat& x, IntegerVector left, IntegerVector var,
                          NumericVector val, LogicalVector na_left) {
  std::vector<Node> nd(left.size());
  for (R_xlen_t k = 0; k < left.size(); ++k) {
    Node nk = {-1, left[k], var[k], val[k], na_left[k] != 0, 0};
    nd[k] = nk;
  }
  IntegerVector out(x.n_rows);
  for (arma::uword i = 0; i < x.n_rows; ++i) out[i] = descend(nd, x.memptr() + i, x.n_rows);
  return out;
}

// [[Rcpp::export(name = ".leaf_loglik_cpp")]]
double leaf_loglik_cpp(int cnt, const arma::vec& s, const arma::mat& Omega, double kappa,
                       const arma::vec& mu0) {
  arma::mat U;
  arma::vec lam;
  arma::eig_sym(lam, U, Omega);
  const int p = s.n_elem;
  return leaf_loglik(cnt, s.memptr(), p, kappa, U, lam, kappa * (U.t() * mu0)) +
    0.5 * (p * std::log(kappa) - kappa * arma::dot(mu0, mu0));
}

// One Gaussian update of the missing entries of y (n x p; m is 0 for missing),
// where row i has density proportional to exp(-y'Qy/2 + b_i'y).
// [[Rcpp::export(name = ".draw_ymis_cpp")]]
arma::mat draw_ymis_cpp(const arma::mat& y, const IntegerMatrix& m, const arma::mat& Q,
                        const arma::mat& b, bool conditional) {
  const int n = y.n_rows, p = y.n_cols;
  std::vector<Pattern> patterns;
  for (int i = 0; i < n; ++i) {
    Pattern pt;
    std::vector<arma::uword> mis, obs;
    for (int j = 0; j < p; ++j) { if (m(i, j) == 0) mis.push_back(j); else obs.push_back(j); }
    if (mis.empty()) continue;
    pt.mis = arma::uvec(mis);
    pt.obs = arma::uvec(obs);
    pt.rows.push_back(i);
    patterns.push_back(pt);
  }
  arma::mat Yt = y.t();
  draw_ymis(Yt, Q, b.t(), patterns, conditional);
  return Yt.t();
}

// ----------------------------------------------------------------------------
// Predictions from stored trees
// ----------------------------------------------------------------------------

// Adds the sum-of-trees prediction for one observation to out[0..p-1]. `tr` is
// a matrix written by Forest::export_trees(), `roots` the rows of the roots.
static inline void stored_predict(const NumericMatrix& tr, const std::vector<int>& roots,
                                  const double* x, std::size_t stride, int p, double* out) {
  const double* left = &tr(0, 1);
  const double* var = &tr(0, 2);
  const double* val = &tr(0, 3);
  const double* nal = &tr(0, 4);
  for (std::size_t t = 0; t < roots.size(); ++t) {
    int k = roots[t];
    while (left[k] > 0) {
      const double v = x[stride * (static_cast<int>(var[k]) - 1)];
      const bool go_left = std::isnan(v) ? (nal[k] > 0) : (v <= val[k]);
      k = static_cast<int>(left[k]) - 1 + (go_left ? 0 : 1);
    }
    for (int j = 0; j < p; ++j) out[j] += tr(k, 5 + j);
  }
}

static std::vector<int> root_rows(const NumericMatrix& tr) {
  std::vector<int> roots;
  double last = 0.0;
  for (int k = 0; k < tr.nrow(); ++k) if (tr(k, 0) != last) { roots.push_back(k); last = tr(k, 0); }
  return roots;
}

// Sum-of-trees predictions (n x p x draws) on the scale used by the sampler.
// [[Rcpp::export(name = ".forest_predict_cpp")]]
arma::cube forest_predict_cpp(const List& trees, const arma::mat& X, int p) {
  const int n = X.n_rows, D = trees.size();
  arma::cube out(n, p, D, arma::fill::zeros);
  std::vector<double> buf(p);
  for (int d = 0; d < D; ++d) {
    const NumericMatrix tr = trees[d];
    const std::vector<int> roots = root_rows(tr);
    for (int i = 0; i < n; ++i) {
      std::fill(buf.begin(), buf.end(), 0.0);
      stored_predict(tr, roots, X.memptr() + i, n, p, buf.data());
      for (int j = 0; j < p; ++j) out(i, j, d) = buf[j];
    }
    if (d % 16 == 0) checkUserInterrupt();
  }
  return out;
}

// Partial dependence: predictions with column `col` (1-based) of X set to each
// value of `grid` in turn. If the same draw-specific values should be used for
// some columns (imputed responses), they are taken from X as given. With
// probit = true the predictions are pnorm(offset + sum of trees). Returns the
// individual curves averaged over draws (n x G x p) and the curve averaged
// over observations for every draw (D x G x p).
// [[Rcpp::export(name = ".forest_pdp_cpp")]]
List forest_pdp_cpp(const List& trees, const arma::mat& X, int col, const arma::vec& grid, int p,
                    const arma::vec& offset, bool probit) {
  const int n = X.n_rows, q = X.n_cols, D = trees.size(), G = grid.n_elem;
  arma::cube ice(n, G, p, arma::fill::zeros), pdp(D, G, p, arma::fill::zeros);
  std::vector<double> row(q), buf(p);
  for (int d = 0; d < D; ++d) {
    const NumericMatrix tr = trees[d];
    const std::vector<int> roots = root_rows(tr);
    for (int i = 0; i < n; ++i) {
      for (int v = 0; v < q; ++v) row[v] = X(i, v);
      for (int g = 0; g < G; ++g) {
        row[col - 1] = grid[g];
        std::fill(buf.begin(), buf.end(), 0.0);
        stored_predict(tr, roots, row.data(), 1, p, buf.data());
        for (int j = 0; j < p; ++j) {
          const double v = probit ? R::pnorm(offset[j] + buf[j], 0.0, 1.0, 1, 0) : buf[j];
          ice(i, g, j) += v / D;
          pdp(d, g, j) += v / n;
        }
      }
    }
    if (d % 4 == 0) checkUserInterrupt();
  }
  return List::create(_["ice"] = ice, _["pdp"] = pdp);
}
