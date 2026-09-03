#pragma once
/**
 * @file
 * @brief The correlation between independent fission yields, imported from a published product.
 * @ingroup nucdata
 */
//
// THE ONE PIECE OF THIS TOOL'S UNCERTAINTY STORY THAT CANNOT BE DERIVED OR STAGED.
//
// ENDF carries no covariance for decay data at all -- the tapes hold MF1 and MF8 and nothing
// else -- and the branching correlation NuSIFT reports is derived rather than imported, because
// "a nuclide's modes sum to one" is a fact about the data model. Fission yields are different in
// both respects. They are correlated through the mass and charge balance the evaluation was
// fitted under, that correlation is not recoverable from the yields themselves, and no
// evaluation publishes it. It exists only as a separate product:
//
//   E.F. Matthews, L.A. Bernstein, W. Younes, "Stochastically estimated covariance matrices for
//   independent and cumulative fission yields in the ENDF/B-VIII.0 and JEFF-3.3 evaluations",
//   Atomic Data and Nuclear Data Tables 140 (2021) 101441, doi:10.1016/j.adt.2021.101441.
//   Matrices published as FYCoM, doi:10.5281/zenodo.5985936.
//
// The method resamples the evaluated yields under mass, charge and normalization conservation,
// so what it estimates is the correlation the EVALUATION PROCESS induces. Four decisions follow
// from what that product is, and each is enforced here rather than left to the caller.
//
// ONLY THE CORRELATION IS IMPORTED, never the covariance. The method's own primary variances are
// wider than the evaluation's -- the paper says so, and supplies a "normalized" covariance,
// correlation times evaluated variance, as its own remedy. Pairing the imported correlation with
// the sigma_Y the store stages IS that construction, one edition later, and it buys an identity
// worth having: the diagonal of the assembled Sigma is exactly the diagonal figure already
// reported, so the correlated and diagonal answers differ by the off-diagonal and by nothing
// else. A reader can therefore attribute the whole difference to the correlation.
//
// THE EDITION PAIRING IS DECLARED, NOT ASSUMED. The matrices are built for ENDF/B-VIII.0 and
// JEFF-3.3; the store is VIII.1. A correlation is dimensionless and structural, which is the
// argument that it transports across an edition better than a variance would -- but it is an
// argument, not a proof, so `library` travels with the matrix into the report, beside the
// citation the licence requires.
//
// THE MATRIX IS NOT POSITIVE SEMI-DEFINITE, and this is measured rather than feared. The
// published U-235 thermal correlation has 280 negative eigenvalues against a trace of 998, the
// most negative -2.95, and rank 308 -- the signature of a stochastic estimate with fewer samples
// than dimensions. Restriction to the products a seed actually carries does not repair it: a
// principal submatrix inherits positive semi-definiteness, but there is none to inherit. So a
// quadratic form g^T Sigma g can in principle come out NEGATIVE, and the propagation checks the
// sign and reports the diagnostic instead of trusting it.
//
// AND IT IS NOT PROJECTED ONTO THE NEAREST PSD MATRIX. That repair is cheap and standard and is
// deliberately not done: the value of this figure is that it is the published product rather
// than NuSIFT's smoothing of it, and a projection would quietly replace an evaluated correlation
// with one no evaluation stands behind.
//
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace nusift {

// Where a matrix came from, carried into every answer built on it. All four fields reach the
// report: a number from this file is not interpretable without knowing which fissioning system
// and which evaluation produced it, and the licence releases the data for use on the condition
// that its author is cited.
struct YieldCovarianceProvenance {
  std::string path;      // the file as the user named it
  std::string system;    // the system code the filename carries, e.g. "U235T"
  std::string library;   // the evaluation the matrix was BUILT FOR, e.g. "ENDF/B-VIII.0"
  std::string citation;  // the paper, spelled out
};

// A symmetric correlation matrix over fission products, keyed by NuSIFT's own ZAI key.
//
// Dense and in memory: the published systems run to about a thousand products, so the largest is
// some 8 MB of double, and every use contracts it against a vector that touches most rows. A
// sparse form would save nothing and cost the symmetry check.
class YieldCorrelation {
public:
  // Read a FYCoM correlation CSV. The format is its own: a header row of product keys preceded
  // by an empty cell, then one row per product, each beginning with its key.
  //
  // KEYS ARE TRANSLATED ON READ, because the two conventions collide silently. FYCoM writes
  // Z*10000 + M*1000 + A and NuSIFT writes Z*10000 + A*10 + I, and 290068 is a legal key in
  // both -- Cu-68 in one, and nothing sensible in the other. Translating at the boundary means
  // no other code has to know FYCoM's convention exists.
  //
  // Throws InputError for a file that cannot be opened, a ragged or non-square matrix, a
  // duplicate key, a diagonal that is not one, or an asymmetry beyond a rounding tolerance --
  // each of which means the file is not what it claims to be rather than that the data is
  // inconvenient.
  static YieldCorrelation read(const std::string& path, const std::string& library);

  std::size_t size() const { return keys_.size(); }
  const std::vector<std::int64_t>& keys() const { return keys_; }
  const YieldCovarianceProvenance& provenance() const { return provenance_; }

  // Row-major element access. Both indices must be in range.
  double at(std::size_t i, std::size_t j) const { return values_[i * keys_.size() + j]; }

  // The row index for a NuSIFT ZAI key, or -1 when the matrix does not carry that product.
  // Negative rather than an optional because it indexes into `values_` at every call site.
  int indexOf(std::int64_t key) const;

  // The most negative eigenvalue of the principal submatrix on `rows`, which is <= 0 exactly
  // when that block is not positive semi-definite. Reported rather than repaired -- see the
  // header note. Costs a dense symmetric eigensolve on the submatrix, so it is computed once per
  // answer and not per contributor.
  double smallestEigenvalue(const std::vector<int>& rows) const;

private:
  std::vector<std::int64_t> keys_;
  std::vector<double> values_;  // n*n, row-major
  std::unordered_map<std::int64_t, int> index_;
  YieldCovarianceProvenance provenance_;
};

// Translate a FYCoM product key (Z*10000 + M*1000 + A) into NuSIFT's (Z*10000 + A*10 + I).
// Exposed for the test that pins the collision, and because a reader of the CSV needs it.
std::int64_t yieldKeyFromFycom(std::int64_t fycomKey);

}  // namespace nusift
