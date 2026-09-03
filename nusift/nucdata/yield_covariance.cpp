#include "nusift/nucdata/yield_covariance.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <ios>
#include <sstream>
#include <string_view>
#include <system_error>

#include <Eigen/Dense>

#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "yield covariance";

// A correlation is symmetric by construction, but it arrives through fifteen decimal digits of
// text. This is a rounding tolerance, not a repair: anything larger means the file is not a
// symmetric matrix and the caller is owed an error rather than an average.
constexpr double kSymmetryTolerance = 1e-9;
// The diagonal of a correlation matrix is one. The same reasoning applies.
constexpr double kDiagonalTolerance = 1e-9;

[[noreturn]] void fail(const std::string& what) {
  throw InputError(tagged(kModule, what));
}

std::string_view trimmed(std::string_view text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string_view::npos) {
    return {};
  }
  return text.substr(first, text.find_last_not_of(" \t\r\n") - first + 1);
}

double requireNumber(std::string_view text, const std::string& where) {
  const std::string_view value = trimmed(text);
  double out = 0.0;
  const auto result = std::from_chars(value.data(), value.data() + value.size(), out);
  if (result.ec != std::errc{} || result.ptr != value.data() + value.size()) {
    fail(where + ": \"" + std::string(value) + "\" is not a number");
  }
  // from_chars spells "inf" and "nan" as numbers, and a correlation is neither.
  if (!std::isfinite(out)) {
    fail(where + ": \"" + std::string(value) + "\" is not a finite correlation");
  }
  return out;
}

std::int64_t requireKey(std::string_view text, const std::string& where) {
  const std::string_view value = trimmed(text);
  std::int64_t out = 0;
  const auto result = std::from_chars(value.data(), value.data() + value.size(), out);
  if (result.ec != std::errc{} || result.ptr != value.data() + value.size()) {
    fail(where + ": \"" + std::string(value) + "\" is not a product key");
  }
  return yieldKeyFromFycom(out);
}

// Split one CSV line into fields, in place over the buffer. The matrices run to a thousand
// columns and a thousand rows, so this is called a million times and does no allocation.
void splitFields(std::string_view line, std::vector<std::string_view>& out) {
  out.clear();
  std::size_t start = 0;
  while (true) {
    const std::size_t comma = line.find(',', start);
    if (comma == std::string_view::npos) {
      out.push_back(line.substr(start));
      return;
    }
    out.push_back(line.substr(start, comma - start));
    start = comma + 1;
  }
}

// The system code a FYCoM filename carries: matrices/ENDF/independent/U235T_corr.csv is U235T.
// Recovered rather than required from the user, because it is already in the path and a
// mistyped one would mislabel every answer.
std::string systemFromPath(const std::string& path) {
  std::size_t begin = path.find_last_of("/\\");
  begin = begin == std::string::npos ? 0 : begin + 1;
  const std::size_t underscore = path.find('_', begin);
  if (underscore == std::string::npos) {
    return {};
  }
  return path.substr(begin, underscore - begin);
}

}  // namespace

std::int64_t yieldKeyFromFycom(std::int64_t fycomKey) {
  const std::int64_t z = fycomKey / 10000;
  const std::int64_t rest = fycomKey % 10000;
  const std::int64_t isomer = rest / 1000;
  const std::int64_t a = rest % 1000;
  return Zai{static_cast<int>(z), static_cast<int>(a), static_cast<int>(isomer)}.key();
}

int YieldCorrelation::indexOf(std::int64_t key) const {
  const auto it = index_.find(key);
  return it == index_.end() ? -1 : it->second;
}

YieldCorrelation YieldCorrelation::read(const std::string& path, const std::string& library) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    fail("cannot open \"" + path + "\"");
  }
  std::ostringstream buffer;
  buffer << file.rdbuf();
  const std::string text = buffer.str();

  YieldCorrelation matrix;
  matrix.provenance_.path = path;
  matrix.provenance_.system = systemFromPath(path);
  matrix.provenance_.library = library;
  matrix.provenance_.citation =
      "Matthews, Bernstein & Younes, At. Data Nucl. Data Tables 140 (2021) 101441; "
      "matrices doi:10.5281/zenodo.5985936";

  std::size_t line = 0;
  std::size_t begin = 0;
  std::vector<std::string_view> fields;
  std::size_t row = 0;
  std::size_t n = 0;

  while (begin <= text.size()) {
    const std::size_t end = text.find('\n', begin);
    const std::string_view raw = std::string_view(text).substr(
        begin, end == std::string::npos ? std::string::npos : end - begin);
    begin = end == std::string::npos ? text.size() + 1 : end + 1;
    ++line;
    if (trimmed(raw).empty()) {
      continue;
    }
    const std::string where = "\"" + path + "\" line " + std::to_string(line);
    splitFields(raw, fields);

    // The header row: an empty leading cell, then one key per column. It fixes the order of
    // everything below it, and the row keys are checked against it rather than trusted.
    if (n == 0) {
      if (fields.size() < 2) {
        fail(where + ": expected a header row of product keys, found " +
             std::to_string(fields.size()) + " field(s). This reads a FYCoM correlation CSV");
      }
      if (!trimmed(fields.front()).empty()) {
        fail(where + ": a FYCoM header row begins with an empty cell, found \"" +
             std::string(trimmed(fields.front())) + "\"");
      }
      n = fields.size() - 1;
      matrix.keys_.reserve(n);
      for (std::size_t c = 1; c < fields.size(); ++c) {
        const std::int64_t key = requireKey(fields[c], where);
        if (!matrix.index_.emplace(key, static_cast<int>(matrix.keys_.size())).second) {
          fail(where + ": product key " + std::string(trimmed(fields[c])) +
               " appears in two columns");
        }
        matrix.keys_.push_back(key);
      }
      matrix.values_.assign(n * n, 0.0);
      continue;
    }

    if (row >= n) {
      fail(where + ": " + std::to_string(row + 1) + " rows for " + std::to_string(n) +
           " columns; a correlation matrix is square");
    }
    if (fields.size() != n + 1) {
      fail(where + ": " + std::to_string(fields.size() - 1) + " values for " + std::to_string(n) +
           " columns");
    }
    // The row label has to be the column label at the same position. A file whose rows and
    // columns are ordered differently would contract silently against the wrong products.
    const std::int64_t key = requireKey(fields.front(), where);
    if (key != matrix.keys_[row]) {
      fail(where + ": row " + std::to_string(row + 1) + " is keyed " +
           std::string(trimmed(fields.front())) +
           ", which is not the key in the same position of the header row. Rows and columns must "
           "carry the same products in the same order");
    }
    for (std::size_t c = 0; c < n; ++c) {
      matrix.values_[row * n + c] = requireNumber(fields[c + 1], where);
    }
    ++row;
  }

  if (n == 0) {
    fail("\"" + path + "\" is empty");
  }
  if (row != n) {
    fail("\"" + path + "\": " + std::to_string(row) + " rows for " + std::to_string(n) +
         " columns; a correlation matrix is square");
  }

  for (std::size_t i = 0; i < n; ++i) {
    if (std::abs(matrix.values_[i * n + i] - 1.0) > kDiagonalTolerance) {
      fail("\"" + path + "\": the diagonal entry for product " + std::to_string(matrix.keys_[i]) +
           " is " + std::to_string(matrix.values_[i * n + i]) +
           ", not 1. This reads a CORRELATION matrix -- a *_corr.csv, not a *_cov.csv, because "
           "only the correlation is imported and the variances come from the store");
    }
    for (std::size_t j = i + 1; j < n; ++j) {
      const double asymmetry = std::abs(matrix.values_[i * n + j] - matrix.values_[j * n + i]);
      if (asymmetry > kSymmetryTolerance) {
        fail("\"" + path + "\": entries (" + std::to_string(i + 1) + "," + std::to_string(j + 1) +
             ") and its transpose differ by " + std::to_string(asymmetry) +
             "; a correlation matrix is symmetric");
      }
    }
  }
  return matrix;
}

double YieldCorrelation::smallestEigenvalue(const std::vector<int>& rows) const {
  if (rows.empty()) {
    return 0.0;
  }
  const auto n = static_cast<Eigen::Index>(rows.size());
  Eigen::MatrixXd block(n, n);
  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = 0; j < n; ++j) {
      block(i, j) = at(static_cast<std::size_t>(rows[static_cast<std::size_t>(i)]),
                       static_cast<std::size_t>(rows[static_cast<std::size_t>(j)]));
    }
  }
  // Symmetrized before the solve for the same reason the read tolerates a rounding asymmetry:
  // SelfAdjointEigenSolver reads one triangle anyway, and this makes which one irrelevant.
  block = 0.5 * (block + block.transpose()).eval();
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(block, Eigen::EigenvaluesOnly);
  if (solver.info() != Eigen::Success) {
    fail("the symmetric eigensolve on the restricted correlation block did not converge");
  }
  return solver.eigenvalues().minCoeff();
}

}  // namespace nusift
