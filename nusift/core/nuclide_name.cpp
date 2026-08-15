#include "nusift/core/nuclide_name.hpp"

#include <cctype>
#include <charconv>
#include <string>

#include "nusift/core/element_symbols.hpp"
#include "nusift/core/error.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "nuclide name";

bool isDigit(char c) {
  return c >= '0' && c <= '9';
}
bool isAlpha(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

std::string_view trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())) != 0) {
    s.remove_prefix(1);
  }
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())) != 0) {
    s.remove_suffix(1);
  }
  return s;
}

// Parse a run of digits into `out`. Returns false on overflow or an empty run, so a
// pathological "Cs-99999999999999999999" is rejected rather than wrapping.
bool parseInt(std::string_view s, int& out) {
  if (s.empty()) {
    return false;
  }
  const auto* begin = s.data();
  const auto* end = s.data() + s.size();
  const auto result = std::from_chars(begin, end, out);
  return result.ec == std::errc{} && result.ptr == end;
}

// A nuclide is plausible if the element exists and the mass number can physically hold that
// many protons. This rejects "Cs-0" and "U-1" without needing a table of known nuclides --
// the store is the authority on which nuclides actually have data, and it reports absence
// with a much better message than a parser could.
//
// The upper bound on A is what the packed key can represent, not a physical claim: past three
// digits the mass number carries into the Z field, so "U-1000" would be accepted and then
// persisted and looked up under the key of a nuclide that is not it. A validator that lets an
// identity be silently rewritten is worse than no validator.
bool isPlausible(const Zai& zai) {
  return zai.z >= 1 && zai.z <= kMaxAtomicNumber && zai.a >= zai.z && zai.a <= kMaxMassNumber &&
         zai.i >= 0 && zai.i <= kMaxIsomericState;
}

// Trailing metastable marker: "m", "M", "m1".."m9". Consumes it from `s` and writes the
// isomeric state to `iso`. A bare "m" means state 1, matching ENDF's LISO convention and
// the way every chart of the nuclides prints it.
//
// No marker is not a failure -- most nuclide names are ground state -- so there is nothing
// here for a caller to check, and this returns void rather than a bool that is always true.
void takeIsomerSuffix(std::string_view& s, int& iso) {
  iso = 0;
  if (s.empty()) {
    return;
  }
  // "m2" .. "m9"
  if (s.size() >= 2 && (s[s.size() - 2] == 'm' || s[s.size() - 2] == 'M') && isDigit(s.back())) {
    iso = s.back() - '0';
    s.remove_suffix(2);
    return;
  }
  if (s.back() == 'm' || s.back() == 'M') {
    iso = 1;
    s.remove_suffix(1);
  }
}

}  // namespace

std::optional<Zai> parseNuclideName(std::string_view name) {
  std::string_view s = trim(name);
  if (s.empty()) {
    return std::nullopt;
  }

  // Raw packed key: all digits, e.g. "551370". Checked first because it is unambiguous --
  // no element symbol is all digits.
  bool allDigits = true;
  for (const char c : s) {
    if (!isDigit(c)) {
      allDigits = false;
      break;
    }
  }
  if (allDigits) {
    long long key = 0;
    const auto result = std::from_chars(s.data(), s.data() + s.size(), key);
    if (result.ec != std::errc{} || result.ptr != s.data() + s.size() || key <= 0) {
      return std::nullopt;
    }
    const Zai zai = Zai::fromKey(key);
    return isPlausible(zai) ? std::optional<Zai>(zai) : std::nullopt;
  }

  int iso = 0;
  takeIsomerSuffix(s, iso);

  // Optional leading "<Z>-" as in the IAEA filename convention "092-U-235". Dropped rather
  // than cross-checked against the symbol: if the two disagree the name is malformed, and
  // trusting the symbol matches what a reader would do.
  if (const auto dash = s.find('-'); dash != std::string_view::npos && dash > 0) {
    bool leadingDigits = true;
    for (std::size_t k = 0; k < dash; ++k) {
      if (!isDigit(s[k])) {
        leadingDigits = false;
        break;
      }
    }
    if (leadingDigits) {
      s.remove_prefix(dash + 1);
    }
  }

  // Split the remaining "<symbol>[-]<mass>" at the first digit.
  std::size_t split = 0;
  while (split < s.size() && isAlpha(s[split])) {
    ++split;
  }
  if (split == 0 || split == s.size()) {
    return std::nullopt;
  }
  const std::string_view symbol = s.substr(0, split);
  std::string_view massField = s.substr(split);
  if (!massField.empty() && massField.front() == '-') {
    massField.remove_prefix(1);
  }

  int a = 0;
  if (!parseInt(massField, a)) {
    return std::nullopt;
  }
  const int z = atomicNumber(symbol);
  if (z == 0) {
    return std::nullopt;
  }

  const Zai zai{z, a, iso};
  return isPlausible(zai) ? std::optional<Zai>(zai) : std::nullopt;
}

Zai requireNuclideName(std::string_view name, std::string_view context) {
  if (const auto zai = parseNuclideName(name)) {
    return *zai;
  }
  std::string message;
  if (!context.empty()) {
    message += std::string(context) + ": ";
  }
  message += "not a nuclide name: \"" + std::string(name) + "\"";
  throw InputError(tagged(kModule, message));
}

std::string formatNuclideName(const Zai& zai) {
  std::string out = elementSymbol(zai.z);
  out += '-';
  out += std::to_string(zai.a);
  if (zai.i == 1) {
    out += 'm';
  } else if (zai.i > 1) {
    out += 'm';
    out += std::to_string(zai.i);
  }
  return out;
}

}  // namespace nusift
