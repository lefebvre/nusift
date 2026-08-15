#include "nusift/io/time_spec.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <string>

#include "nusift/core/error.hpp"
#include "nusift/units.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "time";

[[noreturn]] void bad(std::string_view text, const std::string& why) {
  throw InputError(tagged(kModule, "cannot parse time \"" + std::string(text) + "\": " + why));
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

double unitSeconds(char suffix) {
  switch (std::tolower(static_cast<unsigned char>(suffix))) {
    case 's':
      return 1.0;
    case 'm':
      return units::kSecondsPerMinute;
    case 'h':
      return units::kSecondsPerHour;
    case 'd':
      return units::kSecondsPerDay;
    case 'y':
      return units::kSecondsPerYear;
    default:
      return 0.0;
  }
}

// Neither endpoint may be infinite or NaN. A grid is built by interpolating between them, so
// one non-finite endpoint does not produce one bad time -- it poisons every point on the
// grid, and the failure surfaces much later as a CRAM solve against a nonsense time.
void requireFiniteEndpoints(double start, double stop) {
  if (!std::isfinite(start) || !std::isfinite(stop)) {
    throw InputError(tagged(kModule, "a time grid needs finite endpoints"));
  }
}

std::vector<std::string_view> split(std::string_view text, char sep) {
  std::vector<std::string_view> parts;
  std::size_t start = 0;
  while (true) {
    const std::size_t at = text.find(sep, start);
    if (at == std::string_view::npos) {
      parts.push_back(text.substr(start));
      break;
    }
    parts.push_back(text.substr(start, at - start));
    start = at + 1;
  }
  return parts;
}

}  // namespace

double parseDuration(std::string_view text) {
  const std::string_view s = trim(text);
  if (s.empty()) {
    bad(text, "it is empty");
  }

  std::string_view number = s;
  double scale = 1.0;  // a bare number is seconds
  if (const char last = s.back(); std::isalpha(static_cast<unsigned char>(last)) != 0) {
    scale = unitSeconds(last);
    if (scale == 0.0) {
      bad(text, std::string("unknown unit '") + last + "' (use s, m, h, d, or y)");
    }
    number = s.substr(0, s.size() - 1);
  }
  number = trim(number);
  if (number.empty()) {
    bad(text, "it has a unit but no number");
  }

  double value = 0.0;
  const auto result = std::from_chars(number.data(), number.data() + number.size(), value);
  if (result.ec != std::errc{} || result.ptr != number.data() + number.size()) {
    bad(text, "\"" + std::string(number) + "\" is not a number");
  }
  // Finiteness first, so a NaN is reported as what it is rather than as a negative time:
  // from_chars spells "inf" and "nan" as numbers, and neither is a duration anything
  // downstream can do arithmetic with.
  if (!std::isfinite(value)) {
    bad(text, "a time must be finite");
  }
  if (!(value >= 0.0)) {
    bad(text, "a time cannot be negative");
  }
  return value * scale;
}

std::vector<double> logspace(double start, double stop, int count) {
  if (count <= 0) {
    throw InputError(tagged(kModule, "a time grid needs at least one point"));
  }
  requireFiniteEndpoints(start, stop);
  if (!(start > 0.0) || !(stop > 0.0)) {
    throw InputError(tagged(
        kModule, "a log-spaced grid needs positive endpoints; use a linear grid to include 0"));
  }
  if (count == 1) {
    return {start};
  }
  std::vector<double> times(static_cast<std::size_t>(count));
  const double logStart = std::log(start);
  const double logStop = std::log(stop);
  for (int k = 0; k < count; ++k) {
    const double f = static_cast<double>(k) / static_cast<double>(count - 1);
    times[static_cast<std::size_t>(k)] = std::exp(logStart + f * (logStop - logStart));
  }
  // Assigning the endpoints rather than trusting exp(log(x)) to round-trip: a reported time
  // that is 1e-16 off the one the user asked for is confusing in output and, worse, can
  // reorder against an interval endpoint that was meant to coincide with it.
  times.front() = start;
  times.back() = stop;
  return times;
}

std::vector<double> linspace(double start, double stop, int count) {
  if (count <= 0) {
    throw InputError(tagged(kModule, "a time grid needs at least one point"));
  }
  requireFiniteEndpoints(start, stop);
  if (count == 1) {
    return {start};
  }
  std::vector<double> times(static_cast<std::size_t>(count));
  for (int k = 0; k < count; ++k) {
    const double f = static_cast<double>(k) / static_cast<double>(count - 1);
    times[static_cast<std::size_t>(k)] = start + f * (stop - start);
  }
  times.front() = start;
  times.back() = stop;
  return times;
}

std::vector<double> parseTimeGrid(std::string_view text) {
  const std::vector<std::string_view> parts = split(trim(text), ':');
  if (parts.size() != 4) {
    bad(text, "a grid is start:stop:log|lin:count, e.g. 1h:100y:log:60");
  }

  const double start = parseDuration(parts[0]);
  const double stop = parseDuration(parts[1]);
  if (stop <= start) {
    bad(text, "the stop time must be after the start time");
  }

  std::string spacing(parts[2]);
  std::transform(spacing.begin(), spacing.end(), spacing.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  int count = 0;
  const auto result = std::from_chars(parts[3].data(), parts[3].data() + parts[3].size(), count);
  if (result.ec != std::errc{} || result.ptr != parts[3].data() + parts[3].size()) {
    bad(text, "\"" + std::string(parts[3]) + "\" is not a point count");
  }
  if (count < 2) {
    bad(text, "a grid needs at least 2 points");
  }

  if (spacing == "log") {
    return logspace(start, stop, count);
  }
  if (spacing == "lin") {
    return linspace(start, stop, count);
  }
  bad(text, "spacing must be log or lin, not \"" + spacing + "\"");
}

std::vector<double> mergeTimes(std::vector<double> times) {
  times.erase(std::remove_if(times.begin(), times.end(), [](double t) { return !(t >= 0.0); }),
              times.end());
  std::sort(times.begin(), times.end());
  // Collapse times that differ only in the last few bits. The engine rejects exact
  // duplicates outright, and a pair separated by 1e-16 of a second costs a full
  // factorization while carrying no information the neighbouring point does not.
  times.erase(std::unique(times.begin(), times.end(),
                          [](double a, double b) {
                            return std::abs(b - a) <=
                                   1e-12 * std::max({std::abs(a), std::abs(b), 1.0});
                          }),
              times.end());
  return times;
}

std::string formatDuration(double seconds) {
  struct Unit {
    double size;
    const char* suffix;
  };
  // Descending, so the first unit the value reaches is the largest that keeps it readable.
  static constexpr Unit kUnits[] = {
      {units::kSecondsPerYear, "y"},
      {units::kSecondsPerDay, "d"},
      {units::kSecondsPerHour, "h"},
      {units::kSecondsPerMinute, "m"},
      {1.0, "s"},
  };
  // One format for every magnitude, so a column of times lines up: "%g" drops trailing
  // zeros, giving "30 d" and "1.5 y" rather than "30.00 d" beside "1.50 y". Four significant
  // figures is more than enough for a cooling time and keeps the column narrow.
  char buffer[64];
  for (const Unit& unit : kUnits) {
    if (seconds >= unit.size) {
      std::snprintf(buffer, sizeof(buffer), "%.4g %s", seconds / unit.size, unit.suffix);
      return buffer;
    }
  }
  std::snprintf(buffer, sizeof(buffer), "%.4g s", seconds);
  return buffer;
}

}  // namespace nusift
