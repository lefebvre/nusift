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

namespace {

constexpr const char* kDateModule = "date";

[[noreturn]] void badDate(std::string_view text, const std::string& why) {
  throw InputError(tagged(kDateModule, "cannot read date \"" + std::string(text) + "\": " + why));
}

bool isLeap(long long y) {
  return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
}

int daysInMonth(long long y, int m) {
  static constexpr int kDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if (m == 2 && isLeap(y)) {
    return 29;
  }
  return kDays[m - 1];
}

// Days since 1970-01-01 in the proleptic Gregorian calendar -- Howard Hinnant's civil-date
// algorithm, which is exact over the whole range of a signed 64-bit day count and needs no
// library, no locale and no timezone database. Written out rather than delegated to
// std::chrono's calendar because the arithmetic is ten lines and the dependency is not.
long long daysFromCivil(long long y, int m, int d) {
  y -= m <= 2 ? 1 : 0;
  const long long era = (y >= 0 ? y : y - 399) / 400;
  const long long yoe = y - era * 400;                                   // [0, 399]
  const long long doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;  // [0, 365]
  const long long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;           // [0, 146096]
  return era * 146097 + doe - 719468;
}

void civilFromDays(long long z, long long& y, int& m, int& d) {
  z += 719468;
  const long long era = (z >= 0 ? z : z - 146096) / 146097;
  const long long doe = z - era * 146097;                                       // [0, 146096]
  const long long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;  // [0, 399]
  const long long yr = yoe + era * 400;
  const long long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);  // [0, 365]
  const long long mp = (5 * doy + 2) / 153;                       // [0, 11]
  d = static_cast<int>(doy - (153 * mp + 2) / 5 + 1);             // [1, 31]
  m = static_cast<int>(mp + (mp < 10 ? 3 : -9));                  // [1, 12]
  y = yr + (m <= 2 ? 1 : 0);
}

// Exactly `width` digits, no sign, no spaces. from_chars would accept "+5" and " 5" and stop
// early on "2024-03", none of which is a field of a fixed-width date.
bool digits(std::string_view text, std::size_t at, std::size_t width, long long& out) {
  if (at + width > text.size()) {
    return false;
  }
  long long value = 0;
  for (std::size_t i = 0; i < width; ++i) {
    const unsigned char c = static_cast<unsigned char>(text[at + i]);
    if (std::isdigit(c) == 0) {
      return false;
    }
    value = value * 10 + (c - '0');
  }
  out = value;
  return true;
}

}  // namespace

bool looksLikeCalendarDate(std::string_view text) {
  const std::string_view t = trim(text);
  long long ignored = 0;
  return t.size() >= 10 && digits(t, 0, 4, ignored) && t[4] == '-' && digits(t, 5, 2, ignored) &&
         t[7] == '-' && digits(t, 8, 2, ignored);
}

double parseCalendarDate(std::string_view text) {
  const std::string_view t = trim(text);
  long long year = 0;
  long long month = 0;
  long long day = 0;
  if (t.size() < 10 || !digits(t, 0, 4, year) || t[4] != '-' || !digits(t, 5, 2, month) ||
      t[7] != '-' || !digits(t, 8, 2, day)) {
    badDate(text, "expected YYYY-MM-DD, optionally followed by Thh:mm:ss");
  }
  if (month < 1 || month > 12) {
    badDate(text, "there is no month " + std::to_string(month));
  }
  const int lastDay = daysInMonth(year, static_cast<int>(month));
  if (day < 1 || day > lastDay) {
    badDate(text, std::to_string(year) + "-" + (month < 10 ? "0" : "") + std::to_string(month) +
                      " has " + std::to_string(lastDay) + " days");
  }

  long long hour = 0;
  long long minute = 0;
  long long second = 0;
  if (t.size() > 10) {
    std::string_view rest = t.substr(10);
    // A trailing Z is accepted and means what the whole format already means. Any other
    // offset is refused rather than ignored: silently reading +05:00 as UTC would move the
    // assay by five hours, and an inventory reconciled across two of them would be wrong by
    // exactly that much.
    if (!rest.empty() && (rest.back() == 'Z' || rest.back() == 'z')) {
      rest.remove_suffix(1);
    }
    if (rest.empty()) {
      badDate(text, "a date ending in Z still needs its time, as Thh:mm:ss");
    }
    if (rest.front() != 'T' && rest.front() != 't' && rest.front() != ' ') {
      badDate(text, "the time is separated from the date by T");
    }
    rest.remove_prefix(1);
    if (rest.size() != 8 || !digits(rest, 0, 2, hour) || rest[2] != ':' ||
        !digits(rest, 3, 2, minute) || rest[5] != ':' || !digits(rest, 6, 2, second)) {
      badDate(text, "expected a time as hh:mm:ss");
    }
    // 60 is a leap second, and NuSIFT does not model them; 61 is nothing at all.
    if (hour > 23 || minute > 59 || second > 59) {
      badDate(text, "there is no time " + std::string(rest));
    }
  }

  const long long days = daysFromCivil(year, static_cast<int>(month), static_cast<int>(day));
  return static_cast<double>(days) * units::kSecondsPerDay +
         static_cast<double>(hour * 3600 + minute * 60 + second);
}

std::string formatCalendarDate(double secondsSinceEpoch) {
  // Floor rather than truncate, so a date before 1970 lands on the day that contains it
  // rather than on the one after.
  const double dayCount = std::floor(secondsSinceEpoch / units::kSecondsPerDay);
  long long days = static_cast<long long>(dayCount);
  int within = static_cast<int>(std::llround(secondsSinceEpoch - dayCount * units::kSecondsPerDay));
  // Rounding an instant a fraction of a second before midnight lands on 86400, which is not a
  // time of day. Carrying it into the next date is the truthful reading; clamping it to
  // 23:59:59 would move the instant backwards by a whole day's worth of naming.
  if (within >= 86400) {
    within -= 86400;
    ++days;
  }
  if (within < 0) {
    within = 0;
  }

  long long year = 0;
  int month = 0;
  int day = 0;
  civilFromDays(days, year, month, day);

  const int hour = within / 3600;
  const int minute = (within / 60) % 60;
  const int second = within % 60;

  char buffer[96];
  if (within == 0) {
    std::snprintf(buffer, sizeof(buffer), "%04lld-%02d-%02d", year, month, day);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%04lld-%02d-%02dT%02d:%02d:%02dZ", year, month, day,
                  hour, minute, second);
  }
  return buffer;
}

}  // namespace nusift
