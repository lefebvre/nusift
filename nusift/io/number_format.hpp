#pragma once
/**
 * @file
 * @brief The shortest decimal spelling of a double that reads back as the same value.
 * @ingroup io
 */
#include <cstdio>
#include <cstdlib>
#include <string>

namespace nusift {

// The shortest decimal that reads back as the same double.
//
// Files NuSIFT writes exist to be parsed again -- by itself, by a spreadsheet, by a script --
// so surviving the round trip is the requirement, not a nicety: `--at 1.23456789y` printed at
// six significant digits comes back as a different time than the one the report describes, and
// an inventory written at ten comes back seeding a different atom count than it was given.
//
// The ladder is what keeps the output readable. %.17g always round-trips but renders 0.99999
// as 0.99999000000000005; trying the shorter forms first prints the digits the value actually
// has and falls back only when they are not enough.
inline std::string shortestRoundTrip(double value) {
  char buffer[40];
  for (const int digits : {15, 16, 17}) {
    std::snprintf(buffer, sizeof(buffer), "%.*g", digits, value);
    if (std::strtod(buffer, nullptr) == value) {
      break;
    }
  }
  return buffer;
}

}  // namespace nusift
