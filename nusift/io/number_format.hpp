#pragma once
/**
 * @file
 * @brief What every writer shares: a double's shortest round-tripping spelling, and a
 *        string escaped for JSON.
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

// Quotes and backslashes are the obvious cases; control characters are the ones that actually
// occur. A provenance string is a file path or a command line the user supplied, and a newline
// or a tab in one produces a report no JSON parser will read -- a failure that surfaces in
// whatever consumes the report rather than here, where it was caused.
inline std::string escapeJson(const std::string& text) {
  std::string out;
  out.reserve(text.size());
  for (const char c : text) {
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\b':
        out += "\\b";
        break;
      case '\f':
        out += "\\f";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default: {
        const auto byte = static_cast<unsigned char>(c);
        if (byte < 0x20) {
          char buffer[8];
          std::snprintf(buffer, sizeof(buffer), "\\u%04x", static_cast<unsigned>(byte));
          out += buffer;
        } else {
          // Everything at 0x20 and above passes through unchanged, including the bytes of a
          // multi-byte UTF-8 sequence: JSON strings are UTF-8, so escaping them would corrupt
          // a nuclide name or a path that is already correct.
          out += c;
        }
        break;
      }
    }
  }
  return out;
}

}  // namespace nusift
