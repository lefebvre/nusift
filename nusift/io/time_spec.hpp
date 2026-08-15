#pragma once
/**
 * @file
 * @brief Parsing of human-written times and time grids.
 * @ingroup io
 */
//
// Times arrive from a CLI argument or a config file as strings, and nobody writes cooling
// times in seconds. Accepted forms:
//
//   30d  1.5y  90m  12h  3600s  3600      a duration, bare numbers being seconds
//   1h:100y:log:60                        a grid: start, stop, spacing, count
//   0:1d:lin:25
//
// A YEAR IS 365.25 DAYS, the Julian year. The choice is arbitrary but its consequences are
// not -- over a 100 y decay, a 365 d year differs by more than a year of elapsed time -- so
// it is stated here, in --help, and in the report header.
//
#include <string>
#include <string_view>
#include <vector>

namespace nusift {

// Seconds for a single duration. Throws InputError naming the token if it is unparseable,
// negative, or non-finite -- "inf" and "nan" parse as numbers and are not durations.
double parseDuration(std::string_view text);

// A grid "start:stop:log|lin:count". Both endpoints are hit exactly -- a log grid that
// missed its endpoints would put the reported times somewhere other than where the user
// asked, which matters when one of them is a regulatory decision point.
//
// A log grid requires a positive start, since log spacing from zero is undefined; the error
// says so rather than silently substituting a small number.
std::vector<double> parseTimeGrid(std::string_view text);

// Log-spaced times from `start` to `stop` inclusive, `count` points. count == 1 yields
// {start}; both endpoints are exact rather than accumulated by repeated multiplication.
//
// Both endpoints must be finite and positive. A non-finite endpoint does not produce one bad
// time but a grid of them, since every point is interpolated between the two.
std::vector<double> logspace(double start, double stop, int count);

// Linearly spaced, same endpoint guarantee, and the same requirement that both endpoints be
// finite. Unlike a log grid this one may start at zero.
std::vector<double> linspace(double start, double stop, int count);

// Merge, sort, and de-duplicate times, dropping any that are negative. Duplicates would make
// the engine reject the set outright, and near-duplicates that differ only in the last bit
// cost a full factorization for no information.
std::vector<double> mergeTimes(std::vector<double> times);

// Render seconds the way a person would write them, choosing the largest unit that keeps the
// number readable: "30 d", "1.5 y", "45 m". Used in every report.
std::string formatDuration(double seconds);

}  // namespace nusift
