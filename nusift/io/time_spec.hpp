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

// --- calendar dates ------------------------------------------------------------
//
// Everything above is a DURATION -- an elapsed span, which is what a cooling time is. An assay
// date is not one: it is an instant on a calendar, and the sheet it comes from writes it as
// one. The two are kept apart deliberately, because "2024-03-15" and "30d" answer different
// questions and a field that quietly accepted either would let a date be read as a span.
//
// Only ISO-8601 is accepted, and only in UTC. Every other spelling is ambiguous somewhere:
// 03/04/2024 is two different days depending on the reader's country, and a local time is a
// different instant depending on where it was written down. An inventory that reconciles two
// assays a day apart cannot afford either.

// Seconds since 1970-01-01T00:00:00Z for "YYYY-MM-DD", or "YYYY-MM-DDThh:mm:ss" with an
// optional trailing "Z". Negative for dates before 1970, which is a legitimate assay date and
// not an error. The result is only ever used in differences, so the choice of origin does not
// reach any answer.
//
// Throws InputError for a malformed spelling and for a date that does not exist -- 2023-02-29
// is rejected rather than rolled forward to March, because a rolled date silently moves an
// assay by a day and nothing downstream could tell.
double parseCalendarDate(std::string_view text);

// Whether `text` is shaped like a calendar date at all. Used to tell a dated inventory column
// from an absent one, and to give a better error than "not a date" for a field that was never
// meant to be one.
bool looksLikeCalendarDate(std::string_view text);

// The inverse: "2024-03-15", or "2024-03-15T09:30:00Z" when the instant is not midnight. What
// a report prints for an epoch, so that the date a reconciled inventory refers to is a date
// rather than a count of seconds nobody can check.
std::string formatCalendarDate(double secondsSinceEpoch);

}  // namespace nusift
