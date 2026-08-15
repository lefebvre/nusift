#pragma once
/**
 * @file
 * @brief The exception type every NuSIFT module throws, and the module-tagged helper
 *        that constructs it.
 * @ingroup core
 */
//
// Error handling is uniformly by exception -- no std::expected, no error codes, no
// out-parameters. Libraries throw; main() is the single catch site and maps to an exit
// code (0 success, 1 runtime failure, 2 bad input). The one deliberate exception to
// "throw" is a batch loop over many input files, where a per-file failure is reported to
// stderr and skipped rather than aborting the run.
//
// Every message is prefixed with its module, because a NuSIFT failure is nearly always a
// data or input problem and the user's first question is which stage rejected what:
//
//   nusift: nucdata store: no nusift_store_version attribute in "chain.h5"
//
#include <stdexcept>
#include <string>

namespace nusift {

// Base for every error NuSIFT raises. Deriving from std::runtime_error rather than
// introducing a parallel hierarchy means a consumer that only catches std::exception --
// including the Python binding layer -- still gets a useful message.
class NusiftError : public std::runtime_error {
public:
  explicit NusiftError(const std::string& what) : std::runtime_error(what) {}
};

// Raised for malformed user input (a bad CLI argument, an unparseable inventory row, a
// nuclide name that is not a nuclide name) as opposed to a genuine runtime failure. The
// CLI maps this to exit code 2 and everything else to 1, which is the whole reason the
// distinction exists.
class InputError : public NusiftError {
public:
  explicit InputError(const std::string& what) : NusiftError(what) {}
};

// Construct a module-tagged message. Modules define a local wrapper rather than calling
// this at every throw site:
//
//   namespace {
//   [[noreturn]] void fail(const std::string& what) { throw NusiftError(tagged("nucdata store",
//   what)); }
//   }
//
inline std::string tagged(const char* module, const std::string& what) {
  return std::string(module) + ": " + what;
}

}  // namespace nusift
