// Compiles and links against the installed nusift package only.
//
// The includes matter as much as the assertions: every header pulled in here must have been
// listed in an install(FILES) rule, and none of them may transitively drag in a cram, Eigen,
// or HDF5 header -- those are PRIVATE to the library, and a consumer that is forced to have
// them defeats the PIMPL seam the whole design rests on.
#include <cstdio>
#include <cstdlib>

#include <nusift/core/element_symbols.hpp>
#include <nusift/core/error.hpp>
#include <nusift/core/nuclide.hpp>
#include <nusift/core/nuclide_name.hpp>
#include <nusift/units.hpp>
#include <nusift/version.hpp>

int main() {
  // Header-only paths.
  const nusift::Zai cs137{55, 137, 0};
  if (cs137.key() != 551370) {
    std::fprintf(stderr, "packed key wrong: %lld\n", static_cast<long long>(cs137.key()));
    return 1;
  }

  // A compiled path, which is what actually proves the static library was found and linked
  // rather than merely that the headers were installed.
  const auto parsed = nusift::parseNuclideName("cs137");
  if (!parsed || *parsed != cs137) {
    std::fprintf(stderr, "parseNuclideName failed to link or returned the wrong nuclide\n");
    return 1;
  }
  if (nusift::formatNuclideName(cs137) != "Cs-137") {
    std::fprintf(stderr, "formatNuclideName returned the wrong spelling\n");
    return 1;
  }

  // constexpr paths from the installed headers.
  static_assert(nusift::elementSymbol(55)[0] == 'C', "element table not usable at compile time");
  static_assert(nusift::units::kBqPerCi == 3.7e10, "units header not usable at compile time");

  std::printf("nusift %s -> %s key=%lld\n", nusift::kVersion,
              nusift::formatNuclideName(cs137).c_str(), static_cast<long long>(cs137.key()));
  return 0;
}
