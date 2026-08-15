#pragma once
/**
 * @file
 * @brief Nuclide identity: the (Z, A, I) triple and its packed integer key.
 * @ingroup core
 */
//
// NuSIFT carries its own Zai rather than re-exporting cram::Zai, deliberately. The two are
// layout- and semantically identical, and convert in one line -- but keeping this type
// local means every INSTALLED NuSIFT header depends on nothing but the standard library.
// That is what allows cram (and through it Eigen) to be a PRIVATE link dependency: a
// consumer calling the NuSIFT API never compiles a cram or Eigen header. Conversion happens
// only inside the two translation units that talk to cram at all.
//
// The packed key is the same encoding cram and OpenMC both use:
//
//   key = Z * 10000 + A * 10 + I
//
// so U-235 is 922350, Cs-137 is 551370, and Am-242m is 952421. It is stable across
// versions, sorts sensibly (by element, then mass, then isomeric state), and is what the
// data store persists -- which is why the store never depends on chain ordering.
//
#include <cstdint>
#include <functional>

namespace nusift {

// The widest A and I the packed key can hold: three digits for the mass number and one for the
// isomeric state. A larger value does not overflow, it CARRIES -- A=1000 lands in the Z field
// and the key decodes as a different nuclide entirely, so identity is corrupted rather than
// rejected. Nothing physical comes close (the heaviest evaluated nuclide is near A=300), which
// is exactly why these bounds belong to the validator rather than to the encoding.
inline constexpr int kMaxMassNumber = 999;
inline constexpr int kMaxIsomericState = 9;

// A nuclide: atomic number, mass number, isomeric state (0 = ground, 1 = first metastable).
// Trivially copyable and comparable; safe to pass by value everywhere.
struct Zai {
  int z = 0;
  int a = 0;
  int i = 0;

  constexpr std::int64_t key() const {
    return static_cast<std::int64_t>(z) * 10000 + static_cast<std::int64_t>(a) * 10 + i;
  }

  static constexpr Zai fromKey(std::int64_t key) {
    return Zai{static_cast<int>(key / 10000), static_cast<int>((key / 10) % 1000),
               static_cast<int>(key % 10)};
  }

  // Number of neutrons. Meaningful only for a real nuclide; negative for a malformed one,
  // which is a useful thing for a validator to be able to see.
  constexpr int n() const { return a - z; }

  constexpr bool operator==(const Zai& o) const { return z == o.z && a == o.a && i == o.i; }
  constexpr bool operator!=(const Zai& o) const { return !(*this == o); }
  // Ordering is by key, so a sorted range of nuclides reads as the chart of the nuclides
  // does: by element, then by mass, then ground state before metastable.
  constexpr bool operator<(const Zai& o) const { return key() < o.key(); }
};

// Hash on the packed key. The key is already a dense, collision-free integer encoding, so
// there is nothing to mix -- std::hash<int64_t> is identity on most implementations and
// that is exactly right for keys clustered by element.
struct ZaiHash {
  std::size_t operator()(const Zai& zai) const noexcept {
    return std::hash<std::int64_t>{}(zai.key());
  }
};

// The isobaric mass chain a nuclide belongs to. Named rather than open-coded as `.a`
// because "mass chain" is the domain concept NuSIFT aggregates and ranks by, and the two
// meanings of A diverge the moment a fission-yield or decay-heat weighting is involved.
constexpr int massChain(const Zai& zai) {
  return zai.a;
}

}  // namespace nusift
