#include "nusift/exposure/dose_coefficients.hpp"

#include <cctype>
#include <cmath>
#include <cstddef>

namespace nusift::exposure {
namespace {

// ICRP Publication 116, Table A.1 -- effective dose per fluence, pSv*cm^2, for monoenergetic
// photons in the six irradiation geometries, on ICRP's own energy grid in MeV.
//
// Transcribed exactly as the publication prints it, digits and all: 0.0390 rather than 0.039,
// so the table can be checked against the page line by line. The MeV and the pSv*cm^2 are
// converted once, in the accessor below, for the same reason -- a table that has been through
// a unit conversion is no longer a transcription of anything.
//
// The grid is ICRP's, not a resampling, and it carries the energies that matter to this tool:
// 0.511, 0.662, 1.117 and 1.33 MeV are all tabulated points rather than interpolations.
struct DoseCoefficient {
  double energyMeV;
  double perFluencePSvCm2[6];  // indexed by Irradiation
};

// clang-format off
// The alignment is the point here, not a preference: these rows are meant to be read against
// the page they were transcribed from, and a reflowed table cannot be.
constexpr DoseCoefficient kEffectiveDose[] = {
    { 0.01, {0.0685, 0.0184, 0.0189, 0.0182, 0.0337, 0.0288}},
    {0.015, { 0.156, 0.0155, 0.0416, 0.0390, 0.0664, 0.0560}},
    { 0.02, { 0.225, 0.0260, 0.0655, 0.0573, 0.0986, 0.0812}},
    { 0.03, { 0.313, 0.0940,  0.110, 0.0891,  0.158,  0.127}},
    { 0.04, { 0.351,  0.161,  0.140,  0.114,  0.199,  0.158}},
    { 0.05, { 0.370,  0.208,  0.160,  0.133,  0.226,  0.180}},
    { 0.06, { 0.390,  0.242,  0.177,  0.150,  0.248,  0.199}},
    { 0.07, { 0.413,  0.271,  0.194,  0.167,  0.273,  0.218}},
    { 0.08, { 0.444,  0.301,  0.214,  0.185,  0.297,  0.239}},
    {  0.1, { 0.519,  0.361,  0.259,  0.225,  0.355,  0.287}},
    { 0.15, { 0.748,  0.541,  0.395,  0.348,  0.528,  0.429}},
    {  0.2, {  1.00,  0.741,  0.552,  0.492,  0.721,  0.589}},
    {  0.3, {  1.51,   1.16,  0.888,  0.802,   1.12,  0.932}},
    {  0.4, {  2.00,   1.57,   1.24,   1.13,   1.52,   1.28}},
    {  0.5, {  2.47,   1.98,   1.58,   1.45,   1.92,   1.63}},
    {0.511, {  2.52,   2.03,   1.62,   1.49,   1.96,   1.67}},
    {  0.6, {  2.91,   2.38,   1.93,   1.78,   2.30,   1.97}},
    {0.662, {  3.17,   2.62,   2.14,   1.98,   2.54,   2.17}},
    {  0.8, {  3.73,   3.13,   2.59,   2.41,   3.04,   2.62}},
    {  1.0, {  4.49,   3.83,   3.23,   3.03,   3.72,   3.25}},
    {1.117, {  4.90,   4.22,   3.58,   3.37,   4.10,   3.60}},
    { 1.33, {  5.59,   4.89,   4.20,   3.98,   4.75,   4.20}},
    {  1.5, {  6.12,   5.39,   4.68,   4.45,   5.24,   4.66}},
    {  2.0, {  7.48,   6.75,   5.96,   5.70,   6.55,   5.90}},
    {  3.0, {  9.75,   9.12,   8.21,   7.90,   8.84,   8.08}},
    {  4.0, {  11.7,   11.2,   10.2,   9.86,   10.8,   10.0}},
    {  5.0, {  13.4,   13.1,   12.0,   11.7,   12.7,   11.8}},
    {  6.0, {  15.0,   15.0,   13.7,   13.4,   14.4,   13.5}},
    {6.129, {  15.1,   15.2,   13.9,   13.6,   14.6,   13.7}},
    {  8.0, {  17.8,   18.6,   17.0,   16.6,   17.6,   16.6}},
    { 10.0, {  20.5,   22.0,   20.1,   19.7,   20.6,   19.6}},
    { 15.0, {  26.1,   30.3,   27.4,   27.1,   27.7,   26.8}},
    { 20.0, {  30.8,   38.2,   34.4,   34.4,   34.4,   33.8}},
    { 30.0, {  37.9,   51.4,   47.4,   48.1,   46.1,   46.1}},
    { 40.0, {  43.1,   62.0,   59.2,   60.9,   56.0,   56.9}},
    { 50.0, {  47.1,   70.4,   69.5,   72.2,   64.4,   66.2}},
    { 60.0, {  50.1,   76.9,   78.3,   82.0,   71.2,   74.1}},
    { 80.0, {  54.5,   86.6,   92.4,   97.9,   82.0,   87.2}},
    {  100, {  57.8,   93.2,    103,    110,   89.7,   97.5}},
    {  150, {  63.3,    104,    121,    130,    102,    116}},
    {  200, {  67.3,    111,    133,    143,    111,    130}},
    {  300, {  72.3,    119,    148,    161,    121,    147}},
    {  400, {  75.5,    124,    158,    172,    128,    159}},
    {  500, {  77.5,    128,    165,    180,    133,    168}},
    {  600, {  78.9,    131,    170,    186,    136,    174}},
    {  800, {  80.5,    135,    178,    195,    142,    185}},
    { 1000, {  81.7,    138,    183,    201,    145,    193}},
    { 1500, {  83.8,    142,    193,    212,    152,    208}},
    { 2000, {  85.2,    145,    198,    220,    156,    218}},
    { 3000, {  86.9,    148,    206,    229,    161,    232}},
    { 4000, {  88.1,    150,    212,    235,    165,    243}},
    { 5000, {  88.9,    152,    216,    240,    168,    251}},
    { 6000, {  89.5,    153,    219,    244,    170,    258}},
    { 8000, {  90.2,    155,    224,    251,    172,    268}},
    {10000, {  90.7,    155,    228,    255,    175,    276}},
};
// clang-format on
constexpr std::size_t kCount = sizeof(kEffectiveDose) / sizeof(kEffectiveDose[0]);

// pSv*cm^2 to Sv*m^2: 1e-12 Sv, and 1e-4 m^2 per cm^2.
constexpr double kPSvCm2ToSvM2 = 1.0e-16;
constexpr double kEvPerMeV = 1.0e6;

}  // namespace

const char* irradiationName(Irradiation irradiation) {
  switch (irradiation) {
    case Irradiation::AP:
      return "AP";
    case Irradiation::PA:
      return "PA";
    case Irradiation::LLAT:
      return "LLAT";
    case Irradiation::RLAT:
      return "RLAT";
    case Irradiation::ROT:
      return "ROT";
    case Irradiation::ISO:
      return "ISO";
  }
  return "?";
}

bool parseIrradiation(std::string_view text, Irradiation& out) {
  constexpr Irradiation kAll[] = {Irradiation::AP,   Irradiation::PA,  Irradiation::LLAT,
                                  Irradiation::RLAT, Irradiation::ROT, Irradiation::ISO};
  for (const Irradiation candidate : kAll) {
    const std::string_view name = irradiationName(candidate);
    if (name.size() != text.size()) {
      continue;
    }
    bool same = true;
    for (std::size_t i = 0; i < name.size(); ++i) {
      const char a = static_cast<char>(std::toupper(static_cast<unsigned char>(name[i])));
      const char b = static_cast<char>(std::toupper(static_cast<unsigned char>(text[i])));
      if (a != b) {
        same = false;
        break;
      }
    }
    if (same) {
      out = candidate;
      return true;
    }
  }
  return false;
}

double effectiveDosePerFluence(double energyEv, Irradiation irradiation) {
  const std::size_t column = static_cast<std::size_t>(irradiation);
  const auto value = [column](std::size_t i) {
    return kEffectiveDose[i].perFluencePSvCm2[column] * kPSvCm2ToSvM2;
  };
  const double energyMeV = energyEv / kEvPerMeV;

  // Clamped at both ends rather than extrapolated, exactly as the air coefficients are, and
  // for the same reason: outside the table the curve's shape is not something this file knows.
  if (energyMeV <= kEffectiveDose[0].energyMeV) {
    return value(0);
  }
  if (energyMeV >= kEffectiveDose[kCount - 1].energyMeV) {
    return value(kCount - 1);
  }

  std::size_t hi = 1;
  while (hi < kCount && kEffectiveDose[hi].energyMeV < energyMeV) {
    ++hi;
  }
  const std::size_t lo = hi - 1;
  // Log-log, which is what ICRP itself recommends for interpolating these tables, and what the
  // shape asks for: the coefficient rises by three orders of magnitude across the grid and is
  // close to a power law over any one interval of it.
  const double t =
      (std::log(energyMeV) - std::log(kEffectiveDose[lo].energyMeV)) /
      (std::log(kEffectiveDose[hi].energyMeV) - std::log(kEffectiveDose[lo].energyMeV));
  return std::exp(std::log(value(lo)) + t * (std::log(value(hi)) - std::log(value(lo))));
}

bool isOutsideTabulatedDoseRange(double energyEv) {
  return energyEv < kMinTabulatedDoseEv || energyEv > kMaxTabulatedDoseEv;
}

}  // namespace nusift::exposure
