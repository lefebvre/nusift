#pragma once
/**
 * @file
 * @brief The versioned HDF5 nuclear-data store: reading and writing StoreArrays.
 * @ingroup nucdata
 */
//
// The store is what lets production runs rebuild a full depletion chain -- decay data, photon
// spectra, fission yields -- without ENDFtk, without ENDF tapes, and in milliseconds. The
// staging tool reads ENDF and calls writeStore(); everything else calls readStore().
//
// Schema history:
//   v1 -- decay data, discrete photon lines, independent fission yields, atomic weight
//         ratios, MT457 average decay energies, and the continuum photon energy. Reserves
//         the xs_* datasets for one-group activation cross sections, written empty, so
//         adding activation later forces neither a schema bump nor a restage.
//
// A reader accepts any version <= kStoreVersion and leaves fields it does not know empty; a
// newer store is rejected outright rather than read partially. The root attribute is named
// nusift_store_version specifically so that a store from another tool cannot be misread as
// a NuSIFT one, and vice versa.
//
// Datasets are flat 1-D and CSR-packed (see store_arrays.hpp), which keeps the file readable
// by h5py or any other HDF5 tool with no knowledge of NuSIFT's types.
//
#include <string>

#include "nusift/nucdata/store_arrays.hpp"

namespace nusift {

inline constexpr int kStoreVersion = 1;

// Write `arrays` to `path`, truncating any existing file. Validates first, so a malformed
// store cannot be written in the first place.
void writeStore(const std::string& path, const StoreArrays& arrays);

// Read a store. Throws NusiftError for a missing or unreadable file, a file that carries no
// nusift_store_version attribute, or a store written by a newer NuSIFT.
StoreArrays readStore(const std::string& path);

}  // namespace nusift
