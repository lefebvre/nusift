#include "nusift/nucdata/data_store.hpp"

#include <hdf5.h>

#include <string>
#include <vector>

#include "nusift/core/error.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "nucdata store";
constexpr const char* kVersionAttr = "nusift_store_version";

[[noreturn]] void fail(const std::string& what) {
  throw NusiftError(tagged(kModule, what));
}

// RAII for HDF5's integer handles. Every open here is paired with a close even when a read
// throws partway through a file, which matters because a leaked hid_t keeps the file locked
// on Windows and the next open fails with an error that names nothing useful.
class Handle {
public:
  using Closer = herr_t (*)(hid_t);
  Handle(hid_t id, Closer closer) : id_(id), closer_(closer) {}
  ~Handle() {
    if (id_ >= 0) {
      closer_(id_);
    }
  }
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
  Handle(Handle&& other) noexcept : id_(other.id_), closer_(other.closer_) { other.id_ = -1; }

  hid_t get() const { return id_; }
  bool valid() const { return id_ >= 0; }

private:
  hid_t id_;
  Closer closer_;
};

hid_t nativeType(const double*) {
  return H5T_NATIVE_DOUBLE;
}
hid_t nativeType(const int*) {
  return H5T_NATIVE_INT;
}
hid_t nativeType(const std::int64_t*) {
  return H5T_NATIVE_INT64;
}

// Chunked + deflate for anything that can get large. The photon-line and fission-yield arrays
// dominate the file -- tens of thousands of lines, and a yield set per fissile parent per
// energy -- and compressing them is what keeps a full evaluation small enough to ship inside
// a Python wheel. Small arrays are written contiguous, since a chunk header would cost more
// than the compression saves.
template <typename T>
void writeArray(hid_t file, const char* name, const std::vector<T>& values) {
  if (values.empty()) {
    return;  // absent, not zero-filled -- the reader distinguishes the two
  }
  const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
  Handle space(H5Screate_simple(1, dims, nullptr), H5Sclose);
  Handle plist(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
  constexpr hsize_t kCompressAbove = 4096;
  if (dims[0] >= kCompressAbove) {
    const hsize_t chunk[1] = {dims[0] < 65536 ? dims[0] : 65536};
    H5Pset_chunk(plist.get(), 1, chunk);
    H5Pset_deflate(plist.get(), 4);
  }
  const hid_t type = nativeType(static_cast<const T*>(nullptr));
  Handle dset(H5Dcreate2(file, name, type, space.get(), H5P_DEFAULT, plist.get(), H5P_DEFAULT),
              H5Dclose);
  if (!dset.valid()) {
    fail(std::string("could not create dataset ") + name);
  }
  if (H5Dwrite(dset.get(), type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) < 0) {
    fail(std::string("could not write dataset ") + name);
  }
}

// Reads a dataset if present, leaving `out` empty if it is not. Absence is normal: an older
// store predates a field, and a store staged from a source that does not carry photon lines
// legitimately has none.
template <typename T>
void readArray(hid_t file, const char* name, std::vector<T>& out) {
  out.clear();
  if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
    return;
  }
  Handle dset(H5Dopen2(file, name, H5P_DEFAULT), H5Dclose);
  if (!dset.valid()) {
    fail(std::string("could not open dataset ") + name);
  }
  Handle space(H5Dget_space(dset.get()), H5Sclose);
  const int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank != 1) {
    fail(std::string("dataset ") + name + " has rank " + std::to_string(rank) + ", expected 1");
  }
  hsize_t dims[1] = {0};
  H5Sget_simple_extent_dims(space.get(), dims, nullptr);
  out.resize(static_cast<std::size_t>(dims[0]));
  if (dims[0] == 0) {
    return;
  }
  const hid_t type = nativeType(static_cast<const T*>(nullptr));
  if (H5Dread(dset.get(), type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) < 0) {
    fail(std::string("could not read dataset ") + name);
  }
}

void writeIntAttr(hid_t file, const char* name, int value) {
  Handle space(H5Screate(H5S_SCALAR), H5Sclose);
  Handle attr(H5Acreate2(file, name, H5T_NATIVE_INT, space.get(), H5P_DEFAULT, H5P_DEFAULT),
              H5Aclose);
  H5Awrite(attr.get(), H5T_NATIVE_INT, &value);
}

int readIntAttr(hid_t file, const char* name, int fallback) {
  if (H5Aexists(file, name) <= 0) {
    return fallback;
  }
  Handle attr(H5Aopen(file, name, H5P_DEFAULT), H5Aclose);
  int value = fallback;
  H5Aread(attr.get(), H5T_NATIVE_INT, &value);
  return value;
}

void writeStringAttr(hid_t file, const char* name, const std::string& value) {
  Handle type(H5Tcopy(H5T_C_S1), H5Tclose);
  H5Tset_size(type.get(), value.size() + 1);
  H5Tset_strpad(type.get(), H5T_STR_NULLTERM);
  Handle space(H5Screate(H5S_SCALAR), H5Sclose);
  Handle attr(H5Acreate2(file, name, type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
  H5Awrite(attr.get(), type.get(), value.c_str());
}

std::string readStringAttr(hid_t file, const char* name) {
  if (H5Aexists(file, name) <= 0) {
    return {};
  }
  Handle attr(H5Aopen(file, name, H5P_DEFAULT), H5Aclose);
  Handle type(H5Aget_type(attr.get()), H5Tclose);
  const std::size_t size = H5Tget_size(type.get());
  if (size == 0) {
    return {};
  }
  std::string value(size, '\0');
  if (H5Aread(attr.get(), type.get(), value.data()) < 0) {
    return {};
  }
  if (const auto nul = value.find('\0'); nul != std::string::npos) {
    value.resize(nul);
  }
  return value;
}

}  // namespace

void writeStore(const std::string& path, const StoreArrays& a) {
  validateStoreArrays(a);

  Handle file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
  if (!file.valid()) {
    fail("could not create \"" + path + "\"");
  }

  writeIntAttr(file.get(), kVersionAttr, kStoreVersion);
  writeStringAttr(file.get(), "endf_library", a.provenance.library);
  writeStringAttr(file.get(), "created_utc", a.provenance.createdUtc);
  writeStringAttr(file.get(), "nusift_version", a.provenance.nusiftVersion);
  writeIntAttr(file.get(), "staged_tape_count", a.provenance.stagedTapeCount);
  // Per-field rather than per-file, because the ingestion paths cover different fields: a
  // depletion-chain XML carries decay data and yields but no photon lines and no atomic
  // weights. Recording which is which is what lets the exposure metric refuse with a message
  // naming the missing source instead of quietly reporting zeros.
  writeStringAttr(file.get(), "decay_source", dataSourceName(a.provenance.decaySource));
  writeStringAttr(file.get(), "lines_source", dataSourceName(a.provenance.linesSource));
  writeStringAttr(file.get(), "yields_source", dataSourceName(a.provenance.yieldsSource));

  writeArray(file.get(), "nuclide_key", a.nuclideKey);
  writeArray(file.get(), "nuclide_half_life", a.halfLife);
  writeArray(file.get(), "nuclide_awr", a.awr);
  writeArray(file.get(), "nuclide_em_energy_ev", a.emEnergyEv);
  writeArray(file.get(), "nuclide_lp_energy_ev", a.lpEnergyEv);
  writeArray(file.get(), "nuclide_hp_energy_ev", a.hpEnergyEv);
  writeArray(file.get(), "nuclide_continuum_photon_ev", a.continuumPhotonEv);

  writeArray(file.get(), "mode_offset", a.modeOffset);
  writeArray(file.get(), "mode_rtyp", a.modeRtyp);
  writeArray(file.get(), "mode_branching", a.modeBranching);
  writeArray(file.get(), "mode_final_state", a.modeFinalState);
  writeArray(file.get(), "mode_is_fission", a.modeIsFission);

  writeArray(file.get(), "line_offset", a.lineOffset);
  writeArray(file.get(), "line_energy_ev", a.lineEnergyEv);
  writeArray(file.get(), "line_intensity", a.lineIntensity);
  writeArray(file.get(), "line_styp", a.lineStyp);

  writeArray(file.get(), "nfy_parent_key", a.nfyParentKey);
  writeArray(file.get(), "nfy_energy_ev", a.nfyEnergyEv);
  writeArray(file.get(), "nfy_set_offset", a.nfySetOffset);
  writeArray(file.get(), "nfy_product_key", a.nfyProductKey);
  writeArray(file.get(), "nfy_product_yield", a.nfyProductYield);

  writeArray(file.get(), "xs_target_key", a.xsTargetKey);
  writeArray(file.get(), "xs_offset", a.xsOffset);
  writeArray(file.get(), "xs_reaction_type", a.xsReactionType);
  writeArray(file.get(), "xs_product_key", a.xsProductKey);
  writeArray(file.get(), "xs_sigma_barn", a.xsSigmaBarn);
  writeArray(file.get(), "xs_q_ev", a.xsQEv);
  writeArray(file.get(), "xs_energy_ev", a.xsEnergyEv);
  writeArray(file.get(), "xs_spectrum_id", a.xsSpectrumId);
}

StoreArrays readStore(const std::string& path) {
  // HDF5 prints its own error stack to stderr on any failure, which is noise on top of the
  // message thrown here and tells the user nothing they can act on.
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

  if (H5Fis_hdf5(path.c_str()) <= 0) {
    fail("\"" + path + "\" is missing or is not an HDF5 file");
  }
  Handle file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
  if (!file.valid()) {
    fail("could not open \"" + path + "\"");
  }

  StoreArrays a;
  a.provenance.version = readIntAttr(file.get(), kVersionAttr, 0);
  if (a.provenance.version == 0) {
    // Distinguishing "not a NuSIFT store" from "a corrupt NuSIFT store" matters, because the
    // remedy differs: one is the wrong file, the other needs restaging.
    fail("\"" + path + "\" has no " + kVersionAttr +
         " attribute, so it is not a NuSIFT data store");
  }
  if (a.provenance.version > kStoreVersion) {
    fail("\"" + path + "\" is a version " + std::to_string(a.provenance.version) +
         " store, but this build understands at most version " + std::to_string(kStoreVersion) +
         "; upgrade NuSIFT or restage the data");
  }

  const auto source = [&](const char* name) {
    const std::string text = readStringAttr(file.get(), name);
    if (text == "endf") {
      return DataSource::Endf;
    }
    if (text == "openmc-chain-xml") {
      return DataSource::OpenmcChainXml;
    }
    return DataSource::None;
  };
  a.provenance.library = readStringAttr(file.get(), "endf_library");
  a.provenance.createdUtc = readStringAttr(file.get(), "created_utc");
  a.provenance.nusiftVersion = readStringAttr(file.get(), "nusift_version");
  a.provenance.stagedTapeCount = readIntAttr(file.get(), "staged_tape_count", 0);
  a.provenance.decaySource = source("decay_source");
  a.provenance.linesSource = source("lines_source");
  a.provenance.yieldsSource = source("yields_source");

  readArray(file.get(), "nuclide_key", a.nuclideKey);
  readArray(file.get(), "nuclide_half_life", a.halfLife);
  readArray(file.get(), "nuclide_awr", a.awr);
  readArray(file.get(), "nuclide_em_energy_ev", a.emEnergyEv);
  readArray(file.get(), "nuclide_lp_energy_ev", a.lpEnergyEv);
  readArray(file.get(), "nuclide_hp_energy_ev", a.hpEnergyEv);
  readArray(file.get(), "nuclide_continuum_photon_ev", a.continuumPhotonEv);

  readArray(file.get(), "mode_offset", a.modeOffset);
  readArray(file.get(), "mode_rtyp", a.modeRtyp);
  readArray(file.get(), "mode_branching", a.modeBranching);
  readArray(file.get(), "mode_final_state", a.modeFinalState);
  readArray(file.get(), "mode_is_fission", a.modeIsFission);

  readArray(file.get(), "line_offset", a.lineOffset);
  readArray(file.get(), "line_energy_ev", a.lineEnergyEv);
  readArray(file.get(), "line_intensity", a.lineIntensity);
  readArray(file.get(), "line_styp", a.lineStyp);

  readArray(file.get(), "nfy_parent_key", a.nfyParentKey);
  readArray(file.get(), "nfy_energy_ev", a.nfyEnergyEv);
  readArray(file.get(), "nfy_set_offset", a.nfySetOffset);
  readArray(file.get(), "nfy_product_key", a.nfyProductKey);
  readArray(file.get(), "nfy_product_yield", a.nfyProductYield);

  readArray(file.get(), "xs_target_key", a.xsTargetKey);
  readArray(file.get(), "xs_offset", a.xsOffset);
  readArray(file.get(), "xs_reaction_type", a.xsReactionType);
  readArray(file.get(), "xs_product_key", a.xsProductKey);
  readArray(file.get(), "xs_sigma_barn", a.xsSigmaBarn);
  readArray(file.get(), "xs_q_ev", a.xsQEv);
  readArray(file.get(), "xs_energy_ev", a.xsEnergyEv);
  readArray(file.get(), "xs_spectrum_id", a.xsSpectrumId);

  // An empty mode CSR is legal for a store of nothing but stable nuclides, but validate
  // rejects a missing offset array, so normalize it here rather than special-casing later.
  if (a.modeOffset.empty()) {
    a.modeOffset.assign(a.nuclideKey.size() + 1, 0);
  }

  validateStoreArrays(a);
  return a;
}

}  // namespace nusift
