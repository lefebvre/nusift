# ---------------------------------------------------------------------------
# Single source of truth for every external dependency NuSIFT pins. Nothing else in
# the tree names a repository or a tag, so a bump happens in exactly one place and
# cannot drift between the top-level build, a CI prebuild, and the wheel build.
#
# NOT here: Eigen (reached transitively through cram, which owns its own
# CRAM_EIGEN_* selection) and HDF5 (system- or vcpkg-provided; there is no sane
# source build to pin).
# ---------------------------------------------------------------------------

# --- cram-depletion ---------------------------------------------------------
# Two independent knobs, because cram is consumed two ways.
#
# NUSIFT_CRAM_MIN_VERSION is the find_package() request against an INSTALLED cram.
# cram's package version file uses SameMajorVersion compatibility, so "2.0" is
# satisfied by any 2.x at or above it -- and NOT by a 1.x. The two knobs therefore have
# to cross a major boundary together: leaving the floor at 1.0 while fetching 2.0.0
# would not fail, it would quietly accept an installed cram whose headers no longer
# match the ones the fetch path compiles against.
#
# NUSIFT_CRAM_VERSION is the git tag fetched when no installed cram is found. It is a
# RELEASE TAG, never a branch: the CRAM solver is the numerical core of every number
# NuSIFT reports, so a floating dependency could silently shift the golden baselines.
#
# 2.0.0 is a major bump for a source-compatibility break rather than a behavioral one:
# cram dropped the default member initializers on Zai, DecayMode, DecayData::halfLife,
# and FissionYields::energy, so that -Wmissing-field-initializers demands each field at
# every braced initialization instead of defaulting it to a zero that reads as a valid
# nuclide. NuSIFT supplies all of them (see nusift/nucdata/nuclear_data.cpp); the
# solver results are unchanged.
#
# cram's burnup API (DepletionSystem, Integrator, loadDepletionChainXml) is still not
# tagged and nothing here depends on it yet. If activation work begins while that is
# still untagged, pin an immutable commit SHA here -- never a branch name -- and set
# GIT_SHALLOW FALSE, which a bare SHA requires.
set(NUSIFT_CRAM_REPO        "https://github.com/lefebvre/cram-depletion.git")
set(NUSIFT_CRAM_VERSION     "v2.0.0" CACHE STRING
    "cram-depletion release tag to fetch when no installed cram is found")
set(NUSIFT_CRAM_MIN_VERSION "2.0" CACHE STRING
    "Minimum acceptable version of an installed cram-depletion package")

# --- CLI11: command-line parsing for the nusift driver ----------------------
set(NUSIFT_CLI11_REPO   "https://github.com/CLIUtils/CLI11.git")
set(NUSIFT_CLI11_TAG    "v2.4.0")

# --- nlohmann/json: machine-readable report output --------------------------
set(NUSIFT_JSON_REPO    "https://github.com/nlohmann/json.git")
set(NUSIFT_JSON_TAG     "v3.11.3")

# --- toml++: the run-config format ------------------------------------------
# TOML rather than JSON for config specifically because a config file is hand-authored
# and diffed: comments, unquoted keys, and no trailing-comma trap. JSON stays the
# OUTPUT format, where nested numeric arrays are the shape and every consumer reads it.
set(NUSIFT_TOMLPP_REPO  "https://github.com/marzer/tomlplusplus.git")
set(NUSIFT_TOMLPP_TAG   "v3.4.0")

# --- GoogleTest -------------------------------------------------------------
# Newer than the v1.15.2 cram pins for its own suite. There is no conflict: NuSIFT
# forces CRAM_ENABLE_TESTS=OFF in the sub-build, so cram never fetches GoogleTest at
# all. Flipping that option on would reintroduce the clash -- don't.
set(NUSIFT_GTEST_REPO   "https://github.com/google/googletest.git")
set(NUSIFT_GTEST_TAG    "v1.17.0")

# --- nanobind: the Python extension -----------------------------------------
# Found via find_package(nanobind CONFIG) from the pip-installed package rather than
# fetched, which is how scikit-build-core expects it to be located.
set(NUSIFT_NANOBIND_MIN_VERSION "2.4.0")
