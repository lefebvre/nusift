# License and attributions

## This project

Copyright (c) 2026. Released under the BSD 3-Clause License (below).

## Dependencies and derived material

NuSIFT builds on, and in places derives from, the following. Each remains under its
own license; nothing here relicenses them.

1. **cram-depletion** (BSD 3-Clause) supplies the CRAM solver, the depletion chain, and
   the ENDF decay/fission-yield readers. Consumed as a package or a sub-build.

2. **Eigen** (MPL2) via cram, for sparse linear algebra.

3. **NIST X-Ray Mass Attenuation Coefficients** (Hubbell & Seltzer, NISTIR 5632) — the dry-air
   mass attenuation and mass energy-absorption tables embedded in `nusift/exposure/`. NIST
   data products are not subject to copyright in the United States.

4. **ENDF/B** evaluated nuclear data (IAEA / Brookhaven National Laboratory), read at staging
   time. Nuclear data files are distributed by the evaluators under their own terms; NuSIFT
   ships derived stores, not the tapes.

5. Build-time and test-time dependencies retain their own licenses: **CLI11** (BSD 3-Clause),
   **nlohmann/json** (MIT), **toml++** (MIT), **GoogleTest** (BSD 3-Clause), **nanobind**
   (BSD 3-Clause), **HDF5** (BSD-style), and, for staging builds only, **njoy/ENDFtk**
   (BSD 3-Clause, Copyright (c) njoy contributors).

---

## BSD 3-Clause License

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
   conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of
   conditions and the following disclaimer in the documentation and/or other materials
   provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used
   to endorse or promote products derived from this software without specific prior written
   permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
