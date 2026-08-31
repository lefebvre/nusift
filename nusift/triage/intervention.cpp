#include "nusift/triage/intervention.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "nusift/core/element_symbols.hpp"
#include "nusift/core/error.hpp"
#include "nusift/core/nuclide.hpp"
#include "nusift/core/nuclide_name.hpp"
#include "nusift/engine/adjoint_engine.hpp"
#include "nusift/nucdata/nuclear_data.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "intervention";

std::string_view trimmed(std::string_view text) {
  const std::size_t first = text.find_first_not_of(" \t");
  if (first == std::string_view::npos) {
    return {};
  }
  return text.substr(first, text.find_last_not_of(" \t") - first + 1);
}

// What a selector matches. Held as a resolved predicate rather than re-parsed per nuclide, so
// the spelling is validated once and a typo cannot be diagnosed differently on two rows.
struct Selection {
  enum class Kind { Nuclide, Element, MassChain } kind = Kind::Nuclide;
  std::int64_t value = 0;  // ZAI key, Z, or A

  bool matches(const Zai& zai) const {
    switch (kind) {
      case Kind::Nuclide:
        return zai.key() == value;
      case Kind::Element:
        return zai.z == static_cast<int>(value);
      case Kind::MassChain:
        return zai.a == static_cast<int>(value);
    }
    return false;
  }
};

// Strip a "z=" or "a=" qualifier, reporting whether it was there. The qualified forms exist
// because a bare number is ambiguous between an atomic number and a mass number, and guessing
// which one a user meant would silently remove the wrong material.
bool stripQualifier(std::string_view& text, char letter) {
  if (text.size() >= 2 && (text[0] == letter || text[0] == static_cast<char>(letter - 32)) &&
      text[1] == '=') {
    text = trimmed(text.substr(2));
    return true;
  }
  return false;
}

std::optional<int> wholeNumber(std::string_view text) {
  if (text.empty() ||
      !std::all_of(text.begin(), text.end(), [](char c) { return c >= '0' && c <= '9'; })) {
    return std::nullopt;
  }
  int value = 0;
  for (const char c : text) {
    value = value * 10 + (c - '0');
    if (value > 100000) {
      return std::nullopt;
    }
  }
  return value;
}

Selection resolveSelector(std::string_view raw) {
  const std::string_view whole = trimmed(raw);
  if (whole.empty()) {
    throw InputError(tagged(kModule, "a removal needs something to remove"));
  }

  std::string_view rest = whole;
  if (stripQualifier(rest, 'a')) {
    if (const std::optional<int> number = wholeNumber(rest);
        number.has_value() && *number >= 1 && *number <= 400) {
      return Selection{Selection::Kind::MassChain, *number};
    }
    throw InputError(tagged(kModule, "\"" + std::string(whole) + "\" is not a mass number"));
  }
  rest = whole;
  if (stripQualifier(rest, 'z')) {
    if (const std::optional<int> number = wholeNumber(rest);
        number.has_value() && *number >= 1 && *number <= kMaxAtomicNumber) {
      return Selection{Selection::Kind::Element, *number};
    }
    throw InputError(tagged(kModule, "\"" + std::string(whole) + "\" is not an atomic number"));
  }

  // A nuclide name before a bare element symbol: "Cs-137" and "Cs" are both valid and mean
  // different amounts of material, so the more specific reading wins where both could parse.
  if (const std::optional<Zai> zai = parseNuclideName(whole); zai.has_value()) {
    return Selection{Selection::Kind::Nuclide, zai->key()};
  }
  if (const int z = atomicNumber(whole); z > 0) {
    return Selection{Selection::Kind::Element, z};
  }

  throw InputError(tagged(kModule, "\"" + std::string(whole) +
                                       "\" is not a nuclide, an element symbol, Z=55, or A=137"));
}

void requireUsableTimes(double interventionTime, double responseTime) {
  if (!std::isfinite(interventionTime) || interventionTime < 0.0) {
    throw InputError(tagged(kModule, "the intervention time must be finite and not negative"));
  }
  if (!std::isfinite(responseTime)) {
    throw InputError(tagged(kModule, "the response time must be finite"));
  }
  if (responseTime < interventionTime) {
    // Removing something cannot change what a response already was. Asked backwards this would
    // silently answer a question about the wrong window.
    throw InputError(tagged(kModule,
                            "the response time is before the intervention time; removing "
                            "something on a later date cannot change an earlier response"));
  }
}

// The state at the intervention date: the inventory the adjoint is seeded from, and -- kept
// separately -- every nuclide the original seed's chain reaches.
//
// The two differ, and the difference is load-bearing. forwardClosure() roots itself only at
// non-zero entries, so a nuclide that has decayed to nothing by t0 is pruned out of the adjoint's
// index space; that is exact, since it contributes nothing. But it also makes "the tin is all
// gone by then" indistinguishable from "this inventory has no tin in its chain", and only the
// second is a mistake the user made. Selectors are therefore checked against `reachable`, while
// the benefit is computed over the adjoint's smaller space.
struct StateAtIntervention {
  Inventory inventory;
  std::vector<std::int64_t> reachable;
};

StateAtIntervention stateAt(const NuclearData& data, const Inventory& inventory, double time,
                            const DecayOptions& options) {
  // Run even for t0 = 0, where it returns the seed itself: the closure it reports is the set a
  // selector has to be judged against, and special-casing zero would give that case a different
  // and smaller answer.
  const DecayResult result = decay(data, inventory, std::vector<double>{time}, options);

  StateAtIntervention state;
  state.reachable = result.nuclideKeys;
  const std::span<const double> atoms = result.atomsAt(0);
  for (int i = 0; i < result.nuclideCount(); ++i) {
    // Clamped at zero: CRAM's rational approximation leaves residuals of either sign around a
    // count that has decayed to nothing -- a few atoms either way out of 1e20. An atom count is
    // non-negative by definition and Inventory enforces that, so the noise is resolved here, and
    // resolved DOWNWARDS. Taking the magnitude would invent material that is not there.
    const double count = std::max(0.0, atoms[static_cast<std::size_t>(i)]);
    if (count > 0.0) {
      state.inventory.addKey(result.nuclideKeys[static_cast<std::size_t>(i)], count);
    }
  }
  state.inventory.setProvenance(inventory.provenance());
  return state;
}

}  // namespace

InterventionStudy compareInterventions(const NuclearData& data, const Inventory& inventory,
                                       double interventionTime, double responseTime,
                                       const ResponseSpec& spec,
                                       std::span<const Intervention> interventions,
                                       const DecayOptions& options) {
  requireUsableTimes(interventionTime, responseTime);
  if (interventions.empty()) {
    throw InputError(tagged(kModule, "no interventions given, so there is nothing to compare"));
  }
  if (spec.aggregate == Aggregate::GammaLine) {
    // A line is a way an emitter's decays get out, not a thing an inventory holds, so there is
    // nothing to take out of one.
    throw InputError(tagged(kModule,
                            "a photon line cannot be removed from an inventory -- rank by "
                            "gamma line, but intervene on nuclides"));
  }

  // Resolve and validate every selector before solving anything, so a typo in the last
  // intervention does not surface only after the expensive part.
  std::vector<std::vector<Selection>> selections(interventions.size());
  for (std::size_t q = 0; q < interventions.size(); ++q) {
    if (interventions[q].name.empty()) {
      throw InputError(tagged(kModule, "intervention " + std::to_string(q) +
                                           " has no name; naming the one that wins is the answer"));
    }
    if (interventions[q].removals.empty()) {
      throw InputError(tagged(kModule, "\"" + interventions[q].name + "\" removes nothing"));
    }
    for (const Removal& removal : interventions[q].removals) {
      if (!(removal.fraction >= 0.0 && removal.fraction <= 1.0)) {
        throw InputError(tagged(kModule, "\"" + interventions[q].name + "\" removes a fraction " +
                                             std::to_string(removal.fraction) + " of " +
                                             removal.selector +
                                             "; a removed fraction lies between 0 and 1"));
      }
      selections[q].push_back(resolveSelector(removal.selector));
    }
  }

  const std::vector<double> weight = responseWeights(data, spec);
  const StateAtIntervention state = stateAt(data, inventory, interventionTime, options);
  // The single adjoint solve. Its importance vector is dR(T)/dn_i(t0) for every nuclide still
  // present at t0, and every intervention below is a dot product against it.
  const SeedImportance importance =
      seedImportance(data, state.inventory, weight, responseTime - interventionTime, options);

  // Where each surviving nuclide sits in the adjoint's space. A key the chain reaches but that
  // is absent here decayed to nothing by t0: removing it is worth exactly zero, which is an
  // answer rather than a failure.
  std::map<std::int64_t, int> column;
  for (int i = 0; i < static_cast<int>(importance.nuclideKeys.size()); ++i) {
    column.emplace(importance.nuclideKeys[static_cast<std::size_t>(i)], i);
  }

  InterventionStudy study;
  study.metric = spec.metric;
  study.unit = spec.unit;
  study.interventionTimeSeconds = interventionTime;
  study.responseTimeSeconds = responseTime;
  study.baseline = importance.response;
  study.seedProvenance = inventory.provenance();

  for (std::size_t q = 0; q < interventions.size(); ++q) {
    // The fraction removed of each nuclide, accumulated across this intervention's removals. A
    // nuclide named twice is refused rather than resolved: half of something taken out twice is
    // not a stated quantity, and picking 75% or 100% for the user would be a guess.
    std::map<std::int64_t, double> removedFraction;
    for (std::size_t r = 0; r < selections[q].size(); ++r) {
      const Selection& selection = selections[q][r];
      bool reached = false;
      // Matched against the whole chain rather than the adjoint's space, so that a nuclide the
      // chain reaches but that is gone by t0 still counts as named.
      for (const std::int64_t key : state.reachable) {
        if (!selection.matches(Zai::fromKey(key))) {
          continue;
        }
        reached = true;
        const auto [entry, inserted] =
            removedFraction.try_emplace(key, interventions[q].removals[r].fraction);
        if (!inserted) {
          throw InputError(tagged(
              kModule, "\"" + interventions[q].name + "\" removes " +
                           formatNuclideName(Zai::fromKey(key)) + " twice, through \"" +
                           interventions[q].removals[r].selector +
                           "\" and an earlier removal. Removing part of something twice is not "
                           "a stated quantity"));
        }
      }
      if (!reached) {
        // Refused for the reason a pin naming nothing is refused: a removal that silently
        // matches no nuclide reads as "taking this out is worth nothing", when in fact the
        // question never reached the chain.
        throw InputError(tagged(kModule, "\"" + interventions[q].removals[r].selector +
                                             "\" names nothing this inventory's chain reaches"));
      }
    }

    InterventionEffect effect;
    effect.name = interventions[q].name;
    for (const auto& [key, fraction] : removedFraction) {
      const auto found = column.find(key);
      if (found == column.end()) {
        // In the chain, but decayed to nothing by t0 and pruned from the adjoint. Worth exactly
        // zero to remove, which is a real answer and not a row worth printing.
        continue;
      }
      const std::size_t i = static_cast<std::size_t>(found->second);
      const double atoms = importance.seedAtoms[i] * fraction;
      const double value = atoms * importance.importance[i];
      if (!(atoms > 0.0)) {
        continue;
      }
      effect.removed += value;
      effect.contributors.push_back(RemovedContributor{
          importance.nuclideKeys[i], formatNuclideName(Zai::fromKey(importance.nuclideKeys[i])),
          atoms, value, 0.0});
    }

    effect.response = study.baseline - effect.removed;
    effect.removedFraction = study.baseline > 0.0 ? effect.removed / study.baseline : 0.0;
    for (RemovedContributor& contributor : effect.contributors) {
      contributor.fraction = effect.removed > 0.0 ? contributor.value / effect.removed : 0.0;
    }
    // By what each nuclide's removal was worth, which is the order the answer is read in --
    // ties broken on key so two equally worthless removals do not swap between runs.
    std::sort(effect.contributors.begin(), effect.contributors.end(),
              [](const RemovedContributor& a, const RemovedContributor& b) {
                if (a.value != b.value) {
                  return a.value > b.value;
                }
                return a.key < b.key;
              });
    study.effects.push_back(std::move(effect));
  }
  return study;
}

}  // namespace nusift
