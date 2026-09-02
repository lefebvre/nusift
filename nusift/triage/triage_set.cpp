#include "nusift/triage/triage_set.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "nusift/core/error.hpp"

namespace nusift {
namespace {

constexpr const char* kModule = "triage set";

[[noreturn]] void fail(const std::string& what) {
  throw InputError(tagged(kModule, what));
}

// What makes two columns in two tables the same contributor. The key alone is not enough for a
// gamma-line table, where one emitter holds a column per line and covering its 662 keV line says
// nothing about its 32 keV one.
struct Identity {
  std::int64_t key = 0;
  double lineEnergyEv = 0.0;

  bool operator<(const Identity& other) const {
    if (key != other.key) {
      return key < other.key;
    }
    return lineEnergyEv < other.lineEnergyEv;
  }
};

Identity identityOf(const ContributorId& id) {
  return Identity{id.key, id.lineEnergyEv};
}

// One (requirement, time) pair: the floor to clear and how much of it is cleared so far.
struct Constraint {
  int requirement = 0;
  int timeIndex = 0;
  double total = 0.0;
  double required = 0.0;  // fraction * total, in the table's unit
  double achieved = 0.0;
  bool met = false;
};

void validate(std::span<const CoverageRequirement> requirements) {
  if (requirements.empty()) {
    fail("no requirements: a set has to be good FOR something");
  }
  for (const CoverageRequirement& r : requirements) {
    if (r.table == nullptr) {
      fail("a requirement carries no response table");
    }
    if (r.label.empty()) {
      fail(
          "every requirement needs a name; a shortfall reported against requirement 2 of 4 names "
          "nothing a reader can act on");
    }
    if (!(r.fraction > 0.0) || r.fraction > 1.0) {
      fail("a coverage fraction is in (0, 1], got " + std::to_string(r.fraction));
    }
    // A nuclide column and a gamma-line column are not the same kind of thing, and a set mixing
    // them would answer neither question. Refused rather than reconciled: which of the two the
    // caller meant is not something this code can know.
    if (r.table->aggregate != requirements.front().table->aggregate) {
      fail("requirement \"" + r.label + "\" ranks by " +
           std::string(aggregateName(r.table->aggregate)) + " where \"" +
           requirements.front().label + "\" ranks by " +
           std::string(aggregateName(requirements.front().table->aggregate)) +
           "; one set cannot cover both");
    }
  }
}

}  // namespace

TriageSet robustTriageSet(std::span<const CoverageRequirement> requirements) {
  validate(requirements);

  // --- the shared contributor space ------------------------------------------
  //
  // Built from every table at once, because a contributor absent from one is not absent from the
  // problem: it contributes zero there and may be the whole answer somewhere else.
  std::map<Identity, int> indexOf;
  std::vector<ContributorId> ids;
  std::vector<std::string> labels;
  std::vector<std::vector<int>> globalOf(requirements.size());

  for (std::size_t r = 0; r < requirements.size(); ++r) {
    const ResponseTable& table = *requirements[r].table;
    globalOf[r].resize(static_cast<std::size_t>(table.contributorCount()));
    for (int c = 0; c < table.contributorCount(); ++c) {
      const ContributorId& id = table.contributors[static_cast<std::size_t>(c)];
      const auto [it, inserted] = indexOf.emplace(identityOf(id), static_cast<int>(ids.size()));
      if (inserted) {
        ids.push_back(id);
        labels.push_back(table.labels[static_cast<std::size_t>(c)]);
      }
      globalOf[r][static_cast<std::size_t>(c)] = it->second;
    }
  }
  const std::size_t nC = ids.size();

  // --- the constraints -------------------------------------------------------
  std::vector<Constraint> constraints;
  for (std::size_t r = 0; r < requirements.size(); ++r) {
    const ResponseTable& table = *requirements[r].table;
    for (int k = 0; k < table.timeCount(); ++k) {
      Constraint constraint;
      constraint.requirement = static_cast<int>(r);
      constraint.timeIndex = k;
      constraint.total = table.totals[static_cast<std::size_t>(k)];
      constraint.required = requirements[r].fraction * constraint.total;
      // Nothing to cover is covered. A time at which a metric is identically zero -- an exposure
      // before any photon emitter exists -- constrains no set, and treating it as unmeetable
      // would make every answer report a shortfall it does not have.
      constraint.met = !(constraint.required > 0.0);
      constraints.push_back(constraint);
    }
  }

  TriageSet result;
  result.candidateCount = static_cast<int>(nC);
  result.constraintCount = static_cast<int>(constraints.size());

  // --- the greedy ------------------------------------------------------------
  std::vector<bool> chosen(nC, false);
  std::vector<double> gain(nC, 0.0);

  while (true) {
    const bool anyUnmet = std::any_of(constraints.begin(), constraints.end(),
                                      [](const Constraint& c) { return !c.met; });
    if (!anyUnmet) {
      break;
    }

    std::fill(gain.begin(), gain.end(), 0.0);
    for (const Constraint& constraint : constraints) {
      if (constraint.met) {
        continue;
      }
      const ResponseTable& table =
          *requirements[static_cast<std::size_t>(constraint.requirement)].table;
      const std::span<const double> values = table.valuesAt(constraint.timeIndex);
      const double remaining = constraint.required - constraint.achieved;
      const std::vector<int>& map = globalOf[static_cast<std::size_t>(constraint.requirement)];
      for (std::size_t c = 0; c < values.size(); ++c) {
        const int global = map[c];
        if (chosen[static_cast<std::size_t>(global)]) {
          continue;
        }
        // Credit is capped at the shortfall it actually closes. Without the cap a contributor
        // that overwhelms one already-nearly-met constraint would outrank one that is the only
        // way to move three others, which is the failure mode plain greedy has here.
        gain[static_cast<std::size_t>(global)] += std::min(values[c], remaining);
      }
    }

    int best = -1;
    double bestGain = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
      if (chosen[c] || !(gain[c] > 0.0)) {
        continue;
      }
      // Ties broken by contributor key, then line energy, so the same tables give the same list
      // on any machine. A monitoring list that depended on iteration order would be a different
      // list every time it was regenerated.
      if (best < 0 || gain[c] > bestGain ||
          (gain[c] == bestGain &&
           identityOf(ids[c]) < identityOf(ids[static_cast<std::size_t>(best)]))) {
        best = static_cast<int>(c);
        bestGain = gain[c];
      }
    }
    if (best < 0) {
      // Nothing left can improve any unmet constraint. The floor is above what these tables can
      // express -- a gamma-line total counts lines below the column floor, so even every column
      // together falls short of it -- and the shortfalls below say by how much.
      break;
    }

    SetMember member;
    member.id = ids[static_cast<std::size_t>(best)];
    member.label = labels[static_cast<std::size_t>(best)];
    member.order = static_cast<int>(result.members.size()) + 1;
    member.closedShortfall = bestGain;
    result.members.push_back(std::move(member));
    chosen[static_cast<std::size_t>(best)] = true;

    for (Constraint& constraint : constraints) {
      if (constraint.met) {
        continue;
      }
      const ResponseTable& table =
          *requirements[static_cast<std::size_t>(constraint.requirement)].table;
      const std::vector<int>& map = globalOf[static_cast<std::size_t>(constraint.requirement)];
      const std::span<const double> values = table.valuesAt(constraint.timeIndex);
      for (std::size_t c = 0; c < values.size(); ++c) {
        if (map[c] == best) {
          constraint.achieved += values[c];
          break;
        }
      }
      constraint.met = constraint.achieved >= constraint.required;
    }
  }

  // --- what the set actually achieves ----------------------------------------
  //
  // Recomputed against the full totals rather than read off the greedy's bookkeeping. The two
  // agree, and that is the point: a coverage figure derived from the same running sum that chose
  // the set could not catch an error in either.
  double bestMargin = std::numeric_limits<double>::infinity();
  for (const Constraint& constraint : constraints) {
    const CoverageRequirement& requirement =
        requirements[static_cast<std::size_t>(constraint.requirement)];
    const ResponseTable& table = *requirement.table;
    const std::span<const double> values = table.valuesAt(constraint.timeIndex);
    const std::vector<int>& map = globalOf[static_cast<std::size_t>(constraint.requirement)];

    double held = 0.0;
    for (std::size_t c = 0; c < values.size(); ++c) {
      if (chosen[static_cast<std::size_t>(map[c])]) {
        held += values[c];
      }
    }

    CoveragePoint point;
    point.requirement = requirement.label;
    point.timeSeconds = table.times[static_cast<std::size_t>(constraint.timeIndex)];
    point.required = requirement.fraction;
    point.achieved = constraint.total > 0.0 ? held / constraint.total : 1.0;

    if (point.achieved + 1.0e-12 < point.required) {
      result.shortfalls.push_back(point);
    }
    const double margin = point.achieved - point.required;
    if (margin < bestMargin) {
      bestMargin = margin;
      result.binding = point;
    }
  }

  // --- what each member is worth ---------------------------------------------
  for (SetMember& member : result.members) {
    const int global = indexOf.at(identityOf(member.id));
    for (std::size_t r = 0; r < requirements.size(); ++r) {
      const ResponseTable& table = *requirements[r].table;
      for (int c = 0; c < table.contributorCount(); ++c) {
        if (globalOf[r][static_cast<std::size_t>(c)] != global) {
          continue;
        }
        for (int k = 0; k < table.timeCount(); ++k) {
          const double total = table.totals[static_cast<std::size_t>(k)];
          if (total > 0.0) {
            member.peakFraction = std::max(member.peakFraction,
                                           table.valuesAt(k)[static_cast<std::size_t>(c)] / total);
          }
        }
      }
    }
  }

  return result;
}

}  // namespace nusift
