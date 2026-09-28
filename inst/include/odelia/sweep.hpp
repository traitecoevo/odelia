// -*-c++-*-
#ifndef ODELIA_SWEEP_HPP_
#define ODELIA_SWEEP_HPP_

#include <odelia/adjoint.hpp>
#include <odelia/ode_solver.hpp>
#include <algorithm>
#include <limits>
#include <span>
#include <utility>
#include <vector>

// Carrying an adjoint back over a recording whose state width changed. One tape,
// held active for the whole descent; one System rebound per width, releasing its
// slots before each recording clears the tape.

namespace odelia {
namespace ode {

// The state a range's first step ran from, and the time it ran at: the insertion
// row's own. It is replayed rather than read off that row because applying the map
// is what leaves the System at the range's width, and a caller reading this state
// wants to run from it. What the row holds is then what the replay is checked
// against.
template <class System>
double state_at_range(System& system,
                      std::span<const step_record<System>> rec,
                      std::size_t range,
                      std::vector<double>& base, std::size_t& start) {
    // The insertion rows, in order, off the kind the run recorded and not off a
    // width that grew: an inference has to refuse a width that shrinks, where a
    // recorded fact cannot be wrong. Row 0 is never one -- an insertion runs on the
    // state the row below it holds.
    std::vector<std::size_t> rows;
    for (std::size_t k = 1; k < rec.size(); ++k) {
        if (rec[k].insertion) {
            rows.push_back(k);
        }
    }
    if (range > rows.size()) {
        util::stop("state_at_range: the recording has no such range");
    }
    if (range == 0) {
        be_at_step(system, rec, 0);
        start = 0;
        base.assign(system.ode_size(), 0.0);
        system.ode_state(base.begin());
        return rec[0].time;
    }
    // The insertion's own map at the row below it, which is the same map the sweep
    // transposes and the same one the run recorded the result of -- so the width is
    // checked against the row the run wrote it into.
    start = rows[range - 1];
    be_at_step(system, rec, start - 1);
    apply_insertion(system, rec[start].time, rec[start - 1].state.begin(), base);
    util::check_length(base.size(), rec[start].state.size());
    return rec[start].time;
}

// Step `forward` over the recording from `from_range` on, at the sizes the run
// took, widening where the run widened. `first` is the recorded row the walk
// begins at, which the caller has already put the System on. The schedule is read
// off the recording, because a recording row is a schedule row plus its state.
//
// Returns the recording's suffix behind a head instruction the caller owns. The
// head's time is the solver's own rather than the row's, which is what lets a
// caller replay a perturbed state at that time. It replaces row `first`, so a
// walk resuming after an insertion names that insertion's row and a walk starting
// at the beginning names row 0, and neither has to say which of them still owes a
// state map.
//
// The step sizes are replayed rather than the times: a size differenced back out
// of two recorded times is not the size that was taken, since fl(fl(t + h) - t)
// is not h, and a walk that chose its own would be differentiating a controller
// the model does not contain.
template <class Record>
std::vector<instruction> program_from(std::span<const Record> rec,
                                      std::size_t first, instruction head) {
    if (first >= rec.size()) {
        util::stop("program_from: the recording has no such row");
    }
    std::vector<instruction> ret;
    ret.reserve(rec.size() - first);
    ret.push_back(head);
    for (std::size_t k = first + 1; k < rec.size(); ++k) {
        ret.push_back(rec[k]);
    }
    return ret;
}

}  // namespace ode
}  // namespace odelia

#endif
