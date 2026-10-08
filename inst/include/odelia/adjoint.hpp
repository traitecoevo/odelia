// -*-c++-*-
#ifndef ODELIA_ADJOINT_HPP_
#define ODELIA_ADJOINT_HPP_

#include <algorithm>
#include <cstddef>
#include <memory>
#include <span>
#include <utility>
#include <vector>
#include <XAD/XAD.hpp>
#include <odelia/ode_interface.hpp>
#include <odelia/ode_util.hpp>

// Reverse mode over one recording at a time. The words used here and in
// ode_solver.hpp:
//
//   tape       a log of arithmetic, appended as a calculation runs and walked
//              backwards for derivatives (XAD's).
//   statement  one tape entry, a single left-hand side over a run of
//              operations; cost is counted in statements.
//   recording  what sits on the tape between clears -- here, one step's
//              arithmetic, or one insertion's.
//   seed       the derivative a backward walk starts from, one per output.
//   row        one output's derivatives against a list of inputs; adjoint_rows
//              is a batch of them, because several seeds share one recording.
//   sweep      the backward walk over a run's rows, last to first.
//   range      consecutive rows at one state width; a System of fixed width has
//              one, a System that widens opens a new one at each insertion.
//
// Entry points: vector_jacobian_product for a free function,
// state_and_parameter_adjoints for a map through a System, and
// Solver::solve_adjoint (ode_solver.hpp) for a whole run.

namespace odelia {
namespace ode {

// One adjoint per seed, every row the same width: what a transpose is seeded with,
// and what it hands back.
//
// Several exist together because they share ONE recording. The recording is a model
// evaluation and a sweep is arithmetic, so a caller wanting several rows pays one
// recording rather than one per row -- which is the whole reason this is a set and
// not a vector. One seed is a batch of one, and nothing takes a single seed
// instead: a second signature over the same recording would be a second place for
// the seam between the recording and the sweep to be got wrong.
//
// The shape is the object's own, so a ragged set and a row of the
// wrong width are not things a caller can build, where a vector of vectors builds
// both and every function on the path then has to test for them.
//
// The width is a run-time number on purpose. A sweep narrows as it descends, and
// how many seeds a caller wants is the caller's; neither is a property of a type.
class adjoint_rows {
public:
  adjoint_rows() = default;
  adjoint_rows(std::size_t rows, std::size_t width) { assign(rows, width); }

  std::size_t rows() const { return rows_; }
  std::size_t width() const { return width_; }
  bool empty() const { return rows_ == 0; }

  // Every row zeroed at a new shape. A sweep hands back rows at the width of the
  // range it swept, so a row the new shape does not reach must not survive as
  // a number from the last one.
  void assign(std::size_t rows, std::size_t width) {
    rows_ = rows;
    width_ = width;
    store_.assign(rows * width, 0.0);
  }

  std::span<double> operator[](std::size_t m) {
    return {store_.data() + m * width_, width_};
  }
  std::span<const double> operator[](std::size_t m) const {
    return {store_.data() + m * width_, width_};
  }

  // The rows `which` names, in the order it names them. Out of range is refused
  // rather than clamped: a caller indexing by position would otherwise be handed
  // a different row than the one it asked for.
  adjoint_rows select(const std::vector<std::size_t>& which) const {
    adjoint_rows ret(which.size(), width_);
    for (std::size_t m = 0; m < which.size(); ++m) {
      if (which[m] >= rows_) {
        util::stop("adjoint_rows::select: row " +
                   util::to_string(static_cast<int>(which[m])) +
                   " is outside a batch of " +
                   util::to_string(static_cast<int>(rows_)));
      }
      const std::span<const double> from = (*this)[which[m]];
      std::copy(from.begin(), from.end(), ret[m].begin());
    }
    return ret;
  }

  // A batch of one, which is how a caller wanting a single transpose row asks for
  // it.
  static adjoint_rows one_row(const std::vector<double>& values) {
    adjoint_rows ret(1, values.size());
    std::copy(values.begin(), values.end(), ret[0].begin());
    return ret;
  }

  // The seeds that ask for every row of a Jacobian: one per output, each picking
  // that output out.
  static adjoint_rows all_rows(std::size_t rows) {
    adjoint_rows ret(rows, rows);
    for (std::size_t m = 0; m < rows; ++m) {
      ret[m][m] = 1.0;
    }
    return ret;
  }

  // Rows of doubles for a caller that hands them on -- the R boundary, which
  // takes nothing else. Not how anything inside works.
  std::vector<std::vector<double>> to_rows() const {
    std::vector<std::vector<double>> ret;
    ret.reserve(rows_);
    for (std::size_t m = 0; m < rows_; ++m) {
      const std::span<const double> row = (*this)[m];
      ret.emplace_back(row.begin(), row.end());
    }
    return ret;
  }

private:
  // Set together and nowhere else, so the extent and the shape it is meant to
  // hold cannot disagree.
  std::size_t rows_ = 0, width_ = 0;
  std::vector<double> store_;
};

// The tape running for the duration, on every exit and exceptions included.
//
// Activated here only if nothing else holds it, and deactivated only where
// activated. Activating a tape twice raises, so a caller holding one across
// several recordings cannot be nested inside by anything that activates; and
// deactivating one an outer scope is holding stops recording mid-recording,
// which reads as a missing derivative term and raises nothing.
template <typename Tape>
struct tape_scope {
  explicit tape_scope(Tape& tape) {
    Tape* const running = Tape::getActive();
    if (running == nullptr) {
      tape.activate();
      held_ = &tape;
    } else if (running != &tape) {
      util::stop("a tape is already active and it is not the one being "
                 "recorded on");
    }
  }
  ~tape_scope() {
    if (held_ != nullptr) {
      held_->deactivate();
    }
  }
  tape_scope(const tape_scope&) = delete;
  tape_scope& operator=(const tape_scope&) = delete;

private:
  Tape* held_ = nullptr;
};

// Hand a value's tape slot back, leaving its number. Move-assigning a freshly
// constructed scalar swaps the two slots, so the temporary carries the old one
// away and its destructor releases it, and nothing is recorded against either.
//
// ⚠️ BEFORE A TAPE RESET AND NOT AFTER. Clearing a tape returns its slot counter
// to zero, and handing back a slot the counter no longer knows about takes it
// below zero; where the library is built to reuse freed slots it then issues the
// same number twice. Released first and cleared second, both stay consistent.
template <class S>
void release_slot(S& x) {
  x = S(util::to_passive(x));
}

// A System rebound to the adjoint scalar, with the addresses of the parameters a
// recording writes through. Both come off one object, so a caller cannot pair a
// System with another System's parameter addresses.
//
// Held across the recordings a walk takes on it, where a copy per recording costs
// every allocation the System owns.
//
// ⚠️ WHAT MAKES THAT LEGAL IS `release`, AND NOTHING ELSE. Writing a member does
// not refresh the slot it carries: assignment keeps the slot the target already
// had, so a member that survives a tape clear writes through a number the next
// recording issues to something else. The rows stay finite and stop being the
// model's, and only once two recordings differ in shape -- which is every sweep
// that narrows. So every active value the System holds is released before the
// clear, and the tape is asked whether any was missed.
template <class System>
struct active_system {
  using scalar = active_scalar<double>;
  using system_type = typename rebound_system<System, scalar>::type;

  // Declared first, so the destructor's tape is initialised before the System it
  // releases. Read back rather than handed to a transpose beside the System: the
  // two arriving separately is a pairing a caller can get wrong.
  //
  // A pointer to a tape someone else holds, never a tape of its own:
  // constructing a `Tape` reserves a large block, so a tape per recording would
  // dominate the cost, while reusing one is cheap.
  adjoint_tape<double>* tape_;
  adjoint_tape<double>& tape() const { return *tape_; }

  active_system(const System& passive, adjoint_tape<double>& tape)
      : tape_(&tape), system(passive.template rebind_from<scalar>()),
        parameters(system.ad_parameters()) {
    static_assert(Rebindable<System, scalar>,
                  "a recording is taken on this System at the adjoint scalar; "
                  "it has no rebind_from()");
  }

  // Handed back on the way out, so a walk that builds one of these per width
  // leaves the tape as it found it. Without this, every width but the last stays
  // registered -- nothing is active when the object goes, so nothing hands its
  // slots back, and the tape zero-fills them on every later recording.
  //
  // Skipped where another tape is active: these slots are not on it, and
  // handing them back would decrement its count instead. Nothing reports the
  // skip; the next release() on this tape counts what was left behind.
  ~active_system() {
    adjoint_tape<double>* const running = adjoint_tape<double>::getActive();
    if (running != nullptr && running != tape_) {
      return;
    }
    tape_scope<adjoint_tape<double>> holding{*tape_};
    system.for_each_active([](scalar& x) { release_slot(x); });
  }

  // The parameter addresses point into `system`, so a copy would hand its own
  // parameters to the object it was copied from. Declaring the destructor above
  // suppresses the move operations, which would leave those addresses pointing
  // into the object moved from.
  active_system(const active_system&) = delete;
  active_system& operator=(const active_system&) = delete;

  // Every active value back to unregistered, so the clear below it is safe and
  // the next recording registers them fresh.
  //
  // Releasing every value `for_each_active` reaches has to return the tape to
  // the number of values it counted before this System existed. So a System
  // holding one the walk does not reach is caught here, rather than producing a
  // row that is finite and is not the model's.
  //
  // Against zero, which the destructor above is what makes true: every clear
  // returns the count to zero, every recording registers this System's values,
  // and the temporaries a recording makes are gone by the time it ends.
  void release() {
    system.for_each_active([](scalar& x) { release_slot(x); });
    const std::size_t still = tape_->getNumVariables();
    if (still != 0) {
      util::stop(
          "active_system::release: " +
          util::to_string(static_cast<int>(still)) +
          " active values of this System are still registered after releasing "
          "every one for_each_active reaches, so it holds some the walk does "
          "not. A recording that reads one of those reads a slot the next clear "
          "issues to something else.");
    }
  }

  system_type system;
  std::vector<scalar*> parameters;
};

namespace internal {

// Seed the recorded outputs and sweep, once per seed row; `read(m)` runs with
// seed m's adjoints standing on the tape.
//
// The derivatives are cleared between sweeps. Without it the previous seed's
// adjoints are still on the slots and every row after the first is the running
// sum of the ones before it. Clearing leaves the recorded operations alone, so
// each sweep is bit-identical to the row a fresh recording would give: one
// recording is substitutable for one per seed rather than an approximation of it.
template <class Read>
void sweep_each_seed(adjoint_tape<double>& tape,
                     std::vector<active_scalar<double>>& outputs,
                     const adjoint_rows& output_adjoints, Read&& read) {
    for (std::size_t m = 0; m < output_adjoints.rows(); ++m) {
        tape.clearDerivatives();
        // The adjoint slots, named directly rather than through an accessor of
        // their own: this is the only place in the family that touches one, and a
        // name for two sites inside the function that owns the tape would say
        // nothing the tape does not.
        for (std::size_t i = 0; i < outputs.size(); ++i) {
            xad::derivative(outputs[i]) = output_adjoints[m][i];
        }
        tape.computeAdjoints();
        read(m);
    }
}

}  // namespace internal

// One block of `f`, recorded ONCE on the tape handed in and swept once per seed:
// `input_adjoints[m]` receives transpose(jacobian) * output_adjoints[m], and the return
// value is the recording's size. `f` is instantiated at the active scalar here, so only
// doubles cross in and out.
//
// The batch carries its own shape (adjoint_rows), so a seed of the wrong width is not
// something a caller can hand in and not something this has to test for.
//
// The tape is the caller's and is reused across calls, so nothing here allocates one.
// Stops if a tape other than this one is active: recording onto a tape this product does
// not own would sweep the block's adjoints twice.
//
// An empty seed is swept anyway rather than skipped: the row is then zeros, which is
// what the caller's accumulator expects, and skipping would make the result depend on
// which seeds happen to vanish at this state.
template <class F>
std::size_t vector_jacobian_product(adjoint_tape<double>& tape,
                                     const std::vector<double>& x,
                                     const adjoint_rows& output_adjoints,
                                     F&& f,
                                     adjoint_rows& input_adjoints) {
    using ad_type = active_scalar<double>;
    using tape_type = adjoint_tape<double>;

    if (x.empty()) {
        util::stop("vector_jacobian_product: 'x' must have at least one entry");
    }
    if (output_adjoints.empty()) {
        util::stop("vector_jacobian_product: needs at least one seed");
    }
    const std::size_t n_out = output_adjoints.width();
    if (n_out == 0) {
        util::stop("vector_jacobian_product: a seed must have at least one entry");
    }

    tape_scope<tape_type> running{tape};
    tape.clearAll();

    std::vector<ad_type> x_active(x.begin(), x.end());
    tape.registerInputs(x_active);
    tape.newRecording();

    std::vector<ad_type> y_active(n_out);
    f(x_active, y_active);
    if (y_active.size() != n_out) {
        util::stop("vector_jacobian_product: 'f' resized the output buffer; it is "
                   "handed one entry per output adjoint and must write in place");
    }
    tape.registerOutputs(y_active);

    input_adjoints.assign(output_adjoints.rows(), x.size());
    internal::sweep_each_seed(tape, y_active, output_adjoints,
                              [&](std::size_t m) -> void {
        for (std::size_t i = 0; i < x_active.size(); ++i) {
            input_adjoints[m][i] = xad::derivative(x_active[i]);
        }
    });

    return tape.getMemory();
}

// A transpose taken with respect to a System's state AND the parameters its
// active-scalar lists: one recording over both, swept once per seed.
//
// The active System is the caller's, held across the recordings a walk takes on
// it. It is released here rather than there, because the release has to sit
// immediately before the clear and the clear is inside the product below.
//
// `evaluate` is handed that System, the state half of the recorded inputs, and
// the output buffer, with the System's parameters ALREADY written from the other
// half. That order is why the two halves are one recording rather than two
// calls: a quantity the state determines reads the parameters while deriving it,
// so a state written first derives it at the previous values.
//
// The sweep splits back along the same seam, and the two halves are handled
// differently. `state_adjoint` is REPLACED: it is resized to one row per seed and
// each row assigned, so rows a shorter batch does not reach do not survive as
// stale numbers. `parameter_adjoint` is ADDED to, because a parameter is reached
// once per step and its gradient is the sum over every step swept, so it is the
// caller's to clear once per sweep and may carry rows past the seeds handed in.
//
// They are therefore different objects, and passing one as both is refused:
// that would resize the accumulator between the check on its rows and the
// writes into them.
//
// One width for the whole batch, so a row of parameter adjoints cannot be a
// different length from its neighbours and there is nothing to test per row.
template <class System, class Evaluate>
std::size_t state_and_parameter_adjoints(
    active_system<System>& active, const std::vector<double>& state,
    const adjoint_rows& output_adjoints, Evaluate&& evaluate,
    adjoint_rows& state_adjoint, adjoint_rows& parameter_adjoint) {
    using scalar = active_scalar<double>;
    adjoint_tape<double>& tape = active.tape();
    const std::vector<scalar*>& parameters = active.parameters;
    const std::size_t n_state = state.size();
    const std::size_t n_parameter = parameters.size();
    const std::size_t n_seed = output_adjoints.rows();

    if (&state_adjoint == &parameter_adjoint) {
        util::stop("state_and_parameter_adjoints: the state adjoints are "
                   "replaced and the parameter adjoints accumulated, so they "
                   "cannot be the same batch");
    }
    if (parameter_adjoint.rows() < n_seed) {
        util::stop("state_and_parameter_adjoints: one row of parameter adjoints "
                   "per seed, to accumulate into");
    }
    util::check_length(parameter_adjoint.width(), n_parameter);
    if (state.empty()) {
        util::stop("state_and_parameter_adjoints: 'state' must have at least one "
                   "entry");
    }
    if (output_adjoints.empty()) {
        util::stop("state_and_parameter_adjoints: needs at least one seed");
    }
    const std::size_t n_out = output_adjoints.width();
    if (n_out == 0) {
        util::stop("state_and_parameter_adjoints: a seed must have at least one "
                   "entry");
    }

    // Held for the release, the recording and the sweeps. The release comes
    // first and the clear second (see release_slot); a slot is handed back by the
    // destructor of the temporary that takes it, which does nothing with no tape
    // active.
    tape_scope<adjoint_tape<double>> running{tape};
    active.release();
    tape.clearAll();

    // The two kinds of input, each registered as what it is. The state is this
    // recording's own and arrives as doubles; the parameters live on the System
    // and are registered where they sit. So no value is copied in and no adjoint
    // is split back out by position.
    std::vector<scalar> state_active(state.begin(), state.end());
    tape.registerInputs(state_active);
    for (scalar* p : parameters) {
        tape.registerInput(*p);
    }
    tape.newRecording();

    std::vector<scalar> outputs(n_out);
    evaluate(active.system, state_active.begin(), outputs);
    if (outputs.size() != n_out) {
        util::stop("state_and_parameter_adjoints: the map resized the output "
                   "buffer; it is handed one entry per output adjoint and must "
                   "write in place");
    }
    tape.registerOutputs(outputs);

    state_adjoint.assign(n_seed, n_state);
    internal::sweep_each_seed(tape, outputs, output_adjoints,
                              [&](std::size_t m) -> void {
        for (std::size_t i = 0; i < n_state; ++i) {
            state_adjoint[m][i] = xad::derivative(state_active[i]);
        }
        for (std::size_t p = 0; p < n_parameter; ++p) {
            parameter_adjoint[m][p] += xad::derivative(*parameters[p]);
        }
    });
    return tape.getMemory();
}

// The same transpose where the caller holds a double System rather than an active
// one and is taking a single recording: the tape and the rebind are made here. A
// walk taking many hands in the one it holds for its descent, and rebinds once per
// width.
template <class System, class Evaluate>
    requires Rebindable<System, active_scalar<double>>
std::size_t state_and_parameter_adjoints(
    const System& passive, const std::vector<double>& state,
    const adjoint_rows& output_adjoints, Evaluate&& evaluate,
    adjoint_rows& state_adjoint, adjoint_rows& parameter_adjoint) {
    adjoint_tape<double> tape(false);
    active_system<System> active{passive, tape};
    return state_and_parameter_adjoints(active, state, output_adjoints,
                                        std::forward<Evaluate>(evaluate),
                                        state_adjoint, parameter_adjoint);
}

// The transpose of one rate evaluation: `state_adjoint[m]` receives
// transpose(d dydt / d y) * rate_adjoints[m], with the System's parameters in the
// same recording, so a rate the parameters reach carries their rows too.
//
// What is recorded is derivs() -- the call the forward pass makes -- so the
// transpose cannot drift from the rates it transposes, and a System that
// restores a recorded field restores it here as well. Everything between the
// state and the rates is an intermediate of that one recording, so nothing
// between them needs a transpose written for it.
//
// The rate is taken at the time, not at a recorded stage. The System recorded on
// is rebound for this recording and holds no recorded field, so there is no stage
// for it to read back -- a stage argument here would name a capability the active
// copy cannot have.
template <class System>
std::size_t rates_adjoint(
    active_system<System>& active,
    const std::vector<double>& state, double time,
    const adjoint_rows& rate_adjoints, adjoint_rows& state_adjoint,
    adjoint_rows& parameter_adjoint) {
    using scalar = active_scalar<double>;
    const std::size_t n_state = state.size();
    auto rates = [&](auto& sys,
                     typename std::vector<scalar>::const_iterator x,
                     std::vector<scalar>& dydt) -> void {
        // The recorded state inputs. The parameters are registered where they
        // sit on the System, so they are not in this buffer.
        std::vector<scalar> y(x, x + static_cast<std::ptrdiff_t>(n_state));
        ode::derivs(sys, y, dydt, time);
    };
    return state_and_parameter_adjoints(active, state, rate_adjoints, rates,
                                        state_adjoint, parameter_adjoint);
}

}
}

#endif
