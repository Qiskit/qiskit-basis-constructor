# `BasisConstructor` plugin for Qiskit transpiler

This package provides a translation-stage plugin for the Qiskit transpiler with the key `constructor-beta`.

This translation plugin is an alternative to the Qiskit standard `translator`, which is an amalgamation of Qiskit's built-in `BasisTranslator` and `GateDirection` passes.
This plugin is an exposure of a single transpiler pass, `BasisConstructor`, which does both jobs.

> [!NOTE]
>
> This package is provided as a public beta version of a plugin that is expected to be rewritten and included in Qiskit proper in the future.
> This repository and the associated PyPI package will be archived and no further development on it will be done once that happens.

## Installation

This is available on PyPI under the name `qiskit-basis-constructor`.

```bash
pip install qiskit-basis-constructor
```

You can also install a development version in editable mode directly from this repository:

```bash
git clone https://github.com/Qiskit/qiskit-basis-constructor
pip install -e ./qiskit-basis-constructor
```

## Usage

For basic uses, there is no need to import the package.
Once it is installed, pass `translation_method="constructor-beta"` to either the `transpile` or `generate_preset_pass_manager` functions of Qiskit:

```python
from qiskit import transpile, QuantumCircuit
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_ibm_runtime.fake_provider import FakeTorino

backend = FakeTorino()
qc = QuantumCircuit(2)
qc.h(0)
qc.cx(0, 1)

# Using `generate_preset_pass_manager`
pm = generate_preset_pass_manager(
    optimization_level=2,
    backend=backend,
    translation_method="constructor-beta",
)
pm.run(qc)

# Using `transpile`
transpile(qc, backend, optimization_level=2, translation_method="constructor-beta")
```

(Note the default `translation_method` is `"translator"`, corresponding to the `BasisTranslator`.)

At optimisation levels 0, the scoring heuristic for translations is `(num_2q_gates, num_gates)`.
At optimisation levels 1, 2 (the default for Qiskit 1.3+) and 3, the scoring heuristic is
`(neg_log_fidelity, num_2q_gates, num_gates)`, where the fidelity homogenisation binning has a width of $\ln(10)$ at O2 and O3, while O1 uses the average over all qubit pairs.
See "Score components" below for more detail on these.


### Advanced use

You can more completely configure the pass, including modifying the scoring components, if you
import it from `qiskit_basis_constructor`.
This is tricky to do correctly in an entire Qiskit transpiler pipeline, though, because the
"translation" stage is included by `generate_preset_pass_manager` twice:

1. as the `.translation` stage of the resulting `StagedPassManager`, which is easy to replace (see
   `qiskit_basis_constructor.plugin` for what you should put in there)
2. in the middle of the `optimization` stage to correct non-ISA-aware optimisation passes.
   This is near entirely opaque, and you can't really configure this without remaking the entire optimisation stage.
   (Improving the API here is something that would need changes to Qiskit.)

The easiest way to see the configuration options is just to look at how the plugin is implemented,
in `qiskit_basis_constructor.plugin`.


### Debug use

The pass outputs debugging information in the standard `logging` format at a level of `logging.DEBUG` on loggers descended from `qiskit_basis_constructor`.
To configure basic logging, for example, do:

```python
import logging

logging.basicConfig()
logging.getLogger("qiskit_basis_constructor").setLevel(logging.DEBUG)
```

## Feedback

If you encounter either problems or significant improvements over Qiskit's default translation, please let me know.
You can either use the GitHub Issues here, or message me on Slack (@jakelishman).

### Expected uses

This pass is likely to be better at fractional-gate translations than `BasisTranslator`.
This is primarily because its default equivalence library is more well joined up than the one in
Qiskit that `BasisTranslator` uses.
However, improving `BasisTranslator` is not as simple as switching it to use the new library; the details of the two algorithms (see the "Comparison to `BasisTranslator`" section for more) mean that `BasisTranslator` as a rule of thumb needs the optimal translations to be provided for it directly in the library to be sure it will find them, whereas `BasisConstructor` will _construct_ the optimal translation based on a configurable definition of "optimal".

This pass is also likely to be better at dealing with heterogeneous ISAs (different qubits/links have different operations defined), and targets with 2q links that have more than one defined 2q operation.
`BasisTranslator` is somewhat at the whims of the shallow depths of its `EquivalenceLibrary`, whereas `BasisConstructor` is agnostic to the "depth" of a translation, if it can prove that it is the "best" translation.

Finally, `BasisConstructor` actually does the job of `BasisTranslator` _and_ `GateDirection` in one go, and the concept of directed gates is baked into the construction routines.
You can supply equivalence rules for a 2q gate in terms of itself with permuted qubits (e.g. `cz(0, 1) == cz(1, 0)`) and the construction algorithm output correctly directed gates, and in fact, will raise an exception if it _cannot_, whereas the `BasisTranslator` will happily output gates in an invalid direction, hoping another pass will clean them up.
For historical reasons in Qiskit right now, `GateDirection` is not considered part of the `translation` stage of a preset pass manager, but is an extra pass that's run afterwards.
This means that, even when using the `constructor-beta` translation method, you'll still see `GateDirection` run, it just should have no effect each time.
In the future, we'll move `GateDirection` into the translation plugins for plugins that require it (like `BasisTranslator`), so those that don't (like `BasisConstructor`) don't pay an unnecessary price.

### Known limitations

It's currently totally expected that this pass will take longer than `BasisTranslator`; this pass is still written in pure Python, and has not been optimised.
Algorithmically, there's nothing fundamental that should make this slower, and in fact, it's actually expected to amortise faster than `BasisTranslator` (no matter how complex the heuristics and `Target` are) if the same `PassManager` is used for many circuits.
`BasisConstructor`'s translations are derived purely from the `Target` without reference to an input DAG, unlike `BasisTranslator`, so the cost constructing the translations is a one-off.
In practice, many users of Qiskit will use `transpile` on single circuits, which means we pay the full construction cost each time.

This pass doesn't currently support `BackendV1`-style loose constraints (i.e. passing `basis_gates` and `coupling_map` and not a `Target`).
If you want this, in the interrim, you can probably make do by using `Target.from_configuration` to build a `Target`, then use that.

The message on a failed translation is not currently particularly enlightening.
It will not be difficult to improve this in the future.

## Repo commands

I was using this project as an excuse to try out `hatch` as a Python project manager.
You'll need to install `hatch` to build, run the linters and the tests.

### Build

```
hatch build
```

This will put an sdist and a wheel in the `dist` directory.
The output of the command is pretty self explanatory.

### Dev environment

```
hatch shell
```

This drops you into a new shell with a virtual environment active that has the project installed.
Run `exit` as a command to leave the shell and return to where you were before.

### Lint and style

```
hatch fmt [--check]
```

This runs the linter and the formatter, including any automatic fixes.
Personally, I don't like `hatch`'s default behaviour here of attempting automatic fixes - I'd rather opt into that.
Commands I actually use:

* `hatch fmt -f` or `hatch fmt --formatter` runs only the formatter (in fix-up mode).
* `hatch fmt --check [--linter]` runs the lint and format checks, without the automatic fix-ups.
  I typically autoformat frequently, so `hatch fmt --check` is what I used most of the time, but if
  you only autoformat infrequently, and want to see lint errors in the interrim, you might want
  `hatch fmt --check --linter`.

### Tests

```
hatch test
```

The test suite is a very basic `unittest`-based suite.
This is in anticipation of merging the package into Qiskit proper at some point in the future; it would be a large nuisance to have to convert a `pytest`-based suite back to a `unittest` one, when it's easy enough just to write it in the lowest common denominator to start with.

## Algorithm

At a high level, the constructor algorithm is similar to a Dijkstra search:

1. start from "no gates available"
2. activate all rules whose requirements are fulfilled, and calculate their scores
3. add the gate constructed by the rule with the lowest score to the set of available gates
4. discard any other rules that construct the same gate from consideration
5. repeat steps 2 to 4 until all necessary gates have been constructed (success), or all rules have been exhausted (failure)

Depending on the exact implementation, step 4 can be either proactive (scan through the current candidates and remove all now-invalid rules) or reactive (skip invalid rules as they are encountered).

A "rule" is either:

- an entry in the `Target` for a gate, in which case it has no requirements
- an "equivalence" in the source `EquivalenceLibrary`, in which case its requirements are all the gates that are used in the rule construction

The score for a rule that comes directly from the `Target` is directly derived.
The currently implemented scores for all other rules are found by summing the scores of all constituent components, though one can imagine different reductions being useful or possible.
Some examples of what can a score may comprise:

- the error rates of the gates used
- the number of gates used
- the number of non-Clifford gates used

Heterogeneity can be handled by taking a "necessary" gate to be a pair of the gate itself, and the physical qubits it acts on.
The most straightforward way to do this is to simply run the entire construction algorithm on each set of qargs that are used by the source interaction graph, taking care that whenever a multi-qubit gate is encountered, every subset of the qubits has already had its available basis constructed.
Effectively, there is just a partial order on the qargs that must be respected; `(0,)` and `(1,)` must be constructed before `(0, 1)`, but `(2,)` can be constructed at any point.

Gate direction is handled in a similar way, but some additional care is needed for proper efficiency because a "direction-flipping" equivalence rule (such as `h(0).h(1).cx(0,1).h(0).h(1)` for `cx(1,0)`) is typically mutually recursive with its flipped-qargs version.
The strategy here is to handle each permutation of a multi-qubit set of qargs simultaneously; this may need more tricks to scale to natively-supported massively-multi-qubit gates, but when targeting current hardware that is predominantly limited to two-qubit gates, nothing is yet implemented.


### Score components

Currently, two score components implemented are:

- `GateCount(min_num_qubits=0)`: count the total number of gates in the rule that use at least
  `min_num_qubits` qubits.
- `LogFidelity(bin=None)`: the sum of the negative natural-log-fidelities of the gates.
   For example, an errorless gate will have `0` for this, since $\ln(1) = 0$, while a gate with 0 fidelity would have a score of $\infty$.

`LogFidelity` has further subtleties related to the homogenisation of the `Target` (see more in the "Efficiency control" section below).
The `bin` argument is used to control when two of the same instructions in the `Target` on different `qargs` (say a `cx` on `(0, 1)` and `(2, 3)`) can be considered equivalent.
For each instruction with a fidelity $f = 1 - \epsilon$, we calculate its "binned" infidelity by rounding the _infidelity_ ($\epsilon$) in logarithmic space to the nearest multiple of `bin_width`, then calculating the fidelity from this.
Concretely, the binned infidelity is:

$$
\epsilon_{\text{binned}} = \begin{cases}
    1 & \text{if $f = 0$}\\
    0 & \text{if $f = 1$ or $\text{bin width} = \infty$}\\
    1 - f & \text{if $\text{bin width} = 0$}\\
    \exp\Bigl[\mathrm{round}\bigl(\frac{\ln\epsilon}{\text{bin width}}\bigr) \times \text{bin width}\Bigr]
        & \text{otherwise}.
\end{cases}
$$

After `Target` homogenisation is done using the binned infidelities, the "error" of the homogenized instruction is the geometric mean of the _actual_ infidelities of the instructions considered equal ($\bar\epsilon$).
The actual value used by the score is $-\ln(1 - \bar\epsilon)$, since this has the additive and minimising properties required of a score component.

A few things to note:

- Sensible human values of `bin_width` are things like $k\ln(10)$ for integer $k$; this rounds fidelities to "same number of nines" geometrically.
  For example, fidelities of $0.993$ and $0.94$ will both be binned to $0.99$, since geometrically
  their infidelities are closest to $1\%$.
- Setting `bin_width` to 0 effectively disables binning, while setting it to $\infty$ causes the
  infidelity to be ignored during binning, but the geometrically averaged values over all of the
  same instruction are still used in the rule scoring.
- The geometric average is used rather than the binned infidelity directly to avoid situations where
  all homogenised instructions are right at one edge of the bin but are scored as if they were
  centred on the bin.

### Efficiency control

Running the above algorithm in full generality can mean re-evaluating the full basis construction once for every physical qubit in the device, then again for every permutation of two-qubit `qargs` where at least one of the permutations is used, and so on.
This is not often entirely necessary; most `Target`s are _mostly_ homogeneous in terms of the gates available on each set of `qargs`, and most heterogeneity (in terms of the heuristic) comes only in the form of differing error rates between qubits.

We can trade off full accuracy in the heuristic to reduce the number of different parallel constructions of the algorithm we run.
This is done by running a "homogenisation" pass over the true `Target`, grouping `qargs` that will have the exact same search path.
The construction algorithm is a pure function of its inputs, and the only ways that the `Target` introduces heterogeneity are:

- different `qargs` having different operations defined
- the operations on different `qargs` having different fidelities to others

The homogenisation initially loosely groups based on the first bullet (e.g. if `(0,)` and `(1,)` both have exactly `sx` and `rz` defined, they are grouped, but if `(2,)` additionally has `x`, it gets its own group).
Depending on the optimisation level and the heuristic being used, a further grouping can take place based on the second bullet.
If error rates are not given or not taken into account, then the loose groups simply become the final groups.
Otherwise, we can take each operation's error rates, bin them (geometrically) according to some homogeneity target, and partition the loose groups into final groups that all have the same binned error rates on each operation.
In the limit of rounding to infinity this reduces the problem to the error-rate unaware case, and in the limit of no rounding the translation is handled considering the full heterogeneity of the system.

The optimal choice for the error-rate homogenisation likely depends on how aggressive and effective subsequent optimisation passes will be.
If the optimisation will be weak, cannot be relied on (such as if there are runtime-parametric gates in the source that will fail to resynthesise), or this translation pass is tidying up a heterogeneity-unaware optimisation routine, a higher error-rate heterogeneity consideration in the constructor may be appropriate.
If the translation will be followed by near-complete resynthesis of all gate chains, there is likely little need to spend time searching.

Unless `definition` fields of custom gates are used, the `BasisConstructor` algorithm proceeds largely independently of the source basis.
The source basis simply provides an early termination condition.
Because of this, a `BasisConstructor` can reuse its calculation state for many different circuits.


### Comparison to `BasisTranslator`

Compared to `BasisTranslator`, this pass is far more effective at handling heterogeneous target backends, and its search heuristic is far more controllable for determining the "best" translation.
This lets it use the same equivalence library to target various categories of backends efficiently: Clifford, Clifford+T, CX-like, controlled-Pauli-rotation-like, and iSwap, for example.
`BasisTranslator` requires additional passes to sort out directional gates, which gets ever more complex when dealing with multi-qubit gates.
Its search heuristic cannot be made to reliably find the minimal 2q-count or the minimal error translation of any given basis gate.

The Qiskit translator algorithm was designed at a time when backends were described in terms of loose, homogeneous constraints, and its search algorithm does not handle directional links or heterogeneity well.
`BasisTranslator`'s search is also formulated in a way that makes it hard or impossible to choose the _best_ translation (for any definition of "best") when many are available, though it is completely reliable at determining whether there is a valid translation, and finding one if so.



## Future work

Assuming everything goes as planned, this may well replace the default `BasisTranslator` as the translation method for all `Target`s at some point in the future (likely Qiskit 3), and potentially for limited subsets of `Target`s sooner (e.g. for fractional-gate backends).

When this does enter Qiskit proper, it will almost certainly be rewritten in Rust.
The current implementation is not optimised, and does not enjoy the low-level fast-path access to
the DAG data structures that other core Qiskit transpiler passes have.

## Copyright and licensing

This is a Qiskit project.
This is licensed under the terms of the Apache 2.0 licence.
