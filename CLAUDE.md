# CLAUDE.md

Guidance for Claude Code working in this repository.

## What this is

**CelineEngine** is the circuit simulator behind Céline and Céline Player: modified
nodal analysis, SPICE-grade device models, and a DK-method solver fast enough to run a
drawn circuit as an audio plugin. It knows nothing about drawing, plugins or the UI —
it is a netlist in, samples out.

It is not built on its own. Like CelineUI, it is pulled into each plugin as a submodule
at `modules/CelineEngine` and compiled inside that plugin's `SharedCode`
(`CelineEngine.cmake`). Its includes are by repository: `#include <CelineEngine/Engine.h>`.

- **Edit a checkout beside the plugins** (`../CelineEngine`) and build a plugin against
  it with `-DCELINE_ENGINE_DIR=../CelineEngine`, never inside a plugin's
  `modules/CelineEngine`: a submodule checks out a commit, not a branch.
- **Commit and push the engine first, then move each plugin's submodule** to that
  commit. A pin to an unpushed commit breaks CI and every fresh clone.
- **The tests live in Céline**, which builds against this checkout: `tests/CelineEngine.cpp`,
  `tests/ComponentModels.cpp` and the benchmarks there. Some of them build drawn sheets,
  which only Céline can do. The datasheets every model is fitted to are in Céline's
  `docs/specs/<kind>/`.
- A change here changes the sound of every sheet in both plugins. Céline's hidden
  `[.sheets]` bench renders real sheets against a reference and a tight-tolerance
  render; run it before and after anything in the solver.

## Layout

The sources are all in `CelineEngine/`; `LIMITATIONS.md` is at the root.

- `Components/` - one header per component type, each holding its parameters and its device models. Nothing here knows how to stamp or solve.
  - `Types.h` - `NodeIndex` / `ComponentId` aliases
  - `Ports.h` - the port abstraction every nonlinear device presents to the solver
  - `Junction.h` - p-n junction maths shared by every semiconductor (gmin, exponential I-V, SPICE voltage limiting)
  - `SpaceCharge.h` - the valve equivalent: 3/2-power Child-Langmuir flow, stable softplus, step limiting
  - `Resistor.h`, `Capacitor.h`, `Inductor.h`, `Potentiometer.h/.cpp` - linear parts
  - `VoltageSource.h` - DC supplies and AC windings (needs an MNA branch row, not a conductance)
  - `Diode.h`, `Transistor.h`, `Jfet.h`, `Mosfet.h` - semiconductors and their models. A `Diode` with `breakdownVoltage` set is a Zener; one with `transitTime` stores charge (integrated backward-Euler inside the device, zeroed for the bias solve); `seriesResistance` and a guard ring are wired by `addDiode()` on an internal node, as a BJT model flagged `darlington` is wired by `addTransistor()` as two devices. Every model added from the MOSFET on is fitted to a datasheet or taken from the vendor's card, and its header comment and `tests/ComponentModels.cpp` say which, and how close it lands. **The device equations are SPICE's**, checked in ngspice with identical parameters to 0.05 mV: the BJT's Early effect is Gummel-Poon's q1 = 1/(1 − Vbc/VAF), not the level-1 (1 + Vce/VAF) it once was (that read every vendor's VAF ~1% wrong), and the JFET has both gate junctions plus CGS/CGD. A vendor card is not automatically right — Fairchild's 2N3904 card gives 0.6× its own sheet's gain, onsemi's J201 card has a −3.95 V pinch-off — so models are fitted to the maker's *measured* curves (digitised from the PDF's vector paths) and the cards are what they are checked against. Every source used lives under Céline's `docs/specs/<kind>/` with a README saying what each file is for.
  - `Vccs.h` - voltage-controlled current source (linear); the primitive that makes amplifiers expressible
  - `IdealOpAmp.h` - a nullor: one constraint row, no ports, no saturation. Cheaper and exact where the op-amp's limits don't matter
  - `Transformer.h` - ideal N-winding transformer (one row per winding). Compose with an inductor across the primary for magnetising inductance, series L/R for leakage and losses — which is what the drawn transformer's **Real** model does for you, in `SchematicBuilder`, rather than here
  - `OpAmp.h` - op-amp parameters. **Not a device**: `addOpAmp()` assembles one from a resistor, two `Vccs`, a capacitor and two diodes (plus two `Vccs` for the input bias current), so the gain is linear and only the clamp diodes reach the Newton loop. See the header for why a single nonlinear device broke DK. A model with a `slewRate` swaps the input `Vccs` for a `DifferentialPair` (`DifferentialPair.h`): the Boyle macromodel's long-tailed pair, whose tail current over the compensation capacitor *is* the slew rate. It is a nonlinear element on the differential input — the shape that once broke DK — so its step limiter matters and a test holds DK and full Newton to one answer on a hard-driven LM308. It costs two Newton ports, so `SchematicBuilder` fits it only where `OpAmpModel::slewsAt()` says the part can slew at the circuit's running rate (a sample can ask at most swing × rate of it): on 9 V at 48 kHz only the LM308 qualifies. Checked against TI's own µA741 Boyle card in ngspice; the input common-mode range and output current limit are deliberately not modelled, for reasons in `LIMITATIONS.md`. Each model also carries a `phaseMargin`, which adds a second pole on a buffered RC node (linear: a node, no Newton port) and raises the dominant pole so the gain crosses unity at `gainBandwidth`. **`gainBandwidth` is the datasheet's unity-gain *crossing*, not a follower's closed-loop bandwidth** — reading the 741's rise time as the latter once left it 8–11% too loud at 10–20 kHz. Céline's `tests/CelineEngine.cpp` holds the TL072 and LM741 to TI's cards on a gain-101 stage, amplitude and phase. **The rail clamp's saturation current is scaled by the input transconductance** and its knee is a quarter of a junction's, so it conducts by the overdrive in volts, not by the gain node's arbitrary current scale. A plain 1N4148 there let the clip level drift 0.4 V with drive. Headroom is quoted at ±15 V into `headroomLoad` and moved per volt of supply by `headroomPerSupplyVolt` only where a sheet plots it (the TL072, the 741). The 741's two sides differ, so `negativeRailHeadroom` gives the negative side its own figure, and its open-loop gain follows the supply through `openLoopGainSupplyExponent` (use `openLoopGainOnRails()`, never the raw ±15 V field). Every new figure has to hold at both 9 V and ±15 V, which the clip-level test runs. **The drawn op-amp takes both rails**: *Supply +* and *Supply −*. The second is a rail (`ElementInfo::valueBIsRail`), so zero is an answer — ground, a single supply — rather than an unfilled value: it is neither drawn red nor refused, a sheet saved before the field reads it as zero, and whatever is typed is taken as that many volts below ground. A ±15 V preamp is drawn as its schematic draws it, references at ground; the old workaround, 30 V around a 15 V virtual ground, is the same circuit to within the 1.5 mV the bias solve's gmin puts between them, which "A drawn op-amp runs on a split supply" in `tests/ComponentModels.cpp` holds. `inputBiasCurrent` is signed by the input pair (into the pins for NPN, out for the 4558's PNP) and drawn by two `Vccs` reading the clamp rails, because a rail of its own cost 7% on the two-op-amp pedal. The 9 V clipping shape is held to TI's cards in Céline's `tests/CelineEngine.cpp` too.
  - `PowerAmp.h` - the LM386, assembled the same way by `addPowerAmp()`: the datasheet's internal resistor network literally, ideal followers holding the pair emitters, and a gain node rolled off by one capacitor. Its gain node sits at a high impedance on purpose — see `gainNodeResistance` — so an overdriven loop hands the clamps milliamps rather than amps.
  - `VacuumDiode.h`, `Triode.h`, `Pentode.h` - valves, Koren models. **Every model is fitted here to a JJ datasheet and reproduces Ia, gm and ra at the quoted operating point** — all four triodes, all seven pentodes, all three rectifiers. The pentodes add one term to Koren's: a plate slope, `Ip *= 1 + Vpk/VA` (`plateSlopeVoltage`), because Koren's single arctangent cannot set both the knee and the plate resistance. Six of the seven are refitted with it — operating point exact, knee and slope fitted to their sheets' plate families, digitised from the vector paths — and `tests/ComponentModels.cpp` embeds those families. The 6V6S keeps Koren's form: its sheet has no pentode family to fit to. `plateSlopeVoltage` sits *after* `gridPerveance` in the struct on purpose — the models use aggregate initialisation in field order. Koren's own published sets were 20–30% out at that point on the triodes and had plate resistance 2–6× too high on the pentodes, so "it's Koren's set" is not evidence of anything; check before trusting. Each model's header comment records what it was fitted to and how well.

  Triodes carry two formulations, selected by `TriodeModel::formulation`. **Koren** fits published curves and is what the datasheet-derived models use. **Dempwolf-Zolzer** (DAFx-11) is fitted to *measurements* of three individual 12AX7s, and is the only path with a real grid current — Koren has none, so this project bolts a guessed perveance onto it, and that guess draws 60% too much at Vg = +1 V. D-Z costs about 2% more per sample. Both share the same ports, so neither needs solver support; `linearise()` branches once. What fitting a point still doesn't buy you is in **`LIMITATIONS.md`** under "Valves" — read it before trusting a valve result.
  - `FastMath.h` - approximations of `exp`, `log` and `pow`, and the process-wide switch that turns them on. Accurate to ~3e-10 relative, four orders tighter than the Newton tolerance they feed. Whether they are *faster* is a platform question: measured a wash on Apple Silicon, because Apple's libm is already excellent. `sqrt` is deliberately absent — it is one instruction.
  - `Components.h` - umbrella include
- `Solver.h/.cpp` - dense LU with partial pivoting, allocation-free once sized. Stores the reciprocal of each pivot so back substitution multiplies rather than divides. **Both substitutions run a column at a time**, from a transposed copy of the factors made at the end of `factorise()`: written a row at a time, as dot products, every multiply-add waited on the one before it, and on a 19-port sheet that made the solve the costliest function in the engine. `factoriseAndSolve()` is for a matrix used once — the Newton Jacobian — and carries the right-hand side through the elimination instead. The DK Newton loop builds its residual and Jacobian a column at a time for the same reason. Measured with the hidden `[.sheets]` bench (`tests/SheetBench.cpp`, git-ignored, reads `../_test schematics`): 11–22% off three real sheets, output within −150 dB of before.
- `Engine.h` - `Circuit`, and the class comment explaining the DK method. Read this first.
- `Engine.cpp` - netlist construction, component updates, node naming, prepare()/reset()
- `EngineSolve.cpp` - stamping, DK precomputation, the Newton loops, per-sample process()
- `Circuits.h/.cpp` - **the engine's test bench, not how the plugin builds circuits.** The plugin builds from the drawing, through `SchematicBuilder`; nothing in Céline's `source/` calls anything here. What keeps it is Céline's `tests/CelineEngine.cpp` and `benchmarks/`, which need circuits whose behaviour is known independently of this codebase — a one-pole whose corner is arithmetic, a tone stack with a published curve, an op-amp overdrive transcribed from a SPICE netlist. A drawn schematic can't play that role: it would only confirm that the builder and the engine agree with each other. Anything here with no caller is dead and should go.
- **`LIMITATIONS.md` - what the simulation does not model, ordered by audible impact. Read before trusting a result or blaming a circuit.**

### Adding a nonlinear device

Nonlinear devices are described entirely by their **ports** (see `Components/Ports.h`) — node pairs whose voltage drives the device and through which its current flows. A device provides four free functions in `CircuitComponents`:

```cpp
constexpr int portCount(const MyDevice&);
void fillPorts(const MyDevice&, Port* out);
bool limitPortVoltages(MyDevice&, double* v);   // damp the Newton step
void linearise(const MyDevice&, const double* v, DeviceLinearisation& out);
```

Then add a `std::vector<MyDevice>` to `Circuit`, an `add*` function, and one line to `forEachNonlinearDevice()` in `EngineSolve.cpp`. That single function defines the port ordering for every pass that walks devices in step (`buildPortList`, `limitPortVoltages`, `linearisePorts`), so there is one place to edit rather than three that have to agree. `clearState()` in `Engine.cpp` still needs its own entry, since each device type resets different members. No stamping or solver code is needed — the solver only ever sees ports.

### Stamping: one description, two systems

The cached per-sample system and the DC operating-point system are the same circuit seen two ways. `stampTopology()` in `EngineSolve.cpp` is written once for everything that stamps identically into both (resistors, `Vccs`, ideal op-amps, transformers, and a voltage source's ±1 pattern), and is templated on what to do with a terminal whose voltage is already known — the input node's coupling in one case, a fixed RHS contribution in the other. Only the genuinely physical differences stay in the callers: capacitors are open at DC, inductors are shorts, and a source's voltage lives in the per-sample RHS so an AC winding can vary without touching the cached matrix. Keep it that way; two copies that must agree exactly is how a bias point quietly stops describing the circuit that actually runs.

### Solver strategy

Two paths, chosen per circuit at `prepare()` time by `SolverStrategy::Auto`:

- **DK** eliminates the linear part once and runs Newton on the `m` ports. Wins when `m < n`, which is the case for anything valve-shaped (lots of linear plumbing, few nonlinear devices).
- **Full Newton** stamps devices into the whole `n × n` and re-factorises each iteration. Wins for small circuits where `m >= n`, e.g. a shunt diode clipper (1 unknown node, 2 ports).

Because DK lifts the nonlinear devices out of the linear matrix, **every port gets `gmin` stamped across it there as well**. Without that, a node whose only company is semiconductor terminals — a transistor's collector wired straight to the next one's base, an ordinary direct-coupled pair — has no conductance in the linear system at all and the factorisation is singular. That failure was silent: the DC system stamps the devices properly, so the bias point solved and looked right, while `process()` saw an invalid factorisation and returned 0 for every sample. This is what SPICE's gmin is for, and the device models already add the same 1e-12 inside their own port currents.

Components that state a *constraint* rather than a current — `VoltageSource`, `IdealOpAmp`, `Transformer` — take extra MNA rows after the node rows, allocated in `prepare()` via `idealOpAmpRowOffset` / `transformerRowOffset`. They are linear, so they cost nothing per sample beyond matrix size.

Force one with `setSolverStrategy()` to benchmark or to check the two agree — they match to ~1e-8 V.

## Code quality

Warnings are errors, in both plugins that build this. Uses `.clang-format`-style
Allman braces, 4-space indentation, no column limit — match the surrounding code.

Everything on the per-sample path is realtime code: nothing allocates after
`prepare()`, nothing locks, nothing grows a container.
