//==============================================================================
/*
    The solving half of Circuit. Netlist construction, component updates and
    lifetime live in Engine.cpp; everything here is about turning that netlist
    into numbers, once per sample.

    Read Engine.h's class comment first -- it explains the DK method these
    functions implement. In short: the linear part of the circuit is eliminated
    ahead of time so that Newton iterates on the handful of nonlinear ports
    rather than on every node.
*/
//==============================================================================

#include "Engine.h"

#include <algorithm>
#include <cmath>

namespace
{
    /** Conductance an inductor is given when solving the DC operating point,
        where it should ideally be a dead short. Small enough to be a short next
        to anything else in an audio circuit, large enough not to wreck the
        matrix conditioning. */
    constexpr double dcInductorConductance = 1.0e3; // 1 milliohm
} // namespace

//==============================================================================
// Stamping
//
// Two systems -- the cached per-sample one and the DC one the operating point
// uses -- are the same circuit seen two ways. Only three things differ, each
// for a physical reason:
//
//     capacitors    a conductance per timestep here, an open circuit at DC
//     inductors     a conductance per timestep here, a short at DC
//     sources       the voltage lands in the per-sample right-hand side, so an
//                   AC winding can vary without touching the cached matrix
//
// Everything else is written once below and used by both, parameterised on the
// one remaining difference: what becomes of a terminal whose voltage is already
// known -- the input node's coupling in one case, a fixed right-hand side
// contribution in the other.
//
// Sharing it is not tidiness. Two copies that must agree exactly is how a bias
// point quietly stops describing the circuit that actually runs, which produces
// a plausible wrong answer rather than an obvious one.
//==============================================================================

template <typename KnownNodeTerm>
void Circuit::stampConductance(double* matrix, NodeIndex a, NodeIndex b, double g, KnownNodeTerm knownNodeTerm) noexcept
{
    const int ra = rowOf(a);
    const int rb = rowOf(b);

    if (ra >= 0)
        matrix[ra * systemSize + ra] += g;
    if (rb >= 0)
        matrix[rb * systemSize + rb] += g;

    if (ra >= 0 && rb >= 0)
    {
        matrix[ra * systemSize + rb] -= g;
        matrix[rb * systemSize + ra] -= g;
    }
    else if (ra >= 0)
    {
        // b is a known voltage, so its contribution is a constant current source.
        knownNodeTerm(ra, b, g);
    }
    else if (rb >= 0)
    {
        knownNodeTerm(rb, a, g);
    }
}

template <typename KnownNodeTerm>
void Circuit::stampTopology(double* matrix, KnownNodeTerm knownNodeTerm) noexcept
{
    for (const auto& r : resistors)
        stampConductance(matrix, r.a, r.b, 1.0 / r.ohms, knownNodeTerm);

    // Controlled sources: the current appears in one pair of rows and the
    // voltage driving it in a different pair of columns. That off-diagonal
    // placement is what makes an amplifier an amplifier, and it is linear.
    for (const auto& g : transconductances)
    {
        const int rFrom = rowOf(g.from);
        const int rTo = rowOf(g.to);

        for (int sign = 0; sign < 2; ++sign)
        {
            const NodeIndex control = sign == 0 ? g.controlPositive : g.controlNegative;
            const double gm = sign == 0 ? g.transconductance : -g.transconductance;
            const int rControl = rowOf(control);

            if (rControl >= 0)
            {
                if (rFrom >= 0)
                    matrix[rFrom * systemSize + rControl] += gm;
                if (rTo >= 0)
                    matrix[rTo * systemSize + rControl] -= gm;
            }
            else
            {
                // The controlling node's voltage is known, so its contribution
                // is a fixed current rather than a matrix entry.
                if (rFrom >= 0)
                    knownNodeTerm(rFrom, control, -gm);
                if (rTo >= 0)
                    knownNodeTerm(rTo, control, gm);
            }
        }
    }

    // Ideal op-amps: the output current is the unknown and the row holds the
    // inputs equal. Asymmetric against a voltage source -- the constraint is on
    // one pair of nodes, the current at a different one.
    for (size_t a = 0; a < idealOpAmps.size(); ++a)
    {
        const auto& amp = idealOpAmps[a];
        const int ampRow = idealOpAmpRowOffset + static_cast<int>(a);
        const int rOut = rowOf(amp.output);

        // The op-amp sources this current into the output node.
        if (rOut >= 0)
            matrix[rOut * systemSize + ampRow] -= 1.0;

        const int rPlus = rowOf(amp.inPlus);
        const int rMinus = rowOf(amp.inMinus);

        if (rPlus >= 0)
            matrix[ampRow * systemSize + rPlus] += 1.0;
        else
            knownNodeTerm(ampRow, amp.inPlus, -1.0);

        if (rMinus >= 0)
            matrix[ampRow * systemSize + rMinus] -= 1.0;
        else
            knownNodeTerm(ampRow, amp.inMinus, 1.0);
    }

    // Transformers. One unknown current per winding, and the same number of
    // equations: every winding shares the volts per turn, and the magnetomotive
    // forces sum to zero.
    int windingRow = transformerRowOffset;

    for (const auto& t : transformers)
    {
        const int firstRow = windingRow;

        for (int w = 0; w < t.windingCount; ++w)
        {
            const int ra = rowOf(t.windings[w].a);
            const int rb = rowOf(t.windings[w].b);

            if (ra >= 0)
                matrix[ra * systemSize + firstRow + w] += 1.0;
            if (rb >= 0)
                matrix[rb * systemSize + firstRow + w] -= 1.0;
        }

        // V(0)/N(0) = V(k)/N(k), written as N(k)*V(0) - N(0)*V(k) = 0 so there's
        // nothing to divide by.
        const auto& first = t.windings[0];

        for (int w = 1; w < t.windingCount; ++w)
        {
            const auto& other = t.windings[w];
            const int equation = firstRow + w - 1;

            auto term = [&](NodeIndex node, double coefficient)
            {
                const int r = rowOf(node);
                if (r >= 0)
                    matrix[equation * systemSize + r] += coefficient;
                else
                    knownNodeTerm(equation, node, -coefficient);
            };

            term(first.a, other.turns);
            term(first.b, -other.turns);
            term(other.a, -first.turns);
            term(other.b, first.turns);
        }

        // Sum of N(k)*I(k) = 0 -- an ideal core takes no magnetising current.
        const int ampereRow = firstRow + t.windingCount - 1;
        for (int w = 0; w < t.windingCount; ++w)
            matrix[ampereRow * systemSize + firstRow + w] += t.windings[w].turns;

        windingRow += t.windingCount;
    }

    // Voltage sources: the +/-1 pattern only. The voltage itself is the caller's
    // business -- see the note at the top of this section.
    for (size_t s = 0; s < voltageSources.size(); ++s)
    {
        const auto& source = voltageSources[s];
        const int sourceRow = numNodeUnknowns + static_cast<int>(s);
        const int rp = rowOf(source.positive);
        const int rn = rowOf(source.negative);

        // The branch current appears in the two node equations...
        if (rp >= 0)
            matrix[rp * systemSize + sourceRow] += 1.0;
        if (rn >= 0)
            matrix[rn * systemSize + sourceRow] -= 1.0;

        // ...and the source's own row states V(positive) - V(negative) = volts.
        if (rp >= 0)
            matrix[sourceRow * systemSize + rp] += 1.0;
        else
            knownNodeTerm(sourceRow, source.positive, -1.0);

        if (rn >= 0)
            matrix[sourceRow * systemSize + rn] -= 1.0;
        else
            knownNodeTerm(sourceRow, source.negative, 1.0);
    }
}

void Circuit::stampCurrentSource(double* rhsVector, NodeIndex a, NodeIndex b, double ieq) noexcept
{
    const int ra = rowOf(a);
    const int rb = rowOf(b);

    if (ra >= 0)
        rhsVector[ra] += ieq;
    if (rb >= 0)
        rhsVector[rb] -= ieq;
}

void Circuit::rebuildLinearSystem()
{
    std::fill(linearMatrix.begin(), linearMatrix.end(), 0.0);
    std::fill(inputCoupling.begin(), inputCoupling.end(), 0.0);

    double* matrix = linearMatrix.data();

    // Ground is 0 V so it contributes nothing at all; the input node changes
    // every sample, so its coefficient is kept and multiplied by vIn in
    // buildRightHandSide().
    auto toInputCoupling = [this](int row, NodeIndex node, double coefficient)
    {
        if (node == inputIndex)
            inputCoupling[static_cast<size_t>(row)] += coefficient;
    };

    for (auto& c : capacitors)
    {
        // A capacitor in series with its ESR: over one trapezoidal step the
        // capacitance looks like dt/2C ohms, so the branch is just those two
        // resistances in series. With no ESR this is the usual 2C/dt.
        c.halfStepResistance = dt / (2.0 * c.farads);
        c.conductance = 1.0 / (c.esrOhms + c.halfStepResistance);
        stampConductance(matrix, c.a, c.b, c.conductance, toInputCoupling);
    }

    for (auto& l : inductors)
    {
        l.conductance = dt / (2.0 * l.henries);
        stampConductance(matrix, l.a, l.b, l.conductance, toInputCoupling);
    }

    stampTopology(matrix, toInputCoupling);

    // gmin across every nonlinear port.
    //
    // DK lifts the nonlinear devices out of this matrix, so a node whose only
    // company is semiconductor terminals has no conductance in the linear system
    // at all and the factorisation is singular. Not an exotic case: a
    // transistor's collector wired to the next one's base is an ordinary
    // direct-coupled pair.
    //
    // The failure was invisible, which is the worse part -- the DC system stamps
    // the devices properly, so the operating point solved and looked right while
    // process() returned silence. This is what SPICE's gmin is for.
    for (const auto& port : ports)
        stampConductance(matrix, port.a, port.b, CircuitComponents::gmin, toInputCoupling);

    // Which rows vIn actually reaches -- typically one or two out of however
    // many the circuit has. Walking those rather than the whole vector is what
    // lets buildRightHandSide() start from a straight fill.
    inputCouplingRows.clear();

    for (int row = 0; row < systemSize; ++row)
        if (inputCoupling[static_cast<size_t>(row)] != 0.0)
            inputCouplingRows.push_back({row, inputCoupling[static_cast<size_t>(row)]});

    // The matrix is now final until the next value change. Factorise it once
    // here and every sample until then is just substitutions -- and, if there
    // are nonlinear devices, so is the DK precomputation that follows.
    factorisationValid = false;

    if (systemSize > 0)
    {
        std::copy(linearMatrix.begin(), linearMatrix.end(), solver.data());
        // Not asserted on. A singular netlist is a drawing somebody made, not a
        // programming error: SchematicBuilder checks hasUsableFactorisation()
        // and names the problem in the console. The assertion that used to be
        // here fired during the test suite's own failure-path cases, which is
        // an assertion that has stopped meaning anything.
        factorisationValid = solver.factorise();

        if (factorisationValid && portCount > 0)
            precomputeDkSystem();
    }

    linearDirty = false;
}

double Circuit::sourceVoltage(const VoltageSource& source) const noexcept
{
    return source.volts + source.acAmplitude * std::sin(source.phase);
}

void Circuit::buildDcSystem(double sourceScale) noexcept
{
    std::fill(dcMatrix.begin(), dcMatrix.end(), 0.0);
    std::fill(dcRhs.begin(), dcRhs.end(), 0.0);

    double* matrix = dcMatrix.data();
    double* rhsVector = dcRhs.data();

    // A known terminal is a fixed voltage, so it becomes a fixed current on the
    // right-hand side. Both known nodes sit at 0 V in practice, so these terms
    // vanish -- written out anyway, since the DC solve has no business assuming
    // it.
    auto toDcRhs = [this, rhsVector](int row, NodeIndex node, double coefficient)
    {
        rhsVector[row] += coefficient * nodeVoltage[static_cast<size_t>(node)];
    };

    // Capacitors are open circuits at DC, so they aren't stamped at all.
    // Inductors are shorts.
    for (const auto& l : inductors)
        stampConductance(matrix, l.a, l.b, dcInductorConductance, toDcRhs);

    stampTopology(matrix, toDcRhs);

    // Each source's row states V(positive) - V(negative) = volts, ramped by
    // sourceScale so the operating point can bring the supplies up gradually.
    //
    // A sinusoidal winding is taken at its peak, not its average: the average is
    // zero, which would leave the reservoir capacitor flat and open playback
    // with a tenth of a second of the rail charging. The peak is where it
    // settles anyway.
    for (size_t s = 0; s < voltageSources.size(); ++s)
    {
        const auto& source = voltageSources[s];
        rhsVector[static_cast<size_t>(numNodeUnknowns) + s] +=
            (source.volts + source.acAmplitude) * sourceScale;
    }

    // Without this, any node reachable only through capacitors would have an
    // all-zero row.
    for (int row = 0; row < numNodeUnknowns; ++row)
        matrix[row * systemSize + row] += CircuitComponents::gmin;
}

//==============================================================================
// Nonlinear ports
//==============================================================================

template <typename Fn>
void Circuit::forEachNonlinearDevice(Fn fn)
{
    // The one place the order of the nonlinear devices is decided, and with it
    // the layout of the port list. Three passes below -- build, limit,
    // linearise -- walk devices in step and line up only because they all come
    // through here. Hand-written, that agreement is an invariant maintained in
    // three places.
    //
    // The running port offset is handed to the callback rather than tracked by
    // it, so the whole thing inlines to what the hand-written loops did.
    int firstPort = 0;

    auto visit = [&fn, &firstPort](auto& device)
    {
        fn(device, firstPort);
        firstPort += CircuitComponents::portCount(device);
    };

    for (auto& d : diodes)
        visit(d);

    for (auto& t : transistors)
        visit(t);

    for (auto& d : vacuumDiodes)
        visit(d);

    for (auto& t : triodes)
        visit(t);

    for (auto& p : pentodes)
        visit(p);

    for (auto& j : jfets)
        visit(j);

    for (auto& m : mosfets)
        visit(m);

    for (auto& d : differentialPairs)
        visit(d);
}

void Circuit::buildPortList()
{
    ports.clear();
    deviceBlocks.clear();

    CircuitComponents::Port scratch[CircuitComponents::maxPortsPerDevice];

    forEachNonlinearDevice(
        [&](const auto& device, int firstPort)
        {
            juce::ignoreUnused(firstPort); // ports.size() is the same thing here
            const int n = CircuitComponents::portCount(device);
            CircuitComponents::fillPorts(device, scratch);
            deviceBlocks.emplace_back(static_cast<int>(ports.size()), n);
            ports.insert(ports.end(), scratch, scratch + n);
        });

    portCount = static_cast<int>(ports.size());

    const auto m = static_cast<size_t>(portCount);
    portRowA.assign(m, -1);
    portRowB.assign(m, -1);
    portInputCoupling.assign(m, 0.0);
    portVoltage.assign(m, 0.0);
    portVoltageStart.assign(m, 0.0);
    portOpenTarget.assign(m, 0.0);
    portOpenStart.assign(m, 0.0);
    portVoltageReached.assign(m, 0.0);
    portCurrent.assign(m, 0.0);
    portOpenVoltage.assign(m, 0.0);
    portResidual.assign(m, 0.0);
    portJacobian.assign(m * m, 0.0);
    portImpedance.assign(m * m, 0.0);
    portImpedanceByColumn.assign(m * m, 0.0);
    portJacobianColumn.assign(m, 0.0);
    portTransfer.assign(static_cast<size_t>(systemSize) * m, 0.0);
    portSolver.resize(portCount);

    portBlockStart.assign(m, 0);
    portBlockSize.assign(m, 0);

    for (const auto& [firstPort, blockSize] : deviceBlocks)
        for (int p = firstPort; p < firstPort + blockSize; ++p)
        {
            portBlockStart[static_cast<size_t>(p)] = firstPort;
            portBlockSize[static_cast<size_t>(p)] = blockSize;
        }

    activeRows.assign(m, 0);
    activeSystem.assign(m * m, 0.0);
    activeRhs.assign(m, 0.0);
    portImpedanceColumnMax.assign(m, 0.0);

    for (int p = 0; p < portCount; ++p)
    {
        const auto& port = ports[static_cast<size_t>(p)];
        portRowA[static_cast<size_t>(p)] = rowOf(port.a);
        portRowB[static_cast<size_t>(p)] = rowOf(port.b);

        // A port terminal sitting on the input node picks up vIn directly rather
        // than through the matrix, since that node isn't an unknown.
        double coupling = 0.0;
        if (port.a == inputIndex)
            coupling += 1.0;
        if (port.b == inputIndex)
            coupling -= 1.0;
        portInputCoupling[static_cast<size_t>(p)] = coupling;
    }
}

bool Circuit::limitPortVoltages() noexcept
{
    bool limited = false;

    forEachNonlinearDevice(
        [&limited, this](auto& device, int firstPort)
        {
            limited |= CircuitComponents::limitPortVoltages(device, portVoltage.data() + firstPort);
        });

    return limited;
}

void Circuit::settleLimiters(const std::vector<double>& target) noexcept
{
    // Each device's limiter damps a step against the voltage it last
    // linearised around, which it remembers. Restarting Newton from a known
    // point needs that memory moved there first: otherwise the first iteration
    // clamps the start towards wherever the abandoned attempt was -- up to ten
    // volts on a grid, a hundred on a plate -- and the restart begins from a
    // point that solves nothing. Measured: continuation failing its very first
    // sub-step, a thousandth of a sample from an exact solution.
    //
    // Limiting towards the target until no device acts moves every memory
    // there by each device's own rules, polarity and reflections included,
    // with nothing device-specific here. A limiter that does not act has
    // recorded exactly the voltage it was given -- which every limiter must
    // reach for Newton to converge at all, so this ends.
    for (int pass = 0; pass < maxLimiterSettlingPasses; ++pass)
    {
        std::copy(target.begin(), target.end(), portVoltage.begin());

        if (! limitPortVoltages())
            break;
    }

    std::copy(target.begin(), target.end(), portVoltage.begin());
}

void Circuit::linearisePorts() noexcept
{
    // portJacobian is deliberately not zeroed: every device overwrites its own
    // diagonal block and nothing reads outside one, since stampPorts() and the
    // K*di/dv product both walk deviceBlocks.
    CircuitComponents::DeviceLinearisation linearised{};

    // Each device fills a square block on the diagonal: its ports interact with
    // each other and with nothing else. The block-diagonal structure is what
    // keeps this cheap however many devices there are.
    forEachNonlinearDevice(
        [&linearised, this](const auto& device, int firstPort)
        {
            const int n = CircuitComponents::portCount(device);
            CircuitComponents::linearise(device, portVoltage.data() + firstPort, linearised);

            for (int r = 0; r < n; ++r)
            {
                portCurrent[static_cast<size_t>(firstPort + r)] = linearised.current[r];

                for (int c = 0; c < n; ++c)
                    portJacobian[static_cast<size_t>((firstPort + r) * portCount + firstPort + c)] =
                        linearised.jacobian[r * n + c];
            }
        });
}

void Circuit::readPortVoltagesFromNodes() noexcept
{
    for (int p = 0; p < portCount; ++p)
    {
        const auto& port = ports[static_cast<size_t>(p)];
        portVoltage[static_cast<size_t>(p)] =
            nodeVoltage[static_cast<size_t>(port.a)] - nodeVoltage[static_cast<size_t>(port.b)];
    }
}

void Circuit::stampPorts(double* matrix, double* rhsVector) noexcept
{
    // Port p carries current i_p from its node a to its node b. Linearised,
    //     i_p = sum_q J[p][q] * v_q + ieq_p
    // and each v_q is itself a node-voltage difference, so one Jacobian entry
    // lands in up to four places in the matrix.
    for (const auto& [firstPort, blockSize] : deviceBlocks)
    for (int p = firstPort; p < firstPort + blockSize; ++p)
    {
        const int ra = portRowA[static_cast<size_t>(p)];
        const int rb = portRowB[static_cast<size_t>(p)];

        double equivalentCurrent = portCurrent[static_cast<size_t>(p)];

        // Only this device's own ports can appear in row p; the rest of the row
        // is structurally zero.
        for (int q = firstPort; q < firstPort + blockSize; ++q)
        {
            const double j = portJacobian[static_cast<size_t>(p * portCount + q)];

            const auto& portQ = ports[static_cast<size_t>(q)];
            const int ca = portRowA[static_cast<size_t>(q)];
            const int cb = portRowB[static_cast<size_t>(q)];

            // Move the tangent's constant term to the right-hand side. Written
            // against the limited port voltage, which is where the tangent was
            // actually taken, not against whatever the nodes currently say.
            equivalentCurrent -= j * portVoltage[static_cast<size_t>(q)];

            if (ra >= 0)
            {
                if (ca >= 0)
                    matrix[ra * systemSize + ca] += j;
                else
                    rhsVector[ra] -= j * nodeVoltage[static_cast<size_t>(portQ.a)];

                if (cb >= 0)
                    matrix[ra * systemSize + cb] -= j;
                else
                    rhsVector[ra] += j * nodeVoltage[static_cast<size_t>(portQ.b)];
            }

            if (rb >= 0)
            {
                if (ca >= 0)
                    matrix[rb * systemSize + ca] -= j;
                else
                    rhsVector[rb] += j * nodeVoltage[static_cast<size_t>(portQ.a)];

                if (cb >= 0)
                    matrix[rb * systemSize + cb] += j;
                else
                    rhsVector[rb] -= j * nodeVoltage[static_cast<size_t>(portQ.b)];
            }
        }

        if (ra >= 0)
            rhsVector[ra] -= equivalentCurrent;
        if (rb >= 0)
            rhsVector[rb] += equivalentCurrent;
    }
}

//==============================================================================
// DK precomputation
//==============================================================================

void Circuit::precomputeDkSystem()
{
    // W = A^-1 N'. One triangular solve per port against the already-factorised
    // linear matrix, so O(n^2) each -- the expensive part, done once per value
    // change rather than once per sample.
    for (int p = 0; p < portCount; ++p)
    {
        double* column = portTransfer.data() + static_cast<size_t>(p) * static_cast<size_t>(systemSize);
        std::fill(column, column + systemSize, 0.0);

        const int ra = portRowA[static_cast<size_t>(p)];
        const int rb = portRowB[static_cast<size_t>(p)];

        if (ra >= 0)
            column[ra] += 1.0;
        if (rb >= 0)
            column[rb] -= 1.0;

        solver.solveInPlace(column);
    }

    // K = N W: the impedance the circuit presents between each pair of ports.
    for (int p = 0; p < portCount; ++p)
    {
        const int ra = portRowA[static_cast<size_t>(p)];
        const int rb = portRowB[static_cast<size_t>(p)];

        for (int q = 0; q < portCount; ++q)
        {
            const double* column = portTransfer.data() + static_cast<size_t>(q) * static_cast<size_t>(systemSize);
            double k = 0.0;

            if (ra >= 0)
                k += column[ra];
            if (rb >= 0)
                k -= column[rb];

            portImpedance[static_cast<size_t>(p * portCount + q)] = k;
            portImpedanceByColumn[static_cast<size_t>(q * portCount + p)] = k;
        }
    }

    for (int q = 0; q < portCount; ++q)
    {
        const double* column = portImpedanceByColumn.data() + static_cast<size_t>(q * portCount);
        double largest = 0.0;

        for (int p = 0; p < portCount; ++p)
            largest = std::max(largest, std::abs(column[p]));

        portImpedanceColumnMax[static_cast<size_t>(q)] = largest;
    }
}

//==============================================================================
// Solving
//==============================================================================

bool Circuit::writeBackSolution(bool checkConvergence) noexcept
{
    bool converged = true;

    for (int row = 0; row < numNodeUnknowns; ++row)
    {
        const auto node = static_cast<size_t>(nodeOfRow[static_cast<size_t>(row)]);
        const double solved = work[static_cast<size_t>(row)];

        // DK decides convergence on the ports, not the nodes, so it asks for this
        // to be skipped -- it's a per-node abs, max and compare for an answer
        // nobody reads.
        if (checkConvergence)
        {
            const double previous = nodeVoltage[node];
            const double tolerance = convergenceAbsTolerance
                                   + convergenceRelTolerance * std::max(std::abs(previous), std::abs(solved));

            if (std::abs(solved - previous) > tolerance)
                converged = false;
        }

        nodeVoltage[node] = solved;
    }

    // The branch currents follow linearly from the node voltages, so they don't
    // get a say in convergence -- once the voltages settle, these have too.
    for (size_t s = 0; s < voltageSources.size(); ++s)
        voltageSources[s].current = work[static_cast<size_t>(numNodeUnknowns) + s];

    for (size_t a = 0; a < idealOpAmps.size(); ++a)
        idealOpAmps[a].current = work[static_cast<size_t>(idealOpAmpRowOffset) + a];

    {
        int windingRow = transformerRowOffset;

        for (auto& t : transformers)
            for (int w = 0; w < t.windingCount; ++w)
                t.windings[w].current = work[static_cast<size_t>(windingRow++)];
    }

    return converged;
}

int Circuit::solveReducedNewtonStep() noexcept
{
    // J = I + K D, with D = di/dv block diagonal. A current row of D that is
    // all but zero -- a grid that is not conducting, a junction that is off --
    // adds nothing to J, so D is really just the k rows that are not, and
    //
    //     J^-1 = I - K_a (I_k + D_a K_a)^-1 D_a        (Woodbury)
    //
    // turns the m x m solve into a k x k one. On a valve amp most grids
    // conduct only on peaks, so k is about two thirds of m and the
    // factorisation a third or less of what it was.
    //
    // A row is left out only when nothing it adds to any entry of J reaches
    // activeRowThreshold, against J's unit diagonal. The step is then a
    // Newton step on a Jacobian that much from exact: still quadratic, still
    // the same root. Measured on a full VH4 against the full step: within
    // -143 dB of peak, the same iteration count to the sample.
    const auto m = static_cast<size_t>(portCount);
    int k = 0;

    for (int r = 0; r < portCount; ++r)
    {
        const auto row = static_cast<size_t>(r);
        const int first = portBlockStart[row];
        double largest = 0.0;

        for (int c = first; c < first + portBlockSize[row]; ++c)
            largest = std::max(largest, std::abs(portJacobian[row * m + static_cast<size_t>(c)]));

        if (largest * portImpedanceColumnMax[row] >= activeRowThreshold)
            activeRows[static_cast<size_t>(k++)] = r;
    }

    // Past about 85% of the rows, the k x k solve and the products around it
    // cost more than factorising J outright (k^3/3 + 2mk against m^3/3), so the
    // full step is taken instead. A transistor circuit, whose junctions all
    // conduct, never gets here.
    if (k * 20 > portCount * 17)
        return 0;

    lastActiveRowCount += k;

    if (k == 0)
        return 1; // J is the identity: the residual is the step

    // S = I + D_a K_a and y = D_a r. Row i of D_a has entries only on its own
    // device's two or three ports.
    const auto kk = static_cast<size_t>(k);
    double* system = activeSystem.data();

    for (size_t i = 0; i < kk; ++i)
    {
        const auto row = static_cast<size_t>(activeRows[i]);
        const int first = portBlockStart[row];
        double* systemRow = system + i * kk;
        double y = 0.0;

        std::fill(systemRow, systemRow + kk, 0.0);

        for (int c = first; c < first + portBlockSize[row]; ++c)
        {
            const double slope = portJacobian[row * m + static_cast<size_t>(c)];
            if (slope == 0.0)
                continue;

            y += slope * portResidual[static_cast<size_t>(c)];

            const double* impedanceRow = portImpedance.data() + static_cast<size_t>(c) * m;
            for (size_t j = 0; j < kk; ++j)
                systemRow[j] += slope * impedanceRow[static_cast<size_t>(activeRows[j])];
        }

        systemRow[i] += 1.0;
        activeRhs[i] = y;
    }

    if (! CelineEngine::DenseSolver::eliminateAndSolve(system, activeRhs.data(), k))
        return -1;

    // step = r - K_a z, a column of K at a time.
    for (size_t j = 0; j < kk; ++j)
    {
        const double z = activeRhs[j];
        if (z == 0.0)
            continue;

        const double* column = portImpedanceByColumn.data() + static_cast<size_t>(activeRows[j]) * m;
        for (size_t p = 0; p < m; ++p)
            portResidual[p] -= column[p] * z;
    }

    return 1;
}

bool Circuit::newtonOnPorts(int maxIterations, int& iterations) noexcept
{
    // Convergence is tested at the top of the loop rather than the bottom, so
    // that when it exits the port currents were evaluated at the very port
    // voltages it's exiting with. Testing at the bottom would leave the two one
    // step out of step, and the recovery in solveWithDk() would then push
    // slightly stale currents back into the circuit.
    int iteration = 0;
    bool converged = false;
    bool stepWasSmall = false;

    while (iteration < maxIterations)
    {
        ++iteration;

        const bool limited = limitPortVoltages();
        linearisePorts();

        if (stepWasSmall && ! limited)
        {
            converged = true;
            break;
        }

        // F(v) = v - vOpen + K i(v),  J = I + K di/dv
        //
        // K i a column at a time: each port's current is spread over every
        // residual at once, rather than each residual summing all the currents
        // in turn -- the same sums in the same order, so the same answer to the
        // bit, but without every multiply-add waiting on the one before it.
        for (int p = 0; p < portCount; ++p)
            portResidual[static_cast<size_t>(p)] = portOpenVoltage[static_cast<size_t>(p)] - portVoltage[static_cast<size_t>(p)];

        for (int q = 0; q < portCount; ++q)
        {
            const double current = portCurrent[static_cast<size_t>(q)];
            if (current == 0.0)
                continue;

            const double* column = portImpedanceByColumn.data() + static_cast<size_t>(q * portCount);
            for (int p = 0; p < portCount; ++p)
                portResidual[static_cast<size_t>(p)] -= column[p] * current;
        }

        // Solve J step = residual: through the current rows that matter when
        // few enough do, otherwise through the whole of J. Either way it is
        // Newton's step to the same root.
        const int reduced = reduceNewtonStep ? solveReducedNewtonStep() : 0;

        if (reduced < 0)
        {
            iterations += iteration;
            return false;
        }

        if (reduced == 0)
        {
            // J = I + K * di/dv, a column at a time. di/dv is block diagonal, so
            // column q only draws on K's columns for the ports inside q's own
            // device -- two or three of them -- and each is a run down every row
            // at once. Summed in the same order the row-at-a-time form used, so
            // the same answer to the bit.
            double* jacobian = portSolver.data();
            double* column = portJacobianColumn.data();
            const auto m = static_cast<size_t>(portCount);

            for (const auto& [firstPort, blockSize] : deviceBlocks)
            {
                for (int q = firstPort; q < firstPort + blockSize; ++q)
                {
                    std::fill(column, column + m, 0.0);

                    for (int r = firstPort; r < firstPort + blockSize; ++r)
                    {
                        const double slope = portJacobian[static_cast<size_t>(r) * m + static_cast<size_t>(q)];
                        if (slope == 0.0)
                            continue;

                        const double* impedance = portImpedanceByColumn.data() + static_cast<size_t>(r) * m;
                        for (size_t p = 0; p < m; ++p)
                            column[p] += impedance[p] * slope;
                    }

                    for (size_t p = 0; p < m; ++p)
                        jacobian[p * m + static_cast<size_t>(q)] = column[p];

                    jacobian[static_cast<size_t>(q) * m + static_cast<size_t>(q)] += 1.0;
                }
            }

            // Built for this one solve and rebuilt next iteration, so factorised
            // and solved in one pass.
            if (! portSolver.factoriseAndSolve(portResidual.data()))
            {
                iterations += iteration;
                return false;
            }
        }

        stepWasSmall = true;
        for (int p = 0; p < portCount; ++p)
        {
            const double step = portResidual[static_cast<size_t>(p)];
            const double v = portVoltage[static_cast<size_t>(p)] + step;

            if (std::abs(step) > convergenceAbsTolerance + convergenceRelTolerance * std::abs(v))
                stepWasSmall = false;

            portVoltage[static_cast<size_t>(p)] = v;
        }
    }

    iterations += iteration;
    return converged;
}

bool Circuit::solveWithDk(double vIn) noexcept
{
    const auto n = static_cast<size_t>(systemSize);

    // Step one: solve the circuit as if the nonlinear devices weren't there.
    std::copy(rhs.begin(), rhs.end(), linearSolution.begin());
    solver.solveInPlace(linearSolution.data());

    for (int p = 0; p < portCount; ++p)
    {
        const int ra = portRowA[static_cast<size_t>(p)];
        const int rb = portRowB[static_cast<size_t>(p)];

        double v = portInputCoupling[static_cast<size_t>(p)] * vIn;
        if (ra >= 0)
            v += linearSolution[static_cast<size_t>(ra)];
        if (rb >= 0)
            v -= linearSolution[static_cast<size_t>(rb)];

        portOpenVoltage[static_cast<size_t>(p)] = v;
    }

    // Where this sample starts from: the last sample's answer, which is also
    // where continuation below walks from if Newton cannot get here directly.
    std::copy(portVoltage.begin(), portVoltage.end(), portVoltageStart.begin());

    // Step two: Newton on the ports alone. v carries over from the previous
    // sample, which is why this usually takes two or three passes.
    int iterations = 0;
    lastActiveRowCount = 0;
    bool converged = newtonOnPorts(firstAttemptIterations, iterations);

    // Continuation. Between one sample and the next, all that changes in the
    // port equations is vOpen -- K and the devices hold -- and the sample
    // began from the answer to the last one. So when Newton cannot get from
    // there to here in one go, it is walked there: back to where the sample
    // started, then towards this sample's vOpen in sub-steps, each starting
    // close to its own answer. The last sub-step is this sample's own
    // equation, so the answer is the one a plain solve converging would have
    // found; only the route differs.
    //
    // A valve driven into grid current turns over -- its plate current's slope
    // against the grid changes sign within a fraction of a volt -- and two such
    // stages coupled through a 100k-odd plate load (a gain stage DC-coupled to
    // a cathode follower, as in most high-gain preamps) leave Newton a
    // straight-line model good for a few tens of millivolts. And where a grid
    // crosses the onset of grid current its conductance falls from tens of
    // kilohms to nothing within millivolts, which with a cascade's gain behind
    // it Newton orbits rather than settles on. On a Diezel and a
    // Dual Rectifier preamp it cycled or crawled to the iteration cap on
    // samples a guitar plays constantly, and each one it gave up on came out
    // as a click. With this, none do. Damping the step instead was measured
    // too: it made every sample on those sheets pay to save the few.
    //
    // The starting point's own vOpen is recomputed rather than remembered: the
    // DK equation v = vOpen - K i(v), read backwards, gives the vOpen the start
    // solves exactly, even if the sample before had itself given up.
    if (! converged)
    {
        const auto m = static_cast<size_t>(portCount);

        std::copy(portOpenVoltage.begin(), portOpenVoltage.end(), portOpenTarget.begin());
        std::copy(portVoltageStart.begin(), portVoltageStart.end(), portVoltage.begin());
        linearisePorts();

        for (size_t p = 0; p < m; ++p)
        {
            double open = portVoltageStart[p];
            for (size_t q = 0; q < m; ++q)
                open += portImpedance[p * m + q] * portCurrent[q];
            portOpenStart[p] = open;
        }

        // Adaptive: a sub-step that converges is kept and the next one is
        // twice as long; one that does not is retried from the last that did,
        // half as long. The small steps go only where the trouble is. No step
        // runs past the end of the sample, so a last step that fails is tried
        // again shorter -- not again from the same point, which can only fail
        // the same way.
        settleLimiters(portVoltageStart);
        std::copy(portVoltageStart.begin(), portVoltageStart.end(), portVoltageReached.begin());

        const int budget = iterations + continuationBudget;
        double reached = 0.0;
        double stride = 0.25;

        while (reached < 1.0 && stride >= minContinuationStride && iterations < budget)
        {
            const bool lastStep = stride >= 1.0 - reached;
            if (lastStep)
                stride = 1.0 - reached;

            const double t = lastStep ? 1.0 : reached + stride;
            for (size_t p = 0; p < m; ++p)
                portOpenVoltage[p] = portOpenStart[p] + t * (portOpenTarget[p] - portOpenStart[p]);

            if (newtonOnPorts(continuationIterations, iterations))
            {
                reached = t;
                stride *= 2.0;
                std::copy(portVoltage.begin(), portVoltage.end(), portVoltageReached.begin());
            }
            else
            {
                stride *= 0.5;
                settleLimiters(portVoltageReached);
            }
        }

        converged = reached >= 1.0;
        std::copy(portOpenTarget.begin(), portOpenTarget.end(), portOpenVoltage.begin());
    }

    lastIterationCount = iterations;

    // Step three: put the converged device currents back into the circuit.
    // z = zLinear - W i.
    std::copy(linearSolution.begin(), linearSolution.end(), work.begin());

    for (int p = 0; p < portCount; ++p)
    {
        const double current = portCurrent[static_cast<size_t>(p)];
        if (current == 0.0)
            continue;

        const double* column = portTransfer.data() + static_cast<size_t>(p) * n;
        for (size_t row = 0; row < n; ++row)
            work[row] -= column[row] * current;
    }

    writeBackSolution(false);
    return converged;
}

bool Circuit::runFullNewton(const double* baseMatrix, const double* baseRhs) noexcept
{
    if (systemSize == 0)
    {
        lastIterationCount = 0;
        return true;
    }

    const auto matrixElements = static_cast<size_t>(systemSize) * static_cast<size_t>(systemSize);
    const auto vectorElements = static_cast<size_t>(systemSize);

    if (! isNonlinear())
    {
        std::copy_n(baseMatrix, matrixElements, newtonSolver.data());
        std::copy_n(baseRhs, vectorElements, work.data());

        // Solved once and rebuilt next time: no factors kept for later solves.
        if (! newtonSolver.factoriseAndSolve(work.data()))
            return false;

        writeBackSolution();
        lastIterationCount = 1;
        return true;
    }

    int iteration = 0;
    bool converged = false;

    while (iteration < maxNewtonIterations && ! converged)
    {
        ++iteration;

        std::copy_n(baseMatrix, matrixElements, newtonSolver.data());
        std::copy_n(baseRhs, vectorElements, work.data());

        readPortVoltagesFromNodes();
        const bool limited = limitPortVoltages();
        linearisePorts();
        stampPorts(newtonSolver.data(), work.data());

        // A Jacobian is used for one solve, so the factors solveInPlace() would
        // need are never built.
        if (! newtonSolver.factoriseAndSolve(work.data()))
            return false;

        // A step the limiter had to damp isn't a converged one, however small the
        // resulting voltage change looked.
        converged = writeBackSolution() && ! limited;
    }

    lastIterationCount = iteration;
    return converged;
}

void Circuit::solveOperatingPoint() noexcept
{
    operatingPointConverged = true;

    if (systemSize == 0)
        return;

    // Nothing moves at the bias point, so no stored charge is moving either:
    // the diodes solve as their static curves, and have their charge seeded
    // from the answer below.
    for (auto& d : diodes)
        d.transitScale = 0.0;

    // No signal at the input while we find the resting state.
    nodeVoltage[static_cast<size_t>(groundIndex)] = 0.0;
    if (inputIndex >= 0)
        nodeVoltage[static_cast<size_t>(inputIndex)] = 0.0;

    buildDcSystem(1.0);

    if (! runFullNewton(dcMatrix.data(), dcRhs.data()))
    {
        // Source stepping. Newton needs a starting guess in the right
        // neighbourhood, and for a biased transistor stage "everything at 0 V"
        // isn't: the supply is the whole reason any current flows. So bring the
        // supplies up from zero in steps, each solve warm-starting from the last.
        // At zero supply the answer is trivially zero volts everywhere, and each
        // step moves the operating point a little, so Newton always starts close.
        clearState();
        nodeVoltage[static_cast<size_t>(groundIndex)] = 0.0;
        if (inputIndex >= 0)
            nodeVoltage[static_cast<size_t>(inputIndex)] = 0.0;

        for (int step = 1; step <= sourceSteppingSteps; ++step)
        {
            buildDcSystem(static_cast<double>(step) / sourceSteppingSteps);
            operatingPointConverged = runFullNewton(dcMatrix.data(), dcRhs.data());
        }

        jassert(operatingPointConverged); // bias network doesn't resolve -- check the wiring
    }

    // Seed the companion models from the operating point: a capacitor holds its
    // DC voltage and passes no current, an inductor drops no voltage and carries
    // whatever current the short was passing.
    reversedCapacitorCount = 0;

    for (auto& c : capacitors)
    {
        c.vCapPrev = nodeVoltage[static_cast<size_t>(c.a)] - nodeVoltage[static_cast<size_t>(c.b)];
        c.iPrev = 0.0; // no current at DC, so no ESR drop to separate out

        // An electrolytic sitting backwards is a build error, not a sound.
        if (c.polarised && c.vCapPrev < -reverseBiasTolerance)
            ++reversedCapacitorCount;
    }

    jassert(reversedCapacitorCount == 0); // a polarised capacitor is wired backwards

    for (auto& l : inductors)
    {
        const double v = nodeVoltage[static_cast<size_t>(l.a)] - nodeVoltage[static_cast<size_t>(l.b)];
        l.vPrev = 0.0;
        l.iPrev = dcInductorConductance * v;
    }

    // Seed the port voltages too, so the first audio sample's Newton starts from
    // the operating point rather than from nothing.
    readPortVoltagesFromNodes();

    // And the diodes' stored charge: whatever they carry at rest is charge they
    // already hold, not charge the first sample has to supply.
    for (auto& d : diodes)
    {
        d.transitScale = d.model.transitTime / dt;

        if (d.transitScale > 0.0)
            d.chargeCurrentPrevious = CircuitComponents::chargeCurrent(
                d, nodeVoltage[static_cast<size_t>(d.anode)] - nodeVoltage[static_cast<size_t>(d.cathode)]);
    }
}

//==============================================================================
// Per-sample processing
//==============================================================================

void Circuit::buildRightHandSide(double vIn) noexcept
{
    // Nothing in the cached system contributes a constant term: a source's
    // voltage is added below so it can vary, and the only other known node is
    // ground at 0 V. So the vector starts empty rather than from a stored one,
    // and the handful of rows vIn reaches are added individually.
    std::fill(rhs.begin(), rhs.end(), 0.0);

    double* rhsVector = rhs.data();

    for (const auto& [row, coefficient] : inputCouplingRows)
        rhsVector[row] += coefficient * vIn;

    for (size_t s = 0; s < voltageSources.size(); ++s)
        rhsVector[static_cast<size_t>(numNodeUnknowns) + s] += sourceVoltage(voltageSources[s]);

    for (const auto& c : capacitors)
        stampCurrentSource(rhsVector, c.a, c.b,
                           c.conductance * (c.vCapPrev + c.halfStepResistance * c.iPrev));

    for (const auto& l : inductors)
        stampCurrentSource(rhsVector, l.a, l.b, -(l.iPrev + l.conductance * l.vPrev));
}

void Circuit::updateReactiveState() noexcept
{
    for (auto& source : voltageSources)
    {
        if (source.acFrequency > 0.0)
        {
            source.phase += 2.0 * juce::MathConstants<double>::pi * source.acFrequency * dt;

            if (source.phase > 2.0 * juce::MathConstants<double>::pi)
                source.phase -= 2.0 * juce::MathConstants<double>::pi;
        }
    }

    for (auto& c : capacitors)
    {
        const double v = nodeVoltage[static_cast<size_t>(c.a)] - nodeVoltage[static_cast<size_t>(c.b)];
        const double vEquivalent = c.vCapPrev + c.halfStepResistance * c.iPrev;

        c.iPrev = c.conductance * (v - vEquivalent);

        // The capacitance itself only sees what's left after the ESR drop.
        c.vCapPrev = v - c.iPrev * c.esrOhms;
    }

    for (auto& l : inductors)
    {
        const double v = nodeVoltage[static_cast<size_t>(l.a)] - nodeVoltage[static_cast<size_t>(l.b)];
        const double ieq = -(l.iPrev + l.conductance * l.vPrev);

        l.iPrev = l.conductance * v - ieq;
        l.vPrev = v;
    }

    // The charge a slow diode holds at the end of this sample is what the next
    // one starts from. One exponential per such diode; the rest skip it.
    for (auto& d : diodes)
        if (d.transitScale > 0.0)
            d.chargeCurrentPrevious = CircuitComponents::chargeCurrent(
                d, nodeVoltage[static_cast<size_t>(d.anode)] - nodeVoltage[static_cast<size_t>(d.cathode)]);
}

float Circuit::process(float vIn)
{
    jassert(inputIndex >= 0 && outputIndex >= 0);    // call setInputNode()/setOutputNode() first
    jassert(rowOfNode.size() == nodeIndices.size()); // netlist grew since prepare() -- call prepare() again

    // Only on an actual value or topology change. Retrying on a failed
    // factorisation instead -- which is what `|| ! factorisationValid` did --
    // rebuilds an identical matrix and fails identically, once per sample: a
    // full O(n^3) stamp and factorise, 48000 times a second, all of it doomed.
    // The guard below already returns silence for that case, so the retry never
    // bought anything. linearDirty starts true and every setter re-sets it, so
    // no legitimate rebuild is skipped by this.
    if (linearDirty)
        rebuildLinearSystem();

    if (systemSize == 0 || ! factorisationValid)
        return 0.0f;

    nodeVoltage[static_cast<size_t>(inputIndex)] = vIn;
    buildRightHandSide(vIn);

    if (portCount == 0)
    {
        // Purely linear: the factorisation is still good, so this is just the
        // two triangular substitutions. Nothing iterates, so there is no
        // convergence to test -- asking for it would cost a per-node abs, max
        // and compare on every sample for an answer that is discarded.
        std::copy(rhs.begin(), rhs.end(), work.begin());
        solver.solveInPlace(work.data());
        writeBackSolution(false);
        lastIterationCount = 1;
    }
    else if (activeStrategy == SolverStrategy::DiscreteK)
    {
        if (! solveWithDk(vIn))
            ++nonConvergenceCount;
    }
    else if (! runFullNewton(linearMatrix.data(), rhs.data()))
    {
        ++nonConvergenceCount;
    }

    updateReactiveState();

    const auto out = static_cast<float>(nodeVoltage[static_cast<size_t>(outputIndex)] - outputOffset);

    // A non-finite sample means the solve came apart; drop the state rather than
    // let a NaN propagate into the host for the rest of the session.
    if (! std::isfinite(out))
    {
        clearState();
        return 0.0f;
    }

    return out;
}

void Circuit::process(float* samples, int numSamples)
{
    for (int i = 0; i < numSamples; ++i)
        samples[i] = process(samples[i]);
}
