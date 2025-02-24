# This code is part of Qiskit.
#
# (C) Copyright IBM 2025.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Library of standard-gate equivalances."""

from math import pi

# "std" is for "standard gates", since we use "lib" to mean the equivalence library.
from qiskit.circuit import EquivalenceLibrary, Parameter, QuantumCircuit, library as std


def standard_equivalence_library() -> EquivalenceLibrary:
    """Construct a new standard-equivalence library.

    Maintainer notes
    ================

    Purpose of equivalences
    -----------------------

    The rules in the standard equivalence library are primarily used to power the basis translator,
    as part of the compilation pipeline.  We do not need to have equivalences between every suitable
    pair of gates (or bases, etc); the translator can combine multiple equivalences into a single
    rule.

    There are a few reasons to add equivalences:

    1. We should be able to translate any (over-)complete parametric bases of standard gates to any
       other, even if the parameters are completely opaque to us (i.e. numerical synthesis cannot
       work).

    2. We should be able to translate any (over-)complete discrete basis of standard into any other,
       without relying on numerical synthesis (e.g. there should be paths to translate sequences of
       Clifford operations into bases that contain only Clifford gates).

    3. If there is a special case that numerical synthesis might struggle to find (especially with
       the higher-degree standard gates, or with parametric 2q+ gates), it is worthwhile adding
       equivalences.

    Remember that basis translation is just one part of the translation pipeline, and while it
    should be sufficient to produce a _valid_ translation (up to gate direction), it need not worry
    about _optimisation_ of those sequences, especially 1q sequences.  Our optimisation and
    resynthesis routines will handle those.

    There is generally no need to have 1q equivalence relations more than a couple of standard gates
    long.  The basis translator can construct arbitrarily long sequences to build up a gate, and the
    simpler the translations are, the more likely we will have a path to the most efficient
    translation for any particular target basis.  The more qubits, the longer the sensible set of
    equivalence relations is likely to be.


    Organisation
    ------------

    This file is organised into separate functions, that each handle one group of cases.  The cases
    are grouped to target a specific component of the reasons in the above section; the use of
    functions is just a grouping construct, and not exposed to users.  This function then calls all
    the others to build a standard equivalence library, with the intention that you can reason about
    the top level "structure" of the equivalence library without needing to know all the individual
    rules going into it.

    You can add more groups if needed.  Try to group things into a coherent shared purposes, and
    write a test that enforces several different cases work within that shared goal, beyond just the
    minimal set of equivalences you added.

    It's not a problem if there's a little bit of duplication between groups if it helps keep things
    looking more structured and understandable (though try not to go crazy)."""
    out = EquivalenceLibrary()
    hierarchical_definitions(out)
    # 1q relations.
    within_discrete_1q(out)
    discrete_to_pauli_rotation_1q(out)
    within_single_parameter_rotation_1q(out)
    single_to_two_parameter_rotation_1q(out)
    single_u_1q(out)
    single_multi_parameter_to_multi_single_parameter_1q(out)
    # 2q relations.
    local_equivalence_cx(out)
    local_equivalence_csx(out)
    local_equivalence_iswap(out)
    local_equivalence_ising(out)
    reverse_gate_direction_2q(out)
    efficient_swap(out)
    between_locally_equivalent_groups_2q(out)
    # Higher order relations.
    multi_q(out)
    return out


def hierarchical_definitions(lib: EquivalenceLibrary):
    """Add each standard gate's hierarchical "definition" field to the equivalence library."""
    alpha, beta, gamma, delta = (Parameter("α"), Parameter("β"), Parameter("γ"), Parameter("δ"))
    for gate in (
        std.HGate(),
        std.IGate(),
        std.XGate(),
        std.YGate(),
        std.ZGate(),
        std.SGate(),
        std.SdgGate(),
        std.SXGate(),
        std.SXdgGate(),
        std.TGate(),
        std.TdgGate(),
        std.CHGate(),
        std.CXGate(),
        std.CYGate(),
        std.CZGate(),
        std.DCXGate(),
        std.ECRGate(),
        std.SwapGate(),
        std.iSwapGate(),
        std.CSGate(),
        std.CSdgGate(),
        std.CSXGate(),
        std.CCXGate(),
        std.CCZGate(),
        std.CSwapGate(),
        std.RCCXGate(),
        # std.C3XGate(),  # skipped because its name 'mcx' is way overloaded and messes things up.
        std.C3SXGate(),
        std.RC3XGate(),
        std.GlobalPhaseGate(alpha),
        std.PhaseGate(alpha),
        std.RXGate(alpha),
        std.RYGate(alpha),
        std.RZGate(alpha),
        std.U1Gate(alpha),
        std.CPhaseGate(alpha),
        std.CRXGate(alpha),
        std.CRYGate(alpha),
        std.CRZGate(alpha),
        std.CU1Gate(alpha),
        std.RXXGate(alpha),
        std.RYYGate(alpha),
        std.RZZGate(alpha),
        std.RZXGate(alpha),
        std.RGate(alpha, beta),
        std.U2Gate(alpha, beta),
        std.XXMinusYYGate(alpha, beta),
        std.XXPlusYYGate(alpha, beta),
        std.UGate(alpha, beta, gamma),
        std.U3Gate(alpha, beta, gamma),
        std.CU3Gate(alpha, beta, gamma),
        std.CUGate(alpha, beta, gamma, delta),
    ):
        if gate.definition is not None:
            lib.add_equivalence(gate, gate.definition)


def within_discrete_1q(lib: EquivalenceLibrary):
    r"""Add equivalences for the discrete 1q gates in the standard library.

    There is a hierarchy of "powers" of gates in our 1q discrete set, and users may have interest in
    having the transpiler not spuriously introduce gates of a higher "power".  We categorise these
    from "weakest" to "strongest" as:

    * Paulis
    * Cliffords
    * Clifford+T

    We arrange the rules so it's always possible to translate any generating set of one group to any
    generating set of the same group or any more powerful group."""
    # Intra-Pauli-gate conversions.

    x_pauli = QuantumCircuit(1, global_phase=pi / 2)
    x_pauli.y(0)
    x_pauli.z(0)
    lib.add_equivalence(std.XGate(), x_pauli)

    y_pauli = QuantumCircuit(1, global_phase=pi / 2)
    y_pauli.z(0)
    y_pauli.x(0)
    lib.add_equivalence(std.YGate(), y_pauli)

    z_pauli = QuantumCircuit(1, global_phase=pi / 2)
    z_pauli.x(0)
    z_pauli.y(0)
    lib.add_equivalence(std.ZGate(), z_pauli)

    # Now bring H into the mix for the Paulis.  The previous rules take care of the rest of
    # converting Y to either X+H or Z+H.

    x_to_z = QuantumCircuit(1)
    x_to_z.h(0)
    x_to_z.z(0)
    x_to_z.h(0)
    lib.add_equivalence(std.XGate(), x_to_z)

    z_to_x = QuantumCircuit(1)
    z_to_x.h(0)
    z_to_x.x(0)
    z_to_x.h(0)
    lib.add_equivalence(std.ZGate(), z_to_x)

    # Similarly, X and Z in terms of their square roots.

    x_to_sx = QuantumCircuit(1)
    x_to_sx.sx(0)
    x_to_sx.sx(0)
    lib.add_equivalence(std.XGate(), x_to_sx)

    x_to_sxdg = QuantumCircuit(1)
    x_to_sxdg.sxdg(0)
    x_to_sxdg.sxdg(0)
    lib.add_equivalence(std.XGate(), x_to_sxdg)

    z_to_s = QuantumCircuit(1)
    z_to_s.s(0)
    z_to_s.s(0)
    lib.add_equivalence(std.ZGate(), z_to_s)

    z_to_sdg = QuantumCircuit(1)
    z_to_sdg.sdg(0)
    z_to_sdg.sdg(0)
    lib.add_equivalence(std.ZGate(), z_to_sdg)

    # The non-Pauli Clifford gates in terms of each other, so the minimal relations are reachable.

    h_to_sqrt = QuantumCircuit(1, global_phase=-pi / 4)
    h_to_sqrt.s(0)
    h_to_sqrt.sx(0)
    h_to_sqrt.s(0)
    lib.add_equivalence(std.HGate(), h_to_sqrt)

    h_to_sqrtdg = QuantumCircuit(1, global_phase=pi / 4)
    h_to_sqrtdg.sdg(0)
    h_to_sqrtdg.sxdg(0)
    h_to_sqrtdg.sdg(0)
    lib.add_equivalence(std.HGate(), h_to_sqrtdg)

    s_to_sdg = QuantumCircuit(1)
    s_to_sdg.sdg(0)
    s_to_sdg.z(0)
    lib.add_equivalence(std.SGate(), s_to_sdg)

    s_to_sx = QuantumCircuit(1)
    s_to_sx.h(0)
    s_to_sx.sx(0)
    s_to_sx.h(0)
    lib.add_equivalence(std.SGate(), s_to_sx)

    sdg_to_s = QuantumCircuit(1)
    sdg_to_s.s(0)
    sdg_to_s.z(0)
    lib.add_equivalence(std.SdgGate(), sdg_to_s)

    sx_to_s = QuantumCircuit(1)
    sx_to_s.h(0)
    sx_to_s.s(0)
    sx_to_s.h(0)
    lib.add_equivalence(std.SXGate(), sx_to_s)

    sxdg_to_s = QuantumCircuit(1, global_phase=-pi / 4)
    sxdg_to_s.s(0)
    sxdg_to_s.h(0)
    sxdg_to_s.s(0)
    lib.add_equivalence(std.SXdgGate(), sxdg_to_s)

    # ... and now relations that bring `t` and `tdg` in.

    s_to_t = QuantumCircuit(1)
    s_to_t.t(0)
    s_to_t.t(0)
    lib.add_equivalence(std.SGate(), s_to_t)

    sdg_to_tdg = QuantumCircuit(1)
    sdg_to_tdg.tdg(0)
    sdg_to_tdg.tdg(0)
    lib.add_equivalence(std.SdgGate(), sdg_to_tdg)

    t_to_tdg = QuantumCircuit(1)
    t_to_tdg.tdg(0)
    t_to_tdg.sdg(0)
    t_to_tdg.z(0)
    lib.add_equivalence(std.TGate(), t_to_tdg)

    tdg_to_t = QuantumCircuit(1)
    tdg_to_t.t(0)
    tdg_to_t.s(0)
    tdg_to_t.z(0)
    lib.add_equivalence(std.TdgGate(), tdg_to_t)


def discrete_to_pauli_rotation_1q(lib: EquivalenceLibrary):
    """Raise all discrete standard gates to a form that uses discrete gates and a 1q Pauli
    rotation.

    This lets us give typically short rules to target the entire 1q rotation space, since we
    separately have conversions between the 1q rotations."""
    # This function is ordered the same as the Rust-space `StandardGate` enum.

    h_to_ry_x = QuantumCircuit(1)
    h_to_ry_x.ry(pi / 2, 0)
    h_to_ry_x.x(0)
    lib.add_equivalence(std.HGate(), h_to_ry_x)

    h_to_ry_z = QuantumCircuit(1)
    h_to_ry_z.z(0)
    h_to_ry_z.ry(pi / 2, 0)
    lib.add_equivalence(std.HGate(), h_to_ry_z)

    x_to_rx = QuantumCircuit(1, global_phase=pi / 2)
    x_to_rx.rx(pi, 0)
    lib.add_equivalence(std.XGate(), x_to_rx)

    y_to_ry = QuantumCircuit(1, global_phase=pi / 2)
    y_to_ry.ry(pi, 0)
    lib.add_equivalence(std.YGate(), y_to_ry)

    z_to_rz = QuantumCircuit(1, global_phase=pi / 2)
    z_to_rz.rz(pi, 0)
    lib.add_equivalence(std.ZGate(), z_to_rz)

    s_to_rz = QuantumCircuit(1, global_phase=pi / 4)
    s_to_rz.rz(pi / 2, 0)
    lib.add_equivalence(std.SGate(), s_to_rz)

    sdg_to_rz = QuantumCircuit(1, global_phase=-pi / 4)
    sdg_to_rz.rz(-pi / 2, 0)
    lib.add_equivalence(std.SdgGate(), sdg_to_rz)

    sx_to_rx = QuantumCircuit(1, global_phase=pi / 4)
    sx_to_rx.rx(pi / 2, 0)
    lib.add_equivalence(std.SXGate(), sx_to_rx)

    sxdg_to_rx = QuantumCircuit(1, global_phase=-pi / 4)
    sxdg_to_rx.rx(-pi / 2, 0)
    lib.add_equivalence(std.SXdgGate(), sxdg_to_rx)

    t_to_rz = QuantumCircuit(1, global_phase=pi / 8)
    t_to_rz.rz(pi / 4, 0)
    lib.add_equivalence(std.TGate(), t_to_rz)

    tdg_to_rz = QuantumCircuit(1, global_phase=-pi / 8)
    tdg_to_rz.rz(-pi / 4, 0)
    lib.add_equivalence(std.TdgGate(), tdg_to_rz)


def within_single_parameter_rotation_1q(lib: EquivalenceLibrary):
    """Add rules that convert between different 1q rotations.

    In general, we prefer to use the "weakest" set of additional gates needed in the equivalence (as
    described in the discrete 1q rules), so that the other parts of the rules can handle conversion
    to more powerful basis sets."""
    alpha = Parameter("α")

    # First, mark `p` and `u1` as exactly equivalent, then we can stop dealing with `u1` as the
    # `p` rules will always take care of it.  Similarly, `p` and `rz` are equivalent up to phase, so
    # we can deal only with `rz`.

    phase_to_u1 = QuantumCircuit(1)
    phase_to_u1.append(std.U1Gate(alpha), [0], [])
    lib.add_equivalence(std.PhaseGate(alpha), phase_to_u1)

    u1_to_phase = QuantumCircuit(1)
    u1_to_phase.p(alpha, 0)
    lib.add_equivalence(std.U1Gate(alpha), u1_to_phase)

    p_to_rz = QuantumCircuit(1, global_phase=alpha / 2)
    p_to_rz.rz(alpha, 0)
    lib.add_equivalence(std.PhaseGate(alpha), p_to_rz)

    rz_to_p = QuantumCircuit(1, global_phase=-alpha / 2)
    rz_to_p.p(alpha, 0)
    lib.add_equivalence(std.RZGate(alpha), rz_to_p)

    # Now we want to link all the Pauli rotations to each other.  Ideally we use the "weakest" set
    # of additional gates necessary for the transformation, as the other discrete rules will take
    # care of mapping these up to more powerful gates (like fixed-angle rotations).

    rx_to_ry = QuantumCircuit(1)
    rx_to_ry.s(0)
    rx_to_ry.ry(alpha, 0)
    rx_to_ry.sdg(0)
    lib.add_equivalence(std.RXGate(alpha), rx_to_ry)

    rx_to_rz = QuantumCircuit(1)
    rx_to_rz.h(0)
    rx_to_rz.rz(alpha, 0)
    rx_to_rz.h(0)
    lib.add_equivalence(std.RXGate(alpha), rx_to_rz)

    ry_to_rx = QuantumCircuit(1)
    ry_to_rx.sdg(0)
    ry_to_rx.rx(alpha, 0)
    ry_to_rx.s(0)
    lib.add_equivalence(std.RYGate(alpha), ry_to_rx)

    ry_to_rz = QuantumCircuit(1)
    ry_to_rz.sx(0)
    ry_to_rz.rz(alpha, 0)
    ry_to_rz.sxdg(0)
    lib.add_equivalence(std.RYGate(alpha), ry_to_rz)

    rz_to_rx = QuantumCircuit(1)
    rz_to_rx.h(0)
    rz_to_rx.rx(alpha, 0)
    rz_to_rx.h(0)
    lib.add_equivalence(std.RZGate(alpha), rz_to_rx)

    rz_to_ry = QuantumCircuit(1)
    rz_to_ry.sxdg(0)
    rz_to_ry.ry(alpha, 0)
    rz_to_ry.sx(0)
    lib.add_equivalence(std.RZGate(alpha), rz_to_ry)


def single_to_two_parameter_rotation_1q(lib: EquivalenceLibrary):
    """Equivalence rules to take single-parameter rotations to two-parameter rotations.

    The rest of the multi-parameter rotations are handled by raising the two-parameter forms to the
    three-parameter gates (only `u3` and `u`, which are equivalent)."""
    # This is a bit of an empty function, because `u2` is pretty useless for describing our
    # one-parameter rotations.  It's a historical artefact from when IBM's stated basis set was
    # `[u1, u2, u3, cx]`; `u2` gave a way of describing `h` with only a single implied `sx`.  We
    # throw in a couple of those old relations just to maintain that historical component.
    #
    # The only meaningful two-parameter gate we have is `r`.

    alpha = Parameter("α")

    h_to_u2 = QuantumCircuit(1)
    h_to_u2.append(std.U2Gate(0, pi), [0])
    lib.add_equivalence(std.HGate(), h_to_u2)

    sx_to_u2 = QuantumCircuit(1, global_phase=pi / 4)
    sx_to_u2.append(std.U2Gate(-pi / 2, pi / 2), [0])
    lib.add_equivalence(std.SXGate(), sx_to_u2)

    sxdg_to_u2 = QuantumCircuit(1, global_phase=-pi / 4)
    sxdg_to_u2.append(std.U2Gate(pi / 2, -pi / 2), [0])
    lib.add_equivalence(std.SXdgGate(), sxdg_to_u2)

    # Actual single-parameter to two-parameter conversion.

    rx_to_r = QuantumCircuit(1)
    rx_to_r.r(alpha, 0, 0)
    lib.add_equivalence(std.RXGate(alpha), rx_to_r)

    ry_to_r = QuantumCircuit(1)
    ry_to_r.r(alpha, pi / 2, 0)
    lib.add_equivalence(std.RYGate(alpha), ry_to_r)


def single_u_1q(lib: EquivalenceLibrary):
    """Conversions between all 1q gates and the single-`U` form."""
    alpha, beta, gamma = Parameter("α"), Parameter("β"), Parameter("γ")

    # All 1q gates already have conversions to a single gate that is a continuous rotation (even `h`
    # in terms of `u2`), so we just need to lift those rotations to `u` and we're done in the most
    # efficient translations possible.
    rz_to_u = QuantumCircuit(1, global_phase=-alpha / 2)
    rz_to_u.u(0, alpha, 0, 0)
    lib.add_equivalence(std.RZGate(alpha), rz_to_u)

    r_to_u = QuantumCircuit(1)
    r_to_u.u(alpha, beta - (pi / 2), (pi / 2) - beta, 0)
    lib.add_equivalence(std.RGate(alpha, beta), r_to_u)

    u2_to_u = QuantumCircuit(1)
    u2_to_u.u(pi / 2, alpha, beta, 0)
    lib.add_equivalence(std.U2Gate(alpha, beta), u2_to_u)

    u3_to_u = QuantumCircuit(1)
    u3_to_u.u(alpha, beta, gamma, 0)
    lib.add_equivalence(std.U3Gate(alpha, beta, gamma), u3_to_u)

    # Also mark the backwards direction for `u3`.
    u_to_u3 = QuantumCircuit(1)
    u_to_u3.append(std.U3Gate(alpha, beta, gamma), [0])
    lib.add_equivalence(std.UGate(alpha, beta, gamma), u_to_u3)


def single_multi_parameter_to_multi_single_parameter_1q(lib: EquivalenceLibrary):
    """Rules that convert multi-parameter 1q gates to multiple separate one-parameter 1q gates.

    These are useful for going "downwards" in the operation tree; it lets us basis-translate
    higher-order 1q gates back down into the 1q set.  Typically we expect _not_ to do this---our 1q
    optimisation passes should generally handle it---but this makes it easier for a user with their
    own custom 1q gate to define rules, and for us to be able to translate things like `u` into
    them."""
    alpha, beta, gamma = Parameter("α"), Parameter("β"), Parameter("γ")

    u_to_rz = QuantumCircuit(1, global_phase=(gamma / 2) + (beta / 2) + (pi / 2))
    u_to_rz.rz(gamma, 0)
    u_to_rz.sx(0)
    u_to_rz.rz(alpha + pi, 0)
    u_to_rz.sx(0)
    u_to_rz.rz(beta + pi, 0)
    lib.add_equivalence(std.UGate(alpha, beta, gamma), u_to_rz)

    u2_to_rz = QuantumCircuit(1, global_phase=(alpha / 2) + (beta / 2) + (7 * pi / 4))
    u2_to_rz.rz(beta - (pi / 2), 0)
    u2_to_rz.sx(0)
    u2_to_rz.rz(alpha + (pi / 2), 0)
    lib.add_equivalence(std.U2Gate(alpha, beta), u2_to_rz)

    # For `r` rotations: rotate the cos(beta)X+sin(beta)Y axis to line up with X or Y, rotate around
    # that, then rotate back.

    r_to_rz_rx = QuantumCircuit(1)
    r_to_rz_rx.rz(-beta, 0)
    r_to_rz_rx.rx(alpha, 0)
    r_to_rz_rx.rz(beta, 0)
    lib.add_equivalence(std.RGate(alpha, beta), r_to_rz_rx)

    r_to_rz_ry = QuantumCircuit(1)
    r_to_rz_ry.rz((pi / 2) - beta, 0)
    r_to_rz_ry.ry(alpha, 0)
    r_to_rz_ry.rz(beta - (pi / 2), 0)
    lib.add_equivalence(std.RGate(alpha, beta), r_to_rz_ry)


def local_equivalence_cx(lib: EquivalenceLibrary):
    """Relations between the gates that are locally equivalent to CX.

    This handles `cx`, `cy`, `cz`, `ch` and `ecr`."""

    # At a minimum, we define `cx` into everything and everything into `cx`; that will make it
    # possible to translate between this equivalence class using only one 2q gate, and 1q peephole
    # optimisation will take care of the rest.  `cx` is still a highly common target gate.
    #
    # We can add more specialised rules for more common cases to ease the pressure on subsequent
    # optimisation loops.

    # CX into everything.

    # q_0: ───────■─────────
    #      ┌───┐┌─┴─┐┌─────┐
    # q_1: ┤ S ├┤ Y ├┤ Sdg ├
    #      └───┘└───┘└─────┘
    cx_to_cy = QuantumCircuit(2)
    cx_to_cy.s(1)
    cx_to_cy.cy(0, 1)
    cx_to_cy.sdg(1)
    lib.add_equivalence(std.CXGate(), cx_to_cy)

    # q_0: ──────■──────
    #      ┌───┐ │ ┌───┐
    # q_1: ┤ H ├─■─┤ H ├
    #      └───┘   └───┘
    cx_to_cz = QuantumCircuit(2)
    cx_to_cz.h(1)
    cx_to_cz.cz(0, 1)
    cx_to_cz.h(1)
    lib.add_equivalence(std.CXGate(), cx_to_cz)

    # q_0: ──────────────■─────────────
    #      ┌──────────┐┌─┴─┐┌─────────┐
    # q_1: ┤ Ry(-π/4) ├┤ H ├┤ Ry(π/4) ├
    #      └──────────┘└───┘└─────────┘
    cx_to_ch_ry = QuantumCircuit(2)
    cx_to_ch_ry.ry(-pi / 4, 1)
    cx_to_ch_ry.ch(0, 1)
    cx_to_ch_ry.ry(pi / 4, 1)
    lib.add_equivalence(std.CXGate(), cx_to_ch_ry)

    # This is the same as the `ry`+`ch` rule above, but explicitly decomposed into a discrete basis.
    # q_0: ───────────────────────■─────────────────────
    #      ┌────┐┌─────┐┌──────┐┌─┴─┐┌────┐┌───┐┌──────┐
    # q_1: ┤ √X ├┤ Tdg ├┤ √Xdg ├┤ H ├┤ √X ├┤ T ├┤ √Xdg ├
    #      └────┘└─────┘└──────┘└───┘└────┘└───┘└──────┘
    cx_to_ch_discrete = QuantumCircuit(2)
    cx_to_ch_discrete.sx(1)
    cx_to_ch_discrete.tdg(1)
    cx_to_ch_discrete.sxdg(1)
    cx_to_ch_discrete.ch(0, 1)
    cx_to_ch_discrete.sx(1)
    cx_to_ch_discrete.t(1)
    cx_to_ch_discrete.sxdg(1)
    lib.add_equivalence(std.CXGate(), cx_to_ch_discrete)

    #      ┌─────┐ ┌──────┐┌───┐
    # q_0: ┤ Sdg ├─┤0     ├┤ X ├
    #      ├─────┴┐│  Ecr │└───┘
    # q_1: ┤ √Xdg ├┤1     ├─────
    #      └──────┘└──────┘
    cx_to_ecr = QuantumCircuit(2, global_phase=pi / 4)
    cx_to_ecr.sdg(0)
    cx_to_ecr.sxdg(1)
    cx_to_ecr.ecr(0, 1)
    cx_to_ecr.x(0)
    lib.add_equivalence(std.CXGate(), cx_to_ecr)

    # CY into others.

    # q_0: ─────────■───────
    #      ┌─────┐┌─┴─┐┌───┐
    # q_1: ┤ Sdg ├┤ X ├┤ S ├
    #      └─────┘└───┘└───┘
    cy_to_cx = QuantumCircuit(2)
    cy_to_cx.sdg(1)
    cy_to_cx.cx(0, 1)
    cy_to_cx.s(1)
    lib.add_equivalence(std.CYGate(), cy_to_cx)

    # q_0: ───────■─────────
    #      ┌────┐ │ ┌──────┐
    # q_1: ┤ √X ├─■─┤ √Xdg ├
    #      └────┘   └──────┘
    cy_to_cz = QuantumCircuit(2)
    cy_to_cz.sx(1)
    cy_to_cz.cz(0, 1)
    cy_to_cz.sxdg(1)
    lib.add_equivalence(std.CYGate(), cy_to_cz)

    # CZ into everything.

    # q_0: ───────■───────
    #      ┌───┐┌─┴─┐┌───┐
    # q_1: ┤ H ├┤ X ├┤ H ├
    #      └───┘└───┘└───┘
    cz_to_cx = QuantumCircuit(2)
    cz_to_cx.h(1)
    cz_to_cx.cx(0, 1)
    cz_to_cx.h(1)
    lib.add_equivalence(std.CZGate(), cz_to_cx)

    # q_0: ──────────■────────
    #      ┌──────┐┌─┴─┐┌────┐
    # q_1: ┤ √Xdg ├┤ Y ├┤ √X ├
    #      └──────┘└───┘└────┘
    cz_to_cy = QuantumCircuit(2)
    cz_to_cy.sxdg(1)
    cz_to_cy.cy(0, 1)
    cz_to_cy.sx(1)
    lib.add_equivalence(std.CZGate(), cz_to_cy)

    # q_0: ─────────────■──────────────
    #      ┌─────────┐┌─┴─┐┌──────────┐
    # q_1: ┤ Ry(π/4) ├┤ H ├┤ Ry(-π/4) ├
    #      └─────────┘└───┘└──────────┘
    cz_to_ch_ry = QuantumCircuit(2)
    cz_to_ch_ry.ry(pi / 4, 1)
    cz_to_ch_ry.ch(0, 1)
    cz_to_ch_ry.ry(-pi / 4, 1)
    lib.add_equivalence(std.CZGate(), cz_to_ch_ry)

    # This is the same as the `ry`+`ch` rule above, but explicitly decomposed into a discrete basis.
    # q_0: ─────────────────────■───────────────────────
    #      ┌────┐┌───┐┌──────┐┌─┴─┐┌────┐┌─────┐┌──────┐
    # q_1: ┤ √X ├┤ T ├┤ √Xdg ├┤ H ├┤ √X ├┤ Tdg ├┤ √Xdg ├
    #      └────┘└───┘└──────┘└───┘└────┘└─────┘└──────┘
    cz_to_ch_discrete = QuantumCircuit(2)
    cz_to_ch_discrete.sx(1)
    cz_to_ch_discrete.t(1)
    cz_to_ch_discrete.sxdg(1)
    cz_to_ch_discrete.ch(0, 1)
    cz_to_ch_discrete.sx(1)
    cz_to_ch_discrete.tdg(1)
    cz_to_ch_discrete.sxdg(1)
    lib.add_equivalence(std.CZGate(), cz_to_ch_discrete)

    #                ┌──────┐ ┌───┐  ┌───┐
    # q_0: ──────────┤0     ├─┤ S ├──┤ X ├──
    #      ┌───┐┌───┐│  Ecr │┌┴───┴┐┌┴───┴─┐
    # q_1: ┤ Z ├┤ H ├┤1     ├┤ Sdg ├┤ √Xdg ├
    #      └───┘└───┘└──────┘└─────┘└──────┘
    cz_to_ecr = QuantumCircuit(2)
    cz_to_ecr.z(1)
    cz_to_ecr.h(1)
    cz_to_ecr.ecr(0, 1)
    cz_to_ecr.s(0)
    cz_to_ecr.x(0)
    cz_to_ecr.sdg(1)
    cz_to_ecr.sxdg(1)
    lib.add_equivalence(std.CZGate(), cz_to_ecr)

    # CH into others.

    # q_0: ─────────────■──────────────
    #      ┌─────────┐┌─┴─┐┌──────────┐
    # q_1: ┤ Ry(π/4) ├┤ X ├┤ Ry(-π/4) ├
    #      └─────────┘└───┘└──────────┘
    ch_to_cx_ry = QuantumCircuit(2)
    ch_to_cx_ry.ry(pi / 4, 1)
    ch_to_cx_ry.cx(0, 1)
    ch_to_cx_ry.ry(-pi / 4, 1)
    lib.add_equivalence(std.CHGate(), ch_to_cx_ry)

    # This is the same as the `ry`+`cx` rule above, but explicitly decomposed into a discrete basis.
    # q_0: ─────────────────────■───────────────────────
    #      ┌────┐┌───┐┌──────┐┌─┴─┐┌────┐┌─────┐┌──────┐
    # q_1: ┤ √X ├┤ T ├┤ √Xdg ├┤ X ├┤ √X ├┤ Tdg ├┤ √Xdg ├
    #      └────┘└───┘└──────┘└───┘└────┘└─────┘└──────┘
    ch_to_cx_discrete = QuantumCircuit(2)
    ch_to_cx_discrete.sx(1)
    ch_to_cx_discrete.t(1)
    ch_to_cx_discrete.sxdg(1)
    ch_to_cx_discrete.cx(0, 1)
    ch_to_cx_discrete.sx(1)
    ch_to_cx_discrete.tdg(1)
    ch_to_cx_discrete.sxdg(1)
    lib.add_equivalence(std.CHGate(), ch_to_cx_discrete)

    # q_0: ─────────────■─────────────
    #      ┌──────────┐ │ ┌─────────┐
    # q_1: ┤ Ry(-π/4) ├─■─┤ Ry(π/4) ├
    #      └──────────┘   └─────────┘
    ch_to_cz_ry = QuantumCircuit(2)
    ch_to_cz_ry.ry(-pi / 4, 1)
    ch_to_cz_ry.cz(0, 1)
    ch_to_cz_ry.ry(pi / 4, 1)
    lib.add_equivalence(std.CHGate(), ch_to_cz_ry)

    # This is the same as the `ry`+`cz` rule above, but explicitly decomposed into a discrete basis.
    # q_0: ──────────────────────■────────────────────
    #      ┌────┐┌─────┐┌──────┐ │ ┌────┐┌───┐┌──────┐
    # q_1: ┤ √X ├┤ Tdg ├┤ √Xdg ├─■─┤ √X ├┤ T ├┤ √Xdg ├
    #      └────┘└─────┘└──────┘   └────┘└───┘└──────┘
    ch_to_cz_discrete = QuantumCircuit(2)
    ch_to_cz_discrete.sx(1)
    ch_to_cz_discrete.tdg(1)
    ch_to_cz_discrete.sxdg(1)
    ch_to_cz_discrete.cz(0, 1)
    ch_to_cz_discrete.sx(1)
    ch_to_cz_discrete.t(1)
    ch_to_cz_discrete.sxdg(1)
    lib.add_equivalence(std.CHGate(), ch_to_cz_discrete)

    # ECR into others.

    #      ┌───┐      ┌───┐
    # q_0: ┤ S ├───■──┤ X ├
    #      ├───┴┐┌─┴─┐└───┘
    # q_1: ┤ √X ├┤ X ├─────
    #      └────┘└───┘
    ecr_to_cx = QuantumCircuit(2, global_phase=-pi / 4)
    ecr_to_cx.s(0)
    ecr_to_cx.sx(1)
    ecr_to_cx.cx(0, 1)
    ecr_to_cx.x(0)
    lib.add_equivalence(std.ECRGate(), ecr_to_cx)

    #                   ┌───┐ ┌─────┐
    # q_0: ───────────■─┤ X ├─┤ Sdg ├
    #      ┌───┐┌───┐ │ ├───┴┐├───┬─┘
    # q_1: ┤ H ├┤ Z ├─■─┤ √X ├┤ S ├─
    #      └───┘└───┘   └────┘└───┘
    ecr_to_cz = QuantumCircuit(2)
    ecr_to_cz.h(1)
    ecr_to_cz.z(1)
    ecr_to_cz.cz(0, 1)
    ecr_to_cz.x(0)
    ecr_to_cz.sdg(0)
    ecr_to_cz.sx(1)
    ecr_to_cz.s(1)
    lib.add_equivalence(std.ECRGate(), ecr_to_cz)


def local_equivalence_csx(lib: EquivalenceLibrary):
    """Relations between gates that are locally equivalent to CSX.

    This handles `cs`, `csdg` and `csx`."""
    # It's sufficient to supply "everything to `cs`" and "`cs` to everything".

    # q_0: ───────■────────
    #      ┌───┐┌─┴──┐┌───┐
    # q_1: ┤ H ├┤ Sx ├┤ H ├
    #      └───┘└────┘└───┘
    cs_to_csx = QuantumCircuit(2)
    cs_to_csx.h(1)
    cs_to_csx.csx(0, 1)
    cs_to_csx.h(1)
    lib.add_equivalence(std.CSGate(), cs_to_csx)

    #                  ┌───┐
    # q_0: ────────■───┤ S ├
    #      ┌───┐┌──┴──┐├───┤
    # q_1: ┤ X ├┤ Sdg ├┤ X ├
    #      └───┘└─────┘└───┘
    cs_to_csdg = QuantumCircuit(2)
    cs_to_csdg.x(1)
    cs_to_csdg.csdg(0, 1)
    cs_to_csdg.s(0)
    cs_to_csdg.x(1)
    lib.add_equivalence(std.CSGate(), cs_to_csdg)

    # q_0: ───────■───────
    #      ┌───┐┌─┴─┐┌───┐
    # q_1: ┤ H ├┤ S ├┤ H ├
    #      └───┘└───┘└───┘
    csx_to_cs = QuantumCircuit(2)
    csx_to_cs.h(1)
    csx_to_cs.cs(0, 1)
    csx_to_cs.h(1)
    lib.add_equivalence(std.CSXGate(), csx_to_cs)

    #                ┌─────┐
    # q_0: ───────■──┤ Sdg ├
    #      ┌───┐┌─┴─┐└┬───┬┘
    # q_1: ┤ X ├┤ S ├─┤ X ├─
    #      └───┘└───┘ └───┘
    csdg_to_cs = QuantumCircuit(2)
    csdg_to_cs.x(1)
    csdg_to_cs.cs(0, 1)
    csdg_to_cs.sdg(0)
    csdg_to_cs.x(1)
    lib.add_equivalence(std.CSdgGate(), csdg_to_cs)


def local_equivalence_iswap(lib: EquivalenceLibrary):
    """Relations between the gates that are locally equivalent to iSWAP.

    This handles `iswap` and `dcx`."""
    #       ┌───┐ ┌─────┐┌────────┐
    # q_0: ─┤ H ├─┤ Sdg ├┤0       ├─────
    #      ┌┴───┴┐└─────┘│  Iswap │┌───┐
    # q_1: ┤ Sdg ├───────┤1       ├┤ H ├
    #      └─────┘       └────────┘└───┘
    dcx_to_iswap = QuantumCircuit(2)
    dcx_to_iswap.h(0)
    dcx_to_iswap.sdg(0)
    dcx_to_iswap.sdg(1)
    dcx_to_iswap.iswap(0, 1)
    dcx_to_iswap.h(1)
    lib.add_equivalence(std.DCXGate(), dcx_to_iswap)

    #      ┌───┐┌───┐┌──────┐
    # q_0: ┤ S ├┤ H ├┤0     ├─────
    #      ├───┤└───┘│  Dcx │┌───┐
    # q_1: ┤ S ├─────┤1     ├┤ H ├
    #      └───┘     └──────┘└───┘
    iswap_to_dcx = QuantumCircuit(2)
    iswap_to_dcx.s(0)
    iswap_to_dcx.s(1)
    iswap_to_dcx.h(0)
    iswap_to_dcx.dcx(0, 1)
    iswap_to_dcx.h(1)
    lib.add_equivalence(std.iSwapGate(), iswap_to_dcx)


def local_equivalence_ising(lib: EquivalenceLibrary):
    """Relations between the gates that are locally equivalent to an Ising rotation (rotation around
    some 2q Pauli axis).

    This handles `crx`, `cry`, `crz`, `rxx`, `ryy`, `rzz`, `rzx`, `cp`, `cu1`, `cu3` and `cu`."""

    # We relate each "controlled" gate to a logical choice for its two-qubit Pauli rotation (and
    # vice-versa), then make sure all two-Pauli terms can convert between themselves.  Where
    # reasonable, we put the single-qubit gates _after_ the 2q gate; this is an arbitary choice, but
    # the consistency is meant to allow the 1q optimiser to collapse chains of rotations that stem
    # from chains of translations.

    alpha, beta, gamma, delta = (Parameter("α"), Parameter("β"), Parameter("γ"), Parameter("δ"))

    # `cp`, `cu1`, `crz` into `rzz` (and between each other).

    # q_0: ─■────────────────────
    #       │ZZ(-α/2) ┌─────────┐
    # q_1: ─■─────────┤ Rz(α/2) ├
    #                 └─────────┘
    crz_to_rzz = QuantumCircuit(2)
    crz_to_rzz.rzz(-alpha / 2, 0, 1)
    crz_to_rzz.rz(alpha / 2, 1)
    lib.add_equivalence(std.CRZGate(alpha), crz_to_rzz)

    # q_0: ─────■───────────────
    #      ┌────┴─────┐┌───────┐
    # q_1: ┤ Rz(-2*α) ├┤ Rz(α) ├
    #      └──────────┘└───────┘
    rzz_to_crz = QuantumCircuit(2)
    rzz_to_crz.crz(-2 * alpha, 0, 1)
    rzz_to_crz.rz(alpha, 1)
    lib.add_equivalence(std.RZZGate(alpha), rzz_to_crz)

    #                 ┌─────────┐
    # q_0: ─■─────────┤ Rz(α/2) ├
    #       │ZZ(-α/2) ├─────────┤
    # q_1: ─■─────────┤ Rz(α/2) ├
    #                 └─────────┘
    cp_to_rzz = QuantumCircuit(2, global_phase=alpha / 4)
    cp_to_rzz.rzz(-alpha / 2, 0, 1)
    cp_to_rzz.rz(alpha / 2, 0)
    cp_to_rzz.rz(alpha / 2, 1)
    lib.add_equivalence(std.CPhaseGate(alpha), cp_to_rzz)

    #                ┌───────┐
    # q_0: ─■────────┤ Rz(α) ├
    #       │P(-2*α) ├───────┤
    # q_1: ─■────────┤ Rz(α) ├
    #                └───────┘
    rzz_to_cp = QuantumCircuit(2, global_phase=alpha / 2)
    rzz_to_cp.cp(-2 * alpha, 0, 1)
    rzz_to_cp.rz(alpha, 0)
    rzz_to_cp.rz(alpha, 1)
    lib.add_equivalence(std.RZZGate(alpha), rzz_to_cp)

    # q_0: ─■────
    #       │U1(α)
    # q_1: ─■────
    cp_to_cu1 = QuantumCircuit(2)
    cp_to_cu1.append(std.CU1Gate(alpha), [0, 1])
    lib.add_equivalence(std.CPhaseGate(alpha), cp_to_cu1)

    # q_0: ─■────
    #       │P(α)
    # q_1: ─■────
    cu1_to_cp = QuantumCircuit(2)
    cu1_to_cp.cp(alpha, 0, 1)
    lib.add_equivalence(std.CU1Gate(alpha), cu1_to_cp)

    # `crx` and `cry` into sensible two-Pauli rotations (and back again).

    #      ┌────────────┐
    # q_0: ┤0           ├───────────
    #      │  Rzx(-α/2) │┌─────────┐
    # q_1: ┤1           ├┤ Rx(α/2) ├
    #      └────────────┘└─────────┘
    crx_to_rzx = QuantumCircuit(2)
    crx_to_rzx.rzx(-alpha / 2, 0, 1)
    crx_to_rzx.rx(alpha / 2, 1)
    lib.add_equivalence(std.CRXGate(alpha), crx_to_rzx)

    # q_0: ─────■───────────────
    #      ┌────┴─────┐┌───────┐
    # q_1: ┤ Rx(-2*α) ├┤ Rx(α) ├
    #      └──────────┘└───────┘
    rzx_to_crx = QuantumCircuit(2)
    rzx_to_crx.crx(-2 * alpha, 0, 1)
    rzx_to_crx.rx(alpha, 1)
    lib.add_equivalence(std.RZXGate(alpha), rzx_to_crx)

    #             ┌────────────┐
    # q_0: ───────┤0           ├────────────────
    #      ┌─────┐│  Rzx(-α/2) │┌─────────┐┌───┐
    # q_1: ┤ Sdg ├┤1           ├┤ Rx(α/2) ├┤ S ├
    #      └─────┘└────────────┘└─────────┘└───┘
    cry_to_rzx = QuantumCircuit(2)
    cry_to_rzx.sdg(1)
    cry_to_rzx.rzx(-alpha / 2, 0, 1)
    cry_to_rzx.rx(alpha / 2, 1)
    cry_to_rzx.s(1)
    lib.add_equivalence(std.CRYGate(alpha), cry_to_rzx)

    # q_0: ──────────■──────────────────────
    #      ┌───┐┌────┴─────┐┌───────┐┌─────┐
    # q_1: ┤ S ├┤ Ry(-2*α) ├┤ Ry(α) ├┤ Sdg ├
    #      └───┘└──────────┘└───────┘└─────┘
    rzx_to_cry = QuantumCircuit(2)
    rzx_to_cry.s(1)
    rzx_to_cry.cry(-2 * alpha, 0, 1)
    rzx_to_cry.ry(alpha, 1)
    rzx_to_cry.sdg(1)
    lib.add_equivalence(std.RZXGate(alpha), rzx_to_cry)

    # ... ok, now we've handled convering all `crA` gates (for all `A`) into an equivalent `rBC`
    # (for one particular `BC`) and back again, so all that's left for the single-parameter
    # rotations is conversions between the two-Pauli rotations.  We just need the connection graph
    # to be fully connected between 2q gates, not necessarily direct all-to-all edges.

    #      ┌───┐┌─────────┐┌───┐
    # q_0: ┤ H ├┤0        ├┤ H ├
    #      └───┘│  Rxx(ϴ) │└───┘
    # q_1: ─────┤1        ├─────
    #           └─────────┘
    rzx_to_rxx = QuantumCircuit(2)
    rzx_to_rxx.h(0)
    rzx_to_rxx.rxx(alpha, 0, 1)
    rzx_to_rxx.h(0)
    lib.add_equivalence(std.RZXGate(alpha), rzx_to_rxx)

    #      ┌───┐┌─────────┐┌───┐
    # q_0: ┤ H ├┤0        ├┤ H ├
    #      └───┘│  Rzx(ϴ) │└───┘
    # q_1: ─────┤1        ├─────
    #           └─────────┘
    rxx_to_rzx = QuantumCircuit(2)
    rxx_to_rzx.h(0)
    rxx_to_rzx.rzx(alpha, 0, 1)
    rxx_to_rzx.h(0)
    lib.add_equivalence(std.RXXGate(alpha), rxx_to_rzx)

    #           ┌─────────┐
    # q_0: ─────┤0        ├─────
    #      ┌───┐│  Rzz(ϴ) │┌───┐
    # q_1: ┤ H ├┤1        ├┤ H ├
    #      └───┘└─────────┘└───┘
    rzx_to_rzz = QuantumCircuit(2)
    rzx_to_rzz.h(1)
    rzx_to_rzz.rzz(alpha, 0, 1)
    rzx_to_rzz.h(1)
    lib.add_equivalence(std.RZXGate(alpha), rzx_to_rzz)

    #           ┌─────────┐
    # q_0: ─────┤0        ├─────
    #      ┌───┐│  Rzx(ϴ) │┌───┐
    # q_1: ┤ H ├┤1        ├┤ H ├
    #      └───┘└─────────┘└───┘
    rzz_to_rzx = QuantumCircuit(2)
    rzz_to_rzx.h(1)
    rzz_to_rzx.rzx(alpha, 0, 1)
    rzz_to_rzx.h(1)
    lib.add_equivalence(std.RZZGate(alpha), rzz_to_rzx)

    #      ┌───┐┌─────────┐┌─────┐
    # q_0: ┤ S ├┤0        ├┤ Sdg ├
    #      ├───┤│  Ryy(α) │├─────┤
    # q_1: ┤ S ├┤1        ├┤ Sdg ├
    #      └───┘└─────────┘└─────┘
    rxx_to_ryy = QuantumCircuit(2)
    rxx_to_ryy.s(0)
    rxx_to_ryy.s(1)
    rxx_to_ryy.ryy(alpha, 0, 1)
    rxx_to_ryy.sdg(0)
    rxx_to_ryy.sdg(1)
    lib.add_equivalence(std.RXXGate(alpha), rxx_to_ryy)

    #      ┌─────┐┌─────────┐┌───┐
    # q_0: ┤ Sdg ├┤0        ├┤ S ├
    #      ├─────┤│  Rxx(α) │├───┤
    # q_1: ┤ Sdg ├┤1        ├┤ S ├
    #      └─────┘└─────────┘└───┘
    ryy_to_rxx = QuantumCircuit(2)
    ryy_to_rxx.sdg(0)
    ryy_to_rxx.sdg(1)
    ryy_to_rxx.rxx(alpha, 0, 1)
    ryy_to_rxx.s(0)
    ryy_to_rxx.s(1)
    lib.add_equivalence(std.RYYGate(alpha), ryy_to_rxx)

    # Lastly, `cu` and `cu3`.  These gates actually are locally equivalent to the other Ising gates,
    # but determining the exact equivalence symbolically is a bit trig heavy, which is unfortunate.
    # In principle, you can convert the Euler angles of `u3` to an axis-angle representation, then
    # sandwich an `rxx` between rotations of the true axis to the XX one and back again.  We don't
    # currently have that rule in the default set mostly because it's not clear it would be useful,
    # and it's an ugly amount of additional `ParameterExpression` code we might end up needing to
    # support for longer.

    #                   ┌──────┐
    # q_0: ──────■──────┤ P(δ) ├
    #      ┌─────┴─────┐└──────┘
    # q_1: ┤ U3(α,β,γ) ├────────
    #      └───────────┘
    cu_to_cu3 = QuantumCircuit(2)
    cu_to_cu3.append(std.CU3Gate(alpha, beta, gamma), [0, 1])
    cu_to_cu3.p(delta, 0)
    lib.add_equivalence(std.CUGate(alpha, beta, gamma, delta), cu_to_cu3)

    # q_0: ──────■───────
    #      ┌─────┴──────┐
    # q_1: ┤ U(α,β,γ,0) ├
    #      └────────────┘
    cu3_to_cu = QuantumCircuit(2)
    cu3_to_cu.cu(alpha, beta, gamma, 0, 0, 1)
    lib.add_equivalence(std.CU3Gate(alpha, beta, gamma), cu3_to_cu)


def between_locally_equivalent_groups_2q(lib: EquivalenceLibrary):
    """Relations linking the three different groups of locally equivalence gates.

    SWAP is handled separately because of its importance in routing; it has a large number of
    rules."""
    alpha = Parameter("α")

    # We don't need everything to everything, provided the basis-translator rules can find paths all
    # the way through every group.  We do want to ensure that each group pair has a rule that uses
    # the minimum amount of 2q gates for the conversion, though.

    # q_0: ──■─────■───
    #      ┌─┴──┐┌─┴──┐
    # q_1: ┤ Sx ├┤ Sx ├
    #      └────┘└────┘
    cx_to_csx = QuantumCircuit(2)
    cx_to_csx.csx(0, 1)
    cx_to_csx.csx(0, 1)
    lib.add_equivalence(std.CXGate(), cx_to_csx)

    #      ┌───┐     ┌────────┐┌───┐     ┌────────┐┌───┐┌───┐
    # q_0: ┤ H ├─────┤0       ├┤ X ├─────┤0       ├┤ H ├┤ S ├─────
    #      ├───┤┌───┐│  Iswap │├───┤┌───┐│  Iswap │├───┤├───┤┌───┐
    # q_1: ┤ X ├┤ H ├┤1       ├┤ X ├┤ H ├┤1       ├┤ S ├┤ X ├┤ H ├
    #      └───┘└───┘└────────┘└───┘└───┘└────────┘└───┘└───┘└───┘
    cx_to_iswap = QuantumCircuit(2, global_phase=3 * pi / 4)
    cx_to_iswap.h(0)
    cx_to_iswap.x(1)
    cx_to_iswap.h(1)
    cx_to_iswap.iswap(0, 1)
    cx_to_iswap.x(0)
    cx_to_iswap.x(1)
    cx_to_iswap.h(1)
    cx_to_iswap.iswap(0, 1)
    cx_to_iswap.h(0)
    cx_to_iswap.s(0)
    cx_to_iswap.s(1)
    cx_to_iswap.x(1)
    cx_to_iswap.h(1)
    lib.add_equivalence(std.CXGate(), cx_to_iswap)

    #               ┌───┐
    # q_0: ────■────┤ S ├
    #      ┌───┴───┐└───┘
    # q_1: ┤ Rx(π) ├─────
    #      └───────┘
    cx_to_crx = QuantumCircuit(2)
    cx_to_crx.crx(pi, 0, 1)
    cx_to_crx.s(0)
    lib.add_equivalence(std.CXGate(), cx_to_crx)

    #      ┌───────────┐┌─────┐
    # q_0: ┤0          ├┤ Sdg ├─
    #      │  Rzx(π/2) │├─────┴┐
    # q_1: ┤1          ├┤ √Xdg ├
    #      └───────────┘└──────┘
    cx_to_rzx = QuantumCircuit(2, global_phase=pi / 4)
    cx_to_rzx.rzx(pi / 2, 0, 1)
    cx_to_rzx.sdg(0)
    cx_to_rzx.sxdg(1)
    lib.add_equivalence(std.CXGate(), cx_to_rzx)

    # There's probably a cleaner way of doing this decomposition, but this was the best I came up
    # with in the amount of time I felt comfortable spending on it. ---Jake
    #      ┌───┐┌─────┐     ┌───┐      ┌───┐ ┌─────┐
    # q_0: ┤ Z ├┤ Tdg ├──■──┤ X ├──■───┤ X ├─┤ Sdg ├
    #      ├───┤├─────┤┌─┴─┐├───┤┌─┴─┐┌┴───┴┐└─────┘
    # q_1: ┤ X ├┤ Sdg ├┤ X ├┤ T ├┤ X ├┤ Tdg ├───────
    #      └───┘└─────┘└───┘└───┘└───┘└─────┘
    cs_to_cx = QuantumCircuit(2, global_phase=pi / 4)
    cs_to_cx.z(0)
    cs_to_cx.tdg(0)
    cs_to_cx.x(1)
    cs_to_cx.sdg(1)
    cs_to_cx.cx(0, 1)
    cs_to_cx.x(0)
    cs_to_cx.t(1)
    cs_to_cx.cx(0, 1)
    cs_to_cx.x(0)
    cs_to_cx.sdg(0)
    cs_to_cx.tdg(1)
    lib.add_equivalence(std.CSGate(), cs_to_cx)

    # q_0: ─■───────
    #       │P(π/2)
    # q_1: ─■───────
    cs_to_cp = QuantumCircuit(2)
    cs_to_cp.cp(pi / 2, 0, 1)
    lib.add_equivalence(std.CSGate(), cs_to_cp)

    #           ┌───┐
    # q_0: ──■──┤ X ├
    #      ┌─┴─┐└─┬─┘
    # q_1: ┤ X ├──■──
    #      └───┘
    dcx_to_cx = QuantumCircuit(2)
    dcx_to_cx.cx(0, 1)
    dcx_to_cx.cx(1, 0)
    lib.add_equivalence(std.DCXGate(), dcx_to_cx)

    # q_0: ──■─────────────■──
    #      ┌─┴─┐┌───────┐┌─┴─┐
    # q_1: ┤ X ├┤ Rz(α) ├┤ X ├
    #      └───┘└───────┘└───┘
    rzz_to_cx = QuantumCircuit(2)
    rzz_to_cx.cx(0, 1)
    rzz_to_cx.rz(alpha, 1)
    rzz_to_cx.cx(0, 1)
    lib.add_equivalence(std.RZZGate(alpha), rzz_to_cx)

    #
    # q_0: ─■───────────■─
    #       │ ┌───────┐ │
    # q_1: ─■─┤ Rx(α) ├─■─
    #         └───────┘
    rzx_to_cz = QuantumCircuit(2)
    rzx_to_cz.cz(0, 1)
    rzx_to_cz.rx(alpha, 1)
    rzx_to_cz.cz(0, 1)
    lib.add_equivalence(std.RZXGate(alpha), rzx_to_cz)

    # This is the `rzx_to_cz` rule, but replacing the `cz`s with `swap - iswap - sdg` sequences,
    # then eliding the swaps.  Probably there's a neater way, but this is at least 2q-optimal.
    #      ┌────────┐┌─────┐┌───────┐┌─────┐┌────────┐
    # q_0: ┤0       ├┤ Sdg ├┤ Rx(α) ├┤ Sdg ├┤0       ├
    #      │  Iswap │└┬───┬┘└───────┘└─────┘│  Iswap │
    # q_1: ┤1       ├─┤ Z ├─────────────────┤1       ├
    #      └────────┘ └───┘                 └────────┘
    rzx_to_iswap = QuantumCircuit(2)
    rzx_to_iswap.iswap(0, 1)
    rzx_to_iswap.sdg(0)
    rzx_to_iswap.rx(alpha, 0)
    rzx_to_iswap.sdg(0)
    rzx_to_iswap.z(1)
    rzx_to_iswap.iswap(0, 1)
    lib.add_equivalence(std.RZXGate(alpha), rzx_to_iswap)


def reverse_gate_direction_2q(lib: EquivalenceLibrary):
    """Rules for reversing a 2q gate in terms of itself and 1q gates."""
    # We don't need every possible reversal here.  As long as we've got one from each 2q equivalence
    # group, the other "within local equivalence group" rules will allow direction switching of any
    # operation.  It might not always be the absolute _best_ in all cases, but 1q optimisation would
    # typically be able to tidy most stuff up again.  Since we have a symmetric gate in each
    # equivalence class, we basically get all the rules for free.
    alpha, beta = Parameter("α"), Parameter("β")
    for symmetric in (
        std.CZGate(),
        std.SwapGate(),
        std.RXXGate(alpha),
        std.RYYGate(alpha),
        std.RZZGate(alpha),
        std.XXMinusYYGate(alpha, beta),
        std.iSwapGate(),
        std.CSGate(),
        std.CSdgGate(),
    ):
        qc = QuantumCircuit(2)
        qc.append(symmetric, [1, 0], [])
        lib.add_equivalence(symmetric, qc)


def efficient_swap(lib: EquivalenceLibrary):
    """Efficient rules translating SWAP to various targets."""
    #           ┌───┐
    # q_0: ──■──┤ X ├──■──
    #      ┌─┴─┐└─┬─┘┌─┴─┐
    # q_1: ┤ X ├──■──┤ X ├
    #      └───┘     └───┘
    swap_to_cx = QuantumCircuit(2)
    swap_to_cx.cx(0, 1)
    swap_to_cx.cx(1, 0)
    swap_to_cx.cx(0, 1)
    lib.add_equivalence(std.SwapGate(), swap_to_cx)

    #      ┌────┐   ┌────┐   ┌────┐
    # q_0: ┤ √X ├─■─┤ √X ├─■─┤ √X ├─■─
    #      ├────┤ │ ├────┤ │ ├────┤ │
    # q_1: ┤ √X ├─■─┤ √X ├─■─┤ √X ├─■─
    #      └────┘   └────┘   └────┘
    swap_to_cz = QuantumCircuit(2, global_phase=-pi / 2)
    swap_to_cz.sx(0)
    swap_to_cz.sx(1)
    swap_to_cz.cz(0, 1)
    swap_to_cz.sx(0)
    swap_to_cz.sx(1)
    swap_to_cz.cz(0, 1)
    swap_to_cz.sx(0)
    swap_to_cz.sx(1)
    swap_to_cz.cz(0, 1)
    lib.add_equivalence(std.SwapGate(), swap_to_cz)

    #      ┌───┐ ┌──────┐┌────┐┌──────┐┌───┐ ┌──────┐
    # q_0: ┤ S ├─┤0     ├┤ √X ├┤1     ├┤ S ├─┤0     ├
    #      ├───┴┐│  Ecr │├───┬┘│  Ecr │├───┴┐│  Ecr │
    # q_1: ┤ √X ├┤1     ├┤ S ├─┤0     ├┤ √X ├┤1     ├
    #      └────┘└──────┘└───┘ └──────┘└────┘└──────┘
    swap_to_ecr = QuantumCircuit(2, global_phase=3 * pi / 4)
    swap_to_ecr.s(0)
    swap_to_ecr.sx(1)
    swap_to_ecr.ecr(0, 1)
    swap_to_ecr.s(1)
    swap_to_ecr.sx(0)
    swap_to_ecr.ecr(1, 0)
    swap_to_ecr.s(0)
    swap_to_ecr.sx(1)
    swap_to_ecr.ecr(0, 1)
    lib.add_equivalence(std.SwapGate(), swap_to_ecr)

    #      ┌────────┐      ┌────────┐┌────┐┌────────┐
    # q_0: ┤0       ├──────┤0       ├┤ √X ├┤0       ├──────
    #      │  Iswap │┌────┐│  Iswap │└────┘│  Iswap │┌────┐
    # q_1: ┤1       ├┤ √X ├┤1       ├──────┤1       ├┤ √X ├
    #      └────────┘└────┘└────────┘      └────────┘└────┘
    swap_to_iswap = QuantumCircuit(2, global_phase=-pi / 2)
    swap_to_iswap.iswap(0, 1)
    swap_to_iswap.sx(1)
    swap_to_iswap.iswap(0, 1)
    swap_to_iswap.sx(0)
    swap_to_iswap.iswap(0, 1)
    swap_to_iswap.sx(1)
    lib.add_equivalence(std.SwapGate(), swap_to_iswap)

    #      ┌──────┐
    # q_0: ┤0     ├──■──
    #      │  Dcx │┌─┴─┐
    # q_1: ┤1     ├┤ X ├
    #      └──────┘└───┘
    swap_to_dcx_cx = QuantumCircuit(2)
    swap_to_dcx_cx.dcx(0, 1)
    swap_to_dcx_cx.cx(0, 1)
    lib.add_equivalence(std.SwapGate(), swap_to_dcx_cx)

    # This is, up to 1q gates, the same as the (DCX, CX) decomposition above, but the (iSwap, CZ)
    # pairing is more natural for superconducting hardware, and this form is already 1q optimal.
    #      ┌─────┐┌────────┐
    # q_0: ┤ Sdg ├┤0       ├─■─
    #      ├─────┤│  Iswap │ │
    # q_1: ┤ Sdg ├┤1       ├─■─
    #      └─────┘└────────┘
    swap_to_iswap_cz = QuantumCircuit(2)
    swap_to_iswap_cz.sdg(0)
    swap_to_iswap_cz.sdg(1)
    swap_to_iswap_cz.iswap(0, 1)
    swap_to_iswap_cz.cz(0, 1)
    lib.add_equivalence(std.SwapGate(), swap_to_iswap_cz)


def multi_q(lib: EquivalenceLibrary):
    """Relations for 3+q operators.

    There's not much structure to this group; we don't do a huge amount of synthesis in terms of
    high-degree gates, so we don't have a huge corpus of knowledge to draw on for what we need."""
    # q_0: ────────■─────────────────■────■───
    #            ┌─┴─┐┌─────┐      ┌─┴─┐  │
    # q_1: ──■───┤ X ├┤ Sdg ├──■───┤ X ├──┼───
    #      ┌─┴──┐├───┤└─────┘┌─┴──┐├───┤┌─┴──┐
    # q_2: ┤ Sx ├┤ Z ├───────┤ Sx ├┤ Z ├┤ Sx ├
    #      └────┘└───┘       └────┘└───┘└────┘
    ccx_to_cx_csx = QuantumCircuit(3)
    ccx_to_cx_csx.csx(1, 2)
    ccx_to_cx_csx.cx(0, 1)
    ccx_to_cx_csx.sdg(1)
    ccx_to_cx_csx.z(2)
    ccx_to_cx_csx.csx(1, 2)
    ccx_to_cx_csx.cx(0, 1)
    ccx_to_cx_csx.z(2)
    ccx_to_cx_csx.csx(0, 2)
    lib.add_equivalence(std.CCXGate(), ccx_to_cx_csx)
