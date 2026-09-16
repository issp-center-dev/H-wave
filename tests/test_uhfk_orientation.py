"""Spec 2026-09-16 section 2.4, gate G-UHFk: UHFk END TO END on the
reversed declaration.

An off-site two-body row ``(r, a, b, v)`` is read in the DOCUMENTED
orientation since issue #193 -- ``v n_{j,a} n_{j+r,b}``, orbital ``a`` in
the original cell and ``b`` in the cell displaced by ``r``. The previous
release (H-wave 2.0.0, the reference revision pinned by
``tests/test_flex_second_order_compat.DEVELOP_COMMIT``) read the two cells
the other way round. The compatibility promise that goes with the change
is therefore an EXACT one:

    this branch's UHFk on a declaration ``F``
        == the reference revision's UHFk on ``F^rev``

where ``F^rev`` negates the displacement of every off-site two-body row
and leaves the orbital indices and the on-site rows alone
(``tests/test_flex_hf_scf._reversed_interaction_dir`` writes it). The
in-process, table-level version of the same statement is
``tests/test_flex_second_order_compat`` (``new(F) == legacy(F^rev)``, at
``np.array_equal``); this module is the end-to-end one -- two real
``hwave.qlms.run`` processes, the two source trees, the shipped
``energy.dat`` and ``green.dat.npz``.

CI provisions the reference revision (a detached worktree at
``DEVELOP_COMMIT`` under ``RUNNER_TEMP``; see
``.github/workflows/ci-python39.yml``) and sets
``HWAVE_REQUIRE_DEVELOP_COMPARISON=1``, so this gate RUNS there and a
missing, wrong-revision or dirty reference is a failure rather than a
skip. Locally the flag is unset and the gate skips with a reason when no
such checkout is at hand; ``TestFixturesAreOrientationSensitive`` below
needs no reference tree and runs either way.

WHY THE COMPARISON IS ``np.array_equal`` AND NOT A TOLERANCE
------------------------------------------------------------
Because the two runs do the SAME arithmetic. The orientation step is a
conjugate transpose at fixed ``r`` of the assembled displacement table,
and the reversal of the declaration produces that same table out of the
reference revision's builder -- bit for bit, not to round-off (Task 2's
kernel test pins exactly that with ``np.array_equal``). The tables being
identical, the whole SCF trajectory is identical: same initial Green
(zero, or the COMMITTED seed :data:`_SEED` -- a file both source trees
read, not a draw from a random generator they would each have to agree
about), same mixing, same eigen-decompositions, same iteration count. A
tolerance here would hide a real difference in the mean field. The only
thing that could separate the two processes is a different numpy/BLAS on
the two ``PYTHONPATH``s, which is not a thing on one machine. The
comparison below therefore carries NO tolerance at all,
not even a round-off fallback: there is no measurement such a fallback
could be calibrated from, and it would be exactly the place where a real
difference in the mean field could hide.

WHICH CASES THE ORIENTATION ACTUALLY MOVES (and why the rest cannot)
--------------------------------------------------------------------
The identity above is asserted for every case. Its anti-vacuity control
-- the reference revision on the UNREVERSED ``F``, which must land
somewhere else -- is only meaningful where the orientation is VISIBLE at
all, and two documented classes of case it is not. Both are structural,
derived rather than observed, and asserted here in the coinciding
direction so that a change making them visible cannot pass unnoticed:

* **Hartree-only** (``flag_fock = false``) for ``CoulombInter`` / ``Hund``
  / ``Ising``: the Hartree term contracts the table only through
  ``sum_r J_ab(r)``. All three places that read it do:
  ``hartree_fock.accumulate_hf`` (the normal-mode mean field: ``hh1 =
  einsum('rab, stb -> rsta', jab_r, hh0)``, then a sum over ``r``), the
  spin-orbital branch of ``uhfk.UHFk._make_ham`` (the same ``hh0``/
  ``hh1``/``hh2``/``sum over r`` before ``_virtual_ham_to_so``), and
  ``uhfk.UHFk._calc_energy``, which is what ``energy.dat`` prints
  (Hartree-only: ``ee = einsum('rab, rab->', jab_r, w1b)`` with ``w1b``
  broadcast over ``r``, i.e. again only ``sum_r J_ab(r)``). The
  reversal/Hermitian closure ``_reverse_closed`` makes that sum
  HERMITIAN, and the orientation maps it to its own conjugate transpose,
  i.e. to itself. Nothing downstream can see the change, in either mode.
  ``PairHop`` is NOT in this class: its term (the ``hh6``/``hh7`` block)
  is not gated by ``flag_fock`` and reads the displacement resolved, so
  the off-site complex PairHop moves with the Fock term on or off.

* **A single orbital**: at ``norb = 1`` the transpose is the identity, so
  the oriented table is the COMPLEX CONJUGATE of the original one. With a
  real transfer, ``F^rev``'s whole Hamiltonian is then the complex
  conjugate of ``F``'s, and every real observable ``energy.dat`` carries
  is invariant under complex conjugation -- a single-orbital chain cannot
  make the orientation visible in the energy no matter what phase the
  off-site amplitude carries (the complex PairHop amplitude here does
  move the TABLE, which is what Task 2's kernel gate sees). The
  fixture is kept because the end-to-end identity is worth pinning on the
  single-orbital path too, and because the statement "this one cannot
  move" is itself falsifiable.

The fixtures under ``tests/uhfk_orientation/`` are a ``CellShape =
[4, 1, 1]`` chain with COMPLEX inter-orbital hopping and ``t_01(+x) !=
conj(t_10(+x))`` (the amplitudes of
``tests/test_flex_second_order_ed_chain._fx_complex``), so the band is
neither inversion- nor transpose-symmetric and the two readings of a bond
really are different Hamiltonians. ``coulombinter.dat`` carries an
ASYMMETRIC inter-orbital bond (``v_12(+x) = 0.4 != v_21(+x) = 0.25``),
``hund.dat`` an inter-orbital ``J``, ``pairhop.dat`` a COMPLEX off-site
amplitude; every case also runs an on-site ``coulombintra.dat`` so the
mean field is nontrivial. ``ising.dat``, ``exchange.dat`` and
``pairlift.dat`` carry the same inter-orbital bond for the three types the
bond and PairHop cases never reach -- their mean-field expressions are
separate lines of ``hartree_fock.accumulate_hf``, duplicated again in the
spin-orbital branch of ``uhfk.UHFk._make_ham``. ``pairlift.dat``'s
amplitude is deliberately LARGE (3.5, against a bandwidth of order 1):
PairLift contracts only the SPIN-OFF-DIAGONAL part of the density, UHFk's
spin-collinear solution is a fixed point of the loop, and below the
symmetry-breaking threshold the term contributes exactly zero however it
is oriented. Those two cases therefore start from a COMMITTED
spin-mixing initial Green function (:data:`_SEED`, built by
:func:`_spin_mixed_seed`) rather than from UHFk's zero density -- and the
amplitude still has to clear the threshold: MEASURED from that seed,
PairLift contributes 2.3e-22 at 1.0 (the transverse component decays and
the case is vacuous for it), 2.15e-01 at 2.0 and 3.91e-01 at 3.5.
``ising.dat`` is raised for the same reason and a weaker one (2.0): in the
broken-symmetry
state PairLift produces, the Ising contribution is suppressed, and at 0.3 it
cleared the anti-vacuity floor by only 2.6x (2.57e-06, i.e. 4e-07 of the
total energy); at 2.0 it is 5.97e-04, and the reference control's gap on the
unreversed declaration widens from 1.93e-03 to 3.93e-03. The declaration is
a test fixture, not a physical model; what it has to do is make the
contraction run and leave a margin worth measuring.

The ``*_so`` directories hold the same band written in spin-orbital
indices (``so = 2a + s``, spin-diagonal) for the
``enable_spin_orbital = true`` twin of each case: UHFk supports every
two-body type there (``uhfk._check_spin_orbital_compatibility``), the
interaction files keep PHYSICAL orbital indices, and ``Ncond`` counts
total electrons.
"""

import collections
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from tests.heavy_tests import heavy

_FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "uhfk_orientation")

#: What the subprocess runs: one real UHFk calculation, driven exactly as a
#: user drives it (``hwave.qlms.run`` on a parsed input), in the current
#: working directory -- which the runner has made a COPY of the fixture, so
#: no ``output/`` is ever written next to the committed inputs.
_RUN_SCRIPT = ("import json, sys; import hwave.qlms; "
               "hwave.qlms.run(input_dict=json.loads(sys.argv[1]))")

Case = collections.namedtuple(
    "Case", "label fixture norb files twobody spin_orbital fock moves why")


_BOND = {"CoulombIntra": "coulombintra.dat",
         "CoulombInter": "coulombinter.dat", "Hund": "hund.dat"}
_BOND_FILES = ("coulombinter.dat", "hund.dat")
_PAIRHOP = {"CoulombIntra": "coulombintra.dat", "PairHop": "pairhop.dat"}
_PAIRHOP_FILES = ("pairhop.dat",)
#: the three types whose spin-orbital contraction is a SEPARATE piece of
#: code from the CoulombInter/Hund one (``uhfk.UHFk._make_ham``'s virtual
#: branch duplicates the normal-mode expressions per type), declared
#: together on inter-orbital off-site bonds so that one case covers them.
_ALL_TYPES = {"CoulombIntra": "coulombintra.dat", "Ising": "ising.dat",
              "Exchange": "exchange.dat", "PairLift": "pairlift.dat"}
_ALL_TYPES_FILES = ("ising.dat", "exchange.dat", "pairlift.dat")
#: Types that contribute NOTHING to a spin-collinear density, so that a case
#: declaring one needs a spin-mixing initial Green to be non-vacuous for it.
#: See the comment in :func:`_params`.
_NEEDS_SPIN_MIXING = frozenset({"PairLift"})

#: The COMMITTED initial Green function those cases start from, one per
#: fixture directory (the two write the same physical state in the two index
#: conventions; see :func:`_spin_mixed_seed`). A file rather than
#: ``initial_mode = "random"``: the random route goes through numpy's global
#: legacy generator, so the state both sides of the comparison start from
#: would be whatever that generator produces on the day -- and the two sides
#: are TWO SOURCE TREES, which need not agree about it. A committed file is
#: read by both, so the end-to-end identity stays exact for a reason that is
#: in the repository.
_SEED = "initial_green_spin_mixed.npz"

#: Transverse amplitude of that seed, in units of the on-site density it is
#: added to (0.25 here). MEASURED: every value from 0.02 to 0.25 sends both
#: all-types fixtures to the same broken solution, so 0.05 sits in the
#: middle of a decade-wide working range rather than on its edge.
_SEED_AMPLITUDE = 0.05

_REGENERATE_SEED_ENV = "HWAVE_REGENERATE_UHFK_SEED"

_HARTREE_ONLY = ("Hartree-only: the mean field reads the table as "
                 "sum_r J_ab(r), which the Hermitian closure makes "
                 "Hermitian and the orientation therefore leaves alone")
_ONE_ORBITAL = ("norb = 1 with a real transfer: F^rev's Hamiltonian is the "
                "complex conjugate of F's, and energy.dat's observables are "
                "conjugation-invariant")

#: The fourteen cases of the gate, keyed by label. ``moves`` says whether
#: the orientation is visible in ``energy_total`` at all -- see the module
#: docstring; ``why`` is the derivation for the ones where it is not.
CASES = {c.label: c for c in (
    Case("bond/norb2/fock=on", "norb2", 2, _BOND, _BOND_FILES,
         False, True, True, None),
    Case("bond/norb2/fock=off", "norb2", 2, _BOND, _BOND_FILES,
         False, False, False, _HARTREE_ONLY),
    Case("bond/norb2/fock=on/so", "norb2_so", 2, _BOND, _BOND_FILES,
         True, True, True, None),
    Case("bond/norb2/fock=off/so", "norb2_so", 2, _BOND, _BOND_FILES,
         True, False, False, _HARTREE_ONLY),
    Case("pairhop/norb2/fock=on", "norb2", 2, _PAIRHOP, _PAIRHOP_FILES,
         False, True, True, None),
    Case("pairhop/norb2/fock=off", "norb2", 2, _PAIRHOP, _PAIRHOP_FILES,
         False, False, True, None),
    Case("pairhop/norb2/fock=on/so", "norb2_so", 2, _PAIRHOP, _PAIRHOP_FILES,
         True, True, True, None),
    Case("pairhop/norb2/fock=off/so", "norb2_so", 2, _PAIRHOP, _PAIRHOP_FILES,
         True, False, True, None),
    Case("pairhop/norb1/fock=on", "norb1", 1, _PAIRHOP, _PAIRHOP_FILES,
         False, True, False, _ONE_ORBITAL),
    Case("pairhop/norb1/fock=off", "norb1", 1, _PAIRHOP, _PAIRHOP_FILES,
         False, False, False, _ONE_ORBITAL),
    Case("pairhop/norb1/fock=on/so", "norb1_so", 1, _PAIRHOP, _PAIRHOP_FILES,
         True, True, False, _ONE_ORBITAL),
    Case("pairhop/norb1/fock=off/so", "norb1_so", 1, _PAIRHOP, _PAIRHOP_FILES,
         True, False, False, _ONE_ORBITAL),
    Case("alltypes/norb2/fock=on", "norb2", 2, _ALL_TYPES, _ALL_TYPES_FILES,
         False, True, True, None),
    Case("alltypes/norb2/fock=on/so", "norb2_so", 2, _ALL_TYPES,
         _ALL_TYPES_FILES, True, True, True, None),
)}


def _spin_mixed_seed(nvol, ns, norb, ncond, amplitude=_SEED_AMPLITUDE):
    """The committed initial Green function of the all-types cases, in the
    layout ``uhfk.UHFk._read_green_from_data`` expects: a ``(nvol, ns, norb,
    ns, norb)`` complex array of ``G_ab(r)``, member name ``green``.

    Two terms, and nothing random:

    * the UNIFORM on-site density, ``Ncond / (nvol * ns * norb)`` on the
      diagonal of ``r = 0`` and zero at every other displacement. That is a
      normalised, featureless density -- it carries the right electron count
      and no structure for the SCF to be steered by;
    * a TRANSVERSE admixture ``amplitude * sigma_y`` on the same
      displacement, i.e. ``-i a`` on the (up, down) entry of each orbital
      and ``+i a`` on its (down, up) partner. Hermitian, so the seed is
      still a density matrix.

    Why sigma_y and not sigma_x: ``PairLift``'s spin table pairs
    ``<c^dag_up c_down>`` with itself rather than with its conjugate, so the
    term is sensitive to the PHASE of the transverse density and the
    instability it drives lies along one transverse direction. MEASURED: the
    real (sigma_x) admixture decays back to the collinear fixed point at
    several amplitudes -- PairLift then contributes 1e-20, i.e. the case is
    vacuous for it -- while the sigma_y admixture reaches the broken
    solution at every amplitude tried, in normal and spin-orbital mode
    alike.

    ``ns == 2`` is the normal mode, where spin is its own axis; ``ns == 1``
    is spin-orbital mode, where the orbital index already carries it as
    ``so = 2 * orbital + spin`` (``uhfk.UHFk._init_lattice``). The two
    fixtures get the SAME physical state written in those two conventions,
    which is why the spin-orbital twin of every case reproduces its normal
    -mode numbers."""
    green = np.zeros((nvol, ns, norb, ns, norb), dtype=np.complex128)
    filling = ncond / float(nvol * ns * norb)
    for s in range(ns):
        for a in range(norb):
            green[0, s, a, s, a] = filling
    if ns == 2:
        for a in range(norb):
            green[0, 0, a, 1, a] = -1j * amplitude
            green[0, 1, a, 0, a] = 1j * amplitude
    else:
        for a in range(norb // 2):
            green[0, 0, 2 * a, 0, 2 * a + 1] = -1j * amplitude
            green[0, 0, 2 * a + 1, 0, 2 * a] = 1j * amplitude
    return green


def _spin_off_diagonal(green, spin_orbital):
    """``max |G|`` over the SPIN-OFF-DIAGONAL entries of a
    ``(nvol, ns, norb, ns, norb)`` Green function -- the part ``PairLift``
    contracts, and the one a collinear state has none of.

    Written for both index conventions because the two fixtures use both:
    in normal mode the spin axes are separate and the block is
    ``G[:, 0, :, 1, :]`` with its partner; in spin-orbital mode spin lives
    in the orbital index as ``so = 2 * orbital + spin``, so the entries
    whose two spin bits differ have to be picked out one by one."""
    green = np.asarray(green)
    ns, norb = green.shape[1], green.shape[2]
    if not spin_orbital:
        return max(float(np.abs(green[:, 0, :, 1, :]).max()),
                   float(np.abs(green[:, 1, :, 0, :]).max()))
    worst = 0.0
    for a in range(norb):
        for b in range(norb):
            if a % 2 != b % 2:                       # different spin bit
                worst = max(worst, float(np.abs(green[:, 0, a, 0, b]).max()))
    return worst


def _seed_path(fixture):
    return os.path.join(_FIXTURES, fixture, _SEED)


def regenerate_seed():
    """Write the committed seed of both all-types fixtures from
    :func:`_spin_mixed_seed`. Deliberate only -- see the module's ``__main__``
    entry point."""
    written = []
    for fixture, ns, norb in (("norb2", 2, 2), ("norb2_so", 1, 4)):
        green = _spin_mixed_seed(nvol=4, ns=ns, norb=norb, ncond=4)
        np.savez_compressed(_seed_path(fixture), green=green)
        written.append(_seed_path(fixture))
    return written


def _params(case):
    """The input a user would write for ``case``, as the parsed dict
    ``hwave.qlms.run`` takes.

    ``flag_fock`` and ``enable_spin_orbital`` sit in the ``mode`` section,
    NOT in ``mode.param``: ``uhfk.UHFk.__init__`` passes ``info_mode`` (not
    ``info_mode["param"]``) to ``_init_mode``, so a ``flag_fock`` written
    one level down is silently ignored and every case would run with the
    default Fock term on.
    """
    interaction = {"path_to_input": ".", "Geometry": "geom.dat",
                   "Transfer": "transfer.dat"}
    interaction.update(case.files)
    green_input = {"path_to_input": "", "interaction": interaction}
    if _NEEDS_SPIN_MIXING & set(case.files):
        # PairLift reads ONLY the spin-off-diagonal part of the density (its
        # spin table is ``spin[0,0,1,1] = spin[1,1,0,0] = 1``, so both legs
        # of its w1/w2 contractions cross the two spins). UHFk's default
        # initial Green is ZERO, and the loop never leaves its
        # spin-off-diagonal block: the term would contribute exactly 0.0 for
        # ever and the case would be vacuous for it -- which is what
        # ``_assert_offsite_energy_is_live`` below measured. The committed
        # seed (:data:`_SEED`, written by :func:`_spin_mixed_seed`) puts the
        # spin mixing there, and puts it there IDENTICALLY on both sides of
        # every comparison -- the two sides are two source trees, and a seed
        # drawn from numpy's global legacy generator is a shared starting
        # point only for as long as the two trees agree about that
        # generator.
        green_input["initial"] = _SEED
    return {
        "log": {"print_level": 0, "print_step": 1000},
        "mode": {"mode": "UHFk",
                 "enable_spin_orbital": case.spin_orbital,
                 "flag_fock": case.fock,
                 "param": {"T": 0.05,
                           "Ncond": 4 if case.norb == 2 else 2,
                           "IterationMax": 1000, "EPS": 12, "Mix": 0.5,
                           "RndSeed": 1,
                           "CellShape": [4, 1, 1], "SubShape": [1, 1, 1]}},
        "file": {"input": green_input,
                 "output": {"path_to_output": "output",
                            "energy": "energy.dat", "green": "green.dat"}},
    }


def _read_energy(path):
    """``energy.dat`` as ``{key: float}`` -- the parse ``tests/test_uhf.py``
    already uses for the reference-data suite, minus its case-insensitive
    container (both sides here are written by the same writer, so the
    spellings match exactly and a difference in the KEY SET is a
    difference worth failing on)."""
    table = {}
    with open(path) as handle:
        lines = handle.read().splitlines()
    for line in lines:
        if not line.strip():
            continue
        key, _, value = line.partition("=")
        table[key.strip()] = float(value)
    return table


def _run(checkout, indir, case):
    """One UHFk run of ``case`` on the declaration in ``indir``, with
    ``checkout``'s source tree, in a COPY of ``indir`` under a temporary
    directory. Returns ``(energy, green)``.

    A subprocess, not an in-process call, for two reasons: the reference
    revision's ``hwave`` cannot be imported into the same interpreter as
    this branch's, and ``hwave.qlms.run`` configures process-wide logging
    and writes relative to the working directory.
    """
    work = tempfile.mkdtemp(prefix="hwave_uhfk_orient_")
    try:
        shutil.copytree(indir, work, dirs_exist_ok=True)
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            (os.path.join(checkout, "src"), checkout))
        proc = subprocess.run(
            [sys.executable, "-B", "-c", _RUN_SCRIPT,
             json.dumps(_params(case))],
            cwd=work, env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                "UHFk failed for case {} on {} with {}:\n{}\n{}".format(
                    case.label, indir, checkout,
                    proc.stdout[-2000:], proc.stderr[-2000:]))
        energy = _read_energy(os.path.join(work, "output", "energy.dat"))
        with np.load(os.path.join(work, "output", "green.dat.npz"),
                     allow_pickle=True) as archive:
            green = {k: np.array(archive[k]) for k in archive.files}
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return energy, green


def _reversed_dir(case):
    """``F^rev`` of ``case``'s fixture in a fresh temporary directory (the
    caller removes it).

    The writer is ``tests/test_flex_hf_scf._reversed_interaction_dir``, the
    one this branch already uses for the FLEX half of the same contract --
    one implementation of "negate the displacement of every off-site row",
    not two. Its ``copy`` argument carries ``coulombintra.dat`` across
    verbatim: an on-site file has no row whose displacement could be
    negated, so passing it as a file to REVERSE would (rightly) be refused
    as vacuous. For the cases that start from the committed seed it also
    carries :data:`_SEED`: ``F`` and ``F^rev`` are two declarations of the
    same Hamiltonian, so the two runs have to start from the same state --
    a reversed directory without the seed would fall back to UHFk's zero
    initial Green and the identity would be comparing two different
    calculations.
    """
    from tests.test_flex_hf_scf import _reversed_interaction_dir
    copy = ["geom.dat", "transfer.dat", "coulombintra.dat"]
    if _NEEDS_SPIN_MIXING & set(case.files):
        copy.append(_SEED)
    return _reversed_interaction_dir(
        os.path.join(_FIXTURES, case.fixture), case.twobody, copy=tuple(copy))


class _OrientationMixin:
    """Shared assertions, mixed into the two gates below (a plain mixin, so
    that no empty ``TestCase`` is discovered for it)."""

    def _assert_offsite_energy_is_live(self, case, energy, what):
        """The OFF-SITE interaction of ``case`` actually contributes to the
        energy this run reports.

        This is the anti-vacuity leg of every comparison below, and it has
        to be the off-site term rather than the Green function: the free
        Green function of a metal is nonzero whatever the interaction does,
        so "the Green function is not zero" would be satisfied by a run
        that never read ``coulombinter.dat`` at all. ``energy.dat`` carries
        one member PER INTERACTION TYPE (``uhfk.UHFk._calc_energy`` fills
        ``physics["Ene"][type]``, ``save_results`` prints it as
        ``Energy_<type>``), so the off-site types of the case can be asked
        for directly.

        Asserted in EVERY case, the structurally blind ones included:
        "blind" means ``F`` and ``F^rev`` give the same answer, not that
        the interaction is absent -- a blind case whose off-site term had
        silently dropped out would agree with everything and prove
        nothing."""
        for t in sorted(case.files):
            if t == "CoulombIntra":
                continue                    # on-site by construction
            key = "Energy_" + t
            self.assertIn(key, energy,
                          "{}: energy.dat carries no {} member, so the "
                          "off-site interaction is not in the energy at all"
                          .format(what, key))
            self.assertGreater(
                abs(energy[key]), 1e-6,
                "{}: the off-site {} contributes {:.3e} to the energy, i.e. "
                "nothing -- the comparison would hold for a run that never "
                "read the interaction".format(what, t, energy[key]))
            print("\nRECORDED off-site energy {}: {} = {:.6e}"
                  .format(what, key, energy[key]))

    def _assert_same_run(self, a, b, what, case):
        """Every member of two runs is EXACTLY equal: the key sets, every
        ``energy.dat`` number and every ``green.dat.npz`` array, at
        ``np.array_equal`` -- no tolerance, not even at round-off.

        The two runs do the same arithmetic on the same machine (see the
        module docstring): the oriented table of ``F`` and the reference
        revision's table of ``F^rev`` are bit-identical, so the whole SCF
        trajectory is. A round-off fallback here would be a place for a
        real difference in the mean field to hide, and there is no
        measurement it could be calibrated from.
        """
        (ea, ga), (eb, gb) = a, b
        self.assertEqual(sorted(ea), sorted(eb),
                         "{}: energy.dat key sets differ".format(what))
        self.assertEqual(sorted(ga), sorted(gb),
                         "{}: green.dat.npz member names differ".format(what))
        self._assert_offsite_energy_is_live(case, ea, what)
        for key in sorted(ea):
            self.assertEqual(
                ea[key], eb[key],
                "{}: energy.dat member {} differs by {:.3e} absolute"
                .format(what, key, abs(ea[key] - eb[key])))
        for key in sorted(ga):
            self.assertTrue(
                np.array_equal(ga[key], gb[key]),
                "{}: green.dat.npz member {} is not bit-identical"
                .format(what, key))

    def _relative_energy_gap(self, a, b):
        ea, eb = a[0], b[0]
        return abs(ea["Energy_Total"] - eb["Energy_Total"]) / max(
            abs(ea["Energy_Total"]), 1e-30)

    def _relative_green_gap(self, a, b):
        """The largest relative difference over the ``green.dat.npz``
        members, ``max|x - y| / max|x|``.

        The second half of "the two runs differ", and the more sensitive
        half: the total energy is a stationary functional of the density,
        so a mean field that moves by 1e-4 shifts the energy by only a few
        1e-6 (measured: 4.4e-05 on the Green function against 6.1e-06 on
        the energy for the same pair of runs). The energy stays the
        headline of the anti-vacuity control because it is the number a
        user compares between releases; this one keeps that control from
        resting on the smallest observable in the output.

        A difference in the member NAMES, or in a non-numeric member such
        as ``momentum_convention``, is a difference too -- infinite, so it
        can never be mistaken for agreement.
        """
        ga, gb = a[1], b[1]
        if sorted(ga) != sorted(gb):
            return float("inf")
        gap = 0.0
        for key in sorted(ga):
            x, y = np.asarray(ga[key]), np.asarray(gb[key])
            if x.dtype.kind not in "fc":
                if not np.array_equal(x, y):
                    return float("inf")
                continue
            scale = max(float(np.abs(x).max()), 1e-30)
            gap = max(gap, float(np.abs(x - y).max()) / scale)
        return gap


class TestUHFkOrientation(_OrientationMixin, unittest.TestCase):
    """G-UHFk: this tree's UHFk on ``F`` equals the reference revision's
    UHFk on ``F^rev``, member for member -- and, wherever the orientation
    is visible at all, differs from the reference revision's UHFk on ``F``.

    Skipped as a whole when the reference checkout is not usable: the
    shared rule of every develop-comparison harness on this branch
    (``tests/test_flex_second_order_compat.develop_checkout``) is that the
    reference tree must be at the NAMED revision and clean, or the
    comparison is against a tree nobody can name.
    """

    def _check(self, *labels):
        from tests.test_flex_second_order_compat import develop_checkout
        here = os.getcwd()
        reference, why = develop_checkout()
        if reference is None:
            self.skipTest(why)
        for label in labels:
            case = CASES[label]
            source = os.path.join(_FIXTURES, case.fixture)
            reversed_dir = _reversed_dir(case)
            try:
                mine = _run(here, source, case)
                theirs = _run(reference, reversed_dir, case)
                control = _run(reference, source, case)
            finally:
                shutil.rmtree(reversed_dir, ignore_errors=True)
            self._assert_same_run(mine, theirs, "case {}".format(label), case)
            gap = self._relative_energy_gap(mine, control)
            green_gap = self._relative_green_gap(mine, control)
            if case.moves:
                self.assertGreater(
                    gap, 1e-6,
                    "case {}: the reference revision's UHFk on the "
                    "UNREVERSED declaration is within {:.3e} relative of "
                    "this tree's total energy, so the identity above says "
                    "nothing".format(label, gap))
                # and the same control on the MEAN FIELD, which carries the
                # difference an order of magnitude more loudly than the
                # stationary energy does
                self.assertGreater(
                    green_gap, 1e-6,
                    "case {}: the reference revision's UHFk on the "
                    "UNREVERSED declaration returns the same Green function "
                    "to {:.3e} relative".format(label, green_gap))
            else:
                # the coinciding direction, asserted rather than assumed:
                # the reason is structural (see the module docstring), so a
                # change that makes the orientation visible here has to come
                # back through this test.
                self.assertLessEqual(
                    gap, 1e-12,
                    "case {}: the orientation moved a case that cannot see "
                    "it ({:.3e} relative in the total energy). The "
                    "derivation is: {}".format(label, gap, case.why))
                # the Green function of a blind case is allowed the SCF
                # residual and no more. The Hartree-only cases agree to the
                # last bit (identical mean field, identical trajectory);
                # the norb = 1 pair solves two CONJUGATE problems whose
                # converged densities meet only at the convergence
                # tolerance, measured 9.5e-10 -- far below anything the
                # orientation moves (4e-05) and far above bit equality.
                self.assertLessEqual(
                    green_gap, 1e-6,
                    "case {}: the orientation moved the Green function of a "
                    "case that cannot see it ({:.3e} relative). The "
                    "derivation is: {}".format(label, green_gap, case.why))

    def test_bond_fock_on(self):
        """The inter-orbital ``CoulombInter`` + ``Hund`` bond with the Fock
        term on -- the one case of this module the FAST gate runs."""
        self._check("bond/norb2/fock=on")

    @heavy
    def test_bond_variants(self):
        """The same bond declaration with the Fock term off and in
        spin-orbital mode."""
        self._check("bond/norb2/fock=off", "bond/norb2/fock=on/so",
                    "bond/norb2/fock=off/so")

    @heavy
    def test_pairhop_complex(self):
        """The complex off-site ``PairHop`` amplitude, two orbitals and one,
        Fock term on and off, normal and spin-orbital mode -- the full
        2 x 2 x 2 matrix of (norb, Fock, spin-orbital)."""
        self._check("pairhop/norb2/fock=on", "pairhop/norb2/fock=off",
                    "pairhop/norb2/fock=on/so", "pairhop/norb2/fock=off/so",
                    "pairhop/norb1/fock=on", "pairhop/norb1/fock=off",
                    "pairhop/norb1/fock=on/so", "pairhop/norb1/fock=off/so")

    @heavy
    def test_all_types_bond(self):
        """``Ising`` + ``Exchange`` + ``PairLift`` on inter-orbital off-site
        bonds, normal and spin-orbital mode.

        Every other case of this module exercises the ``CoulombInter`` /
        ``Hund`` contraction or the ``PairHop`` one. These three are
        separate expressions in ``hartree_fock.accumulate_hf`` (different
        ``spin_table`` patterns) and separate expressions again in the
        spin-orbital branch of ``uhfk.UHFk._make_ham``, which duplicates
        them; this case is what says the orientation reaches those
        duplicates too."""
        self._check("alltypes/norb2/fock=on", "alltypes/norb2/fock=on/so")


class TestFixturesAreOrientationSensitive(_OrientationMixin,
                                          unittest.TestCase):
    """The same statement WITHOUT the reference tree: this tree's UHFk on
    ``F`` against this tree's UHFk on ``F^rev``.

    Its job is to keep the gate above honest when it skips. The identity
    it pins cannot be checked here -- that needs the reference revision --
    but the PREMISE can: that these fixtures are orientation-sensitive at
    all, so that the identity is a statement about the reading of the
    declaration and not about two runs that would agree anyway. One method
    per fixture directory, so a fixture that stops moving is named by the
    failure.
    """

    def _gaps(self, label):
        """``(energy gap, Green-function gap)``, both relative -- the same
        pair of measures the gate above applies to its control."""
        case = CASES[label]
        here = os.getcwd()
        source = os.path.join(_FIXTURES, case.fixture)
        reversed_dir = _reversed_dir(case)
        try:
            mine = _run(here, source, case)
            reversed_run = _run(here, reversed_dir, case)
        finally:
            shutil.rmtree(reversed_dir, ignore_errors=True)
        # the same anti-vacuity leg the gate above applies: an off-site term
        # that contributed nothing would make every gap below meaningless,
        # the "blind" verdicts most of all
        self._assert_offsite_energy_is_live(case, mine[0], "case {}".format(label))
        return (self._relative_energy_gap(mine, reversed_run),
                self._relative_green_gap(mine, reversed_run))

    def _assert_moves(self, label):
        gap, green_gap = self._gaps(label)
        self.assertGreater(
            gap, 1e-6,
            "case {}: reversing every off-site displacement moves this "
            "tree's total energy by only {:.3e} relative, so the fixture "
            "cannot see the orientation".format(label, gap))
        self.assertGreater(
            green_gap, 1e-6,
            "case {}: reversing every off-site displacement moves this "
            "tree's Green function by only {:.3e} relative"
            .format(label, green_gap))

    def _assert_blind(self, label):
        gap, green_gap = self._gaps(label)
        self.assertLessEqual(
            gap, 1e-12,
            "case {}: reversing the declaration moved the total energy of a "
            "case that cannot see the orientation ({:.3e} relative). The "
            "derivation is: {}".format(label, gap, CASES[label].why))
        # the Green function of a blind case is allowed the SCF residual
        # and no more -- see the same assertion in the gate above.
        self.assertLessEqual(
            green_gap, 1e-6,
            "case {}: reversing the declaration moved the Green function of "
            "a case that cannot see the orientation ({:.3e} relative). The "
            "derivation is: {}".format(label, green_gap, CASES[label].why))

    def test_norb2_bond(self):
        self._assert_moves("bond/norb2/fock=on")

    def test_norb2_pairhop(self):
        self._assert_moves("pairhop/norb2/fock=on")

    def test_norb2_so_bond(self):
        self._assert_moves("bond/norb2/fock=on/so")

    def test_norb2_all_types(self):
        """The ``Ising`` + ``Exchange`` + ``PairLift`` declaration, normal
        and spin-orbital mode: three contractions the bond and PairHop
        fixtures never reach."""
        self._assert_moves("alltypes/norb2/fock=on")
        self._assert_moves("alltypes/norb2/fock=on/so")

    def test_hartree_only_is_blind(self):
        """``flag_fock = false``: the Hartree term reads the table only
        through ``sum_r J_ab(r)``, which the orientation cannot move."""
        self._assert_blind("bond/norb2/fock=off")

    def test_one_orbital_is_blind(self):
        """``norb = 1`` with a real transfer: ``F^rev``'s Hamiltonian is the
        complex conjugate of ``F``'s, and ``energy.dat`` is conjugation-
        invariant -- in normal and in spin-orbital mode alike."""
        self._assert_blind("pairhop/norb1/fock=on")
        self._assert_blind("pairhop/norb1/fock=on/so")


class TestCommittedSeed(unittest.TestCase):
    """The committed initial Green function of the all-types cases is what
    it claims to be.

    Needs no reference tree and no UHFk run, so it stands even when the gate
    above skips. Its job is to tell a FIXTURE defect apart from an
    orientation failure: if the seed ever lost its spin-off-diagonal part --
    a regeneration with the wrong amplitude, a file truncated in transit --
    ``PairLift`` would contribute exactly zero, the anti-vacuity leg of the
    gate would fail, and the failure would read like a problem with the
    reading of the declaration. This test names the real cause first."""

    _SEEDED = (("norb2", False, (4, 2, 2, 2, 2)),
               ("norb2_so", True, (4, 1, 4, 1, 4)))

    def test_the_committed_seed_breaks_spin_collinearity(self):
        for fixture, spin_orbital, shape in self._SEEDED:
            with self.subTest(fixture=fixture):
                path = _seed_path(fixture)
                self.assertTrue(os.path.exists(path),
                                "the committed seed {} is missing; it is part of "
                                "the repository -- restore it, or regenerate "
                                "deliberately with {}=1"
                                .format(path, _REGENERATE_SEED_ENV))
                with np.load(path) as archive:
                    self.assertIn("green", archive.files,
                                  "{}: UHFk reads the member named 'green' "
                                  "(uhfk.UHFk._read_green_from_data)".format(path))
                    green = np.array(archive["green"])
                self.assertEqual(green.shape, shape, path)
                self.assertEqual(green.dtype, np.complex128, path)
                off = _spin_off_diagonal(green, spin_orbital)
                self.assertGreater(
                    off, 1e-3,
                    "{}: the seed's spin-off-diagonal density is {:.3e}, i.e. it "
                    "is spin-collinear. PairLift contracts only that block, so "
                    "the all-types case would be vacuous for it -- this is a "
                    "fixture defect, not an orientation failure".format(path, off))
                # and it really is the TRANSVERSE (sigma_y) admixture the
                # generator documents, not a stray real part: a real
                # admixture decays back to the collinear fixed point
                if spin_orbital:
                    entry = green[0, 0, 0, 0, 1]
                else:
                    entry = green[0, 0, 0, 1, 0]
                self.assertAlmostEqual(entry.real, 0.0, places=15, msg=path)
                self.assertAlmostEqual(abs(entry.imag), _SEED_AMPLITUDE, places=15, msg=path)

    def test_the_seed_matches_the_generator(self):
        """The committed bytes are what :func:`_spin_mixed_seed` produces
        today: a pin on the file, so that an edit to the generator that was
        never regenerated (or a regeneration that was never committed) is a
        failure rather than a silent divergence between the recipe in the
        module and the file the runs actually read."""
        for fixture, _so, shape in self._SEEDED:
            with self.subTest(fixture=fixture):
                with np.load(_seed_path(fixture)) as archive:
                    green = np.array(archive["green"])
                want = _spin_mixed_seed(nvol=shape[0], ns=shape[1], norb=shape[2], ncond=4)
                self.assertTrue(
                    np.array_equal(green, want),
                    "{} differs from _spin_mixed_seed() by {:.3e}; regenerate it "
                    "with {}=1".format(_seed_path(fixture),
                                       np.abs(green - want).max(), _REGENERATE_SEED_ENV))

    def test_the_two_fixtures_hold_the_same_physical_state(self):
        """The normal-mode and spin-orbital seeds are the same density in
        the two index conventions (``so = 2 * orbital + spin``).

        This is what entitles the ``/so`` twin of each all-types case to
        reproduce its normal-mode numbers: if the two files drifted apart,
        the twins would be solving different problems and the agreement of
        their energies would be a coincidence."""
        with np.load(_seed_path("norb2")) as archive:
            plain = np.array(archive["green"])
        with np.load(_seed_path("norb2_so")) as archive:
            folded = np.array(archive["green"])
        norb = plain.shape[2]
        want = np.zeros(folded.shape, dtype=np.complex128)
        for a in range(norb):
            for b in range(norb):
                for s in range(2):
                    for t in range(2):
                        want[:, 0, 2 * a + s, 0, 2 * b + t] = plain[:, s, a, t, b]
        self.assertTrue(np.array_equal(folded, want),
                        "the two committed seeds differ by {:.3e} once folded into "
                        "the same convention".format(np.abs(folded - want).max()))
        self.assertGreater(np.abs(want).max(), 1e-3)          # anti-vacuity


def _seed_regeneration_requested():
    """True when :data:`_REGENERATE_SEED_ENV` asks for a rewrite of the
    committed seeds. Read only by the ``__main__`` entry point below --
    never at import time, so that collecting this module with the variable
    exported writes nothing (the same rule
    ``tests/test_flex_orientation_baseline.py`` gives its own vectors)."""
    return os.environ.get(_REGENERATE_SEED_ENV, "").strip() not in (
        "", "0", "false", "no", "off")


if __name__ == "__main__":
    if _seed_regeneration_requested():
        # the ONLY route that rewrites the committed seeds:
        #
        #     HWAVE_REGENERATE_UHFK_SEED=1 PYTHONPATH=src:. \
        #         python3 -B tests/test_uhfk_orientation.py
        #
        # It then exits and runs no test: verifying the fixture against a
        # file it has just written would say nothing.
        print("\n".join(regenerate_seed()))
    else:
        unittest.main()
