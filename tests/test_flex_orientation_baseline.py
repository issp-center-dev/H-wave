"""Fixed-state regression vectors for the paths the interaction-row
orientation change (issue #193, spec 2026-09-16 section 2.4 / D-5)
legitimately MOVED.

Those paths -- the FLEX Hartree-Fock map, the local second-order kernel's
off-site vertex and the bond gate's mixed second-order blocks -- no longer
reproduce develop byte for byte on an input carrying inter-orbital off-site
rows, so ``tests/test_flex_second_order_compat.py``'s develop identity can
no longer carry them. What replaces the identity is this module: the
numerical state of each moved path, recorded ONCE from the reviewed
implementation and committed under ``tests/flex_orientation/``, so that a
later unintended change to any of them fails here with a number rather than
disappearing into a develop comparison that stopped covering it.

These vectors are regression PINS, not oracles. They say "this is what the
reviewed code produced", nothing more. Their CORRECTNESS is established
elsewhere, by the exact-diagonalization and oracle gates of the orientation
work:

* ``tests/test_hartree_fock_orientation.py`` and
  ``tests/test_hartree_fock_kernel.py`` -- the mean-field (Hartree-Fock)
  reading of an off-site row against ED, for every two-body type;
* ``tests/test_flex_second_order_ed_chain.py`` -- the local second-order
  kernel and the bond gate's blocks against a chain ED;
* ``tests/test_second_order_oracle.py`` and
  ``tests/test_second_order_kernel.py`` -- the off-site vertex of the
  local kernel against the pair-space oracle.

A vector that disagrees with the code therefore means one of two things: an
unintended change (fix the code) or an intended one whose ED/oracle gates
above were updated first and pass (regenerate the vectors, in the same
commit as the change, and say so).

Regeneration (deliberate, never automatic) -- this module run as a SCRIPT,
with the variable set::

    HWAVE_REGENERATE_ORIENTATION_VECTORS=1 PYTHONPATH=src:. \\
        python3 -B tests/test_flex_orientation_baseline.py

Nothing else rewrites the vectors. In particular, IMPORTING this module
with the variable set -- which any ``pytest``/``unittest`` collection of
the suite does -- writes nothing: see :func:`_regeneration_requested`. The
script prints the files it wrote and runs no test, because verifying the
code against a file it has just written proves nothing; re-run the suite
afterwards to see the new pins pass.

The fixture is the committed ``tests/rpa/input_2orb`` + ``coulombinter.dat``
at ``CellShape = [4, 4, 1]``, ``norb = 2``, ``Nmat = 16`` (the two batch
vectors record a shorter frequency axis, see :data:`_NMAT_VEC`): its
off-site content includes the INTER-ORBITAL rows ``v_12(-x) = 1`` /
``v_21(+x) = 1``, which is exactly the content the orientation change moves
(it is the identity on orbital-diagonal bonds).

Every vector file also carries a ``metadata`` member (see
:func:`_metadata`): the library versions, the platform and the revision the
numbers were recorded on. It is provenance, never a tolerance -- the
comparisons read it only to PRINT it in a failure message, so that "this
pin no longer holds" arrives together with "and it was recorded here".
"""
import json
import os
import platform
import subprocess
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"
_SHAPE = (4, 4, 1)
_NMAT = 16
#: frequency length of the two BATCH vectors (``w2``, ``mixed_w``). Both
#: kernels act strictly frequency-wise -- ``second_order.accumulate_batch``
#: reads ``chibar[l]`` alone (its ``l0`` argument only names the batch in an
#: error message) and ``flex_bond.dress_and_build_w`` dresses and assembles
#: one frequency at a time -- so a longer axis would repeat identical
#: arithmetic and write a megabyte of binary into the repository for no
#: additional coverage. The frequency-length-dependent parts of the loop
#: (the batching itself) are pinned by
#: ``tests/test_flex_second_order_compat.py::TestCompatibility::
#: test_calc_veff_general_batch_path_equals_the_dense_assembly``.
_NMAT_VEC = 4
_NORB = 2

#: the committed vectors live beside this module, NOT under ``tests/flex/``
#: (which .gitignore excludes as generated run output)
_VECTORS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "flex_orientation")

#: gate-on parameters: the Hartree-Fock term, the bond channels and the
#: dynamic bond archive (#181 Phase B) -- the combination whose develop
#: identity this module's ``smoke`` vector replaces
_GATE = {"flex_hartree_fock": True, "longitudinal_bond_channels": True,
         "longitudinal_bond_output_full": True}

_REGENERATE_ENV = "HWAVE_REGENERATE_ORIENTATION_VECTORS"

#: relative tolerance of the component vectors: the comparison replays the
#: recorded INPUT through the same arithmetic, so only run-to-run
#: reassociation (BLAS threading) can separate the two
_RTOL = 1e-12
#: relative tolerance of the smoke baseline: a converged self-consistent
#: state, so the replay reproduces it to the convergence threshold rather
#: than exactly (measured replay residual: see the module's report)
_RTOL_SMOKE = 1e-10


#: Schema version of the ``metadata`` member. Bumped when the KEY SET
#: below changes, so that a vector recorded under an older generator can be
#: told apart from one whose provenance is simply missing.
_METADATA_SCHEMA = 1

#: The keys :func:`_metadata` always writes (the per-file extras of
#: ``smoke`` come on top). Asserted on every committed file by
#: :meth:`TestOrientationVectors.test_every_vector_carries_its_provenance`.
_METADATA_KEYS = ("schema", "numpy", "scipy", "sparse_ir", "platform",
                  "python", "revision")


def _version(module_name, distribution=None):
    """The version of an OPTIONAL dependency, ``"absent"`` when it is not
    installed and ``"unknown"`` when it is installed but names no version.

    sparse-ir is the reason this is not a plain ``module.__version__``: it
    is optional here (``tests/test_flex_ir_general.py`` skips without it),
    it exports no ``__version__`` attribute at all, and "the vectors were
    recorded on a tree that did not have it" is itself provenance worth
    recording. The installed-distribution metadata is asked second, which
    is where sparse-ir's version actually lives."""
    import importlib.metadata

    try:
        module = __import__(module_name)
    except ImportError:
        return "absent"
    version = getattr(module, "__version__", None)
    if version:
        return str(version)
    try:
        return str(importlib.metadata.version(distribution or module_name))
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _metadata(**extra):
    """The provenance stamp of a vector file, as a JSON string.

    Recorded at REGENERATION time, in the process that writes the file --
    which is a script (see the module docstring), so reading the revision
    out of the repository is a plain subprocess call. A failure to read it
    is not a failure to regenerate: the revision falls back to
    ``"unknown"`` and the rest of the stamp still says which libraries and
    which platform produced the numbers.

    ``extra`` carries the per-file additions (``smoke`` records the
    iteration count and the final residual of each converged run)."""
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        revision = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"],
                                  check=True, capture_output=True,
                                  text=True).stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        revision = "unknown"
    stamp = {"schema": _METADATA_SCHEMA,
             "numpy": np.__version__,
             "scipy": _version("scipy"),
             "sparse_ir": _version("sparse_ir", "sparse-ir"),
             "platform": platform.platform(),
             "python": platform.python_version(),
             "revision": revision}
    stamp.update(extra)
    return json.dumps(stamp, sort_keys=True)


def _provenance(archive):
    """The ``metadata`` member of a loaded vector as a readable one-liner,
    for a failure message. Never raises: a file without the member is a
    finding of its own test, not a reason for another one to blow up with a
    KeyError instead of its own diagnosis."""
    if "metadata" not in getattr(archive, "files", ()):
        return "no metadata member (recorded before the provenance stamp)"
    try:
        return json.dumps(json.loads(str(archive["metadata"])), sort_keys=True)
    except (ValueError, TypeError):
        return "unreadable metadata member"


def _reader(inter=None):
    idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombInter": "coulombinter.dat"}
    idict.update(inter or {})
    return read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})


def _solver(param_extra=None):
    """A general-scheme FLEX solver on the fixture (``Nmat = 16`` unless
    ``param_extra`` says otherwise).

    Built here rather than through ``tests/test_flex_bond_gate._flex`` on
    purpose: the recipe of a committed vector has to be readable in the
    module that owns it, and it must not move when another module's helper
    changes its defaults."""
    import hwave.solver.flex as flex_mod
    par = {"T": 2.0, "filling": 0.5, "CellShape": list(_SHAPE), "SubShape": [1, 1, 1],
           "Nmat": _NMAT, "IterationMax": 3, "Mix": 0.5, "EPS": 1e-12,
           "mixing_scheme": "linear"}
    par.update(param_extra or {})
    r = _reader()
    s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                               "enable_spin_orbital": False,
                                               "calc_scheme": "general"})
    return s, r


def _random_rho_r(seed=20260916):
    """A deterministic per-spin density ``rho_ab(r)`` obeying the convention
    the Hartree-Fock map requires, ``rho_ab(r) = conj(rho_ba(-r))``: the
    real-space transform of a random HERMITIAN k-space density."""
    nvol = int(np.prod(_SHAPE))
    rng = np.random.default_rng(seed)
    rho_k = (rng.normal(size=(nvol, _NORB, _NORB))
             + 1j * rng.normal(size=(nvol, _NORB, _NORB)))
    rho_k = 0.5 * (rho_k + np.conjugate(np.swapaxes(rho_k, -1, -2)))
    return np.fft.ifftn(rho_k.reshape(*_SHAPE, _NORB, _NORB),
                        axes=(0, 1, 2)).reshape(nvol, _NORB, _NORB)


def _random_chibar(seed=20260917):
    """A deterministic ``(_NMAT_VEC, nvol, nd, nd)`` batch for the standalone
    second-order kernel. ``dense_w2`` is LINEAR in its argument and inverts
    nothing, so a random batch is both admissible and a stronger probe than
    a physical bubble: it carries every component of the vertex."""
    nvol, nd = int(np.prod(_SHAPE)), _NORB ** 2
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(_NMAT_VEC, nvol, nd, nd))
            + 1j * rng.normal(size=(_NMAT_VEC, nvol, nd, nd)))


# --- the three moved components -------------------------------------------

def _hf_map_state(rho_r=None):
    """``(rho_r, sigma_hf)``: the Hartree-Fock map of the fixture's rows."""
    from hwave.solver import flex_hf
    rho_r = _random_rho_r() if rho_r is None else np.asarray(rho_r)
    tables = flex_hf.build_flex_hf_tables(_reader().get_param("ham"), _NORB, _SHAPE)
    return rho_r, flex_hf.hf_map(rho_r, tables, _SHAPE, _NORB)


def _w2_state(chibar=None):
    """``(chibar, w2)``: the local kernel's second order at a fixed
    (deterministic random) batch."""
    from hwave.solver.second_order import dense_w2
    chibar = _random_chibar() if chibar is None else np.asarray(chibar)
    s, _ = _solver({"flex_second_order": "local", "IterationMax": 1, "Nmat": _NMAT_VEC})
    assert s._second_order_factors is not None
    assert s._second_order_factors.vpair is not None      # the off-site branch runs
    return chibar, dense_w2(chibar, s._second_order_factors)


def _bond_bubble():
    """The fixture's bare bond bubble ``chibar`` ``(nmat, nvol, ND, ND)``.

    Only used to RECORD the input of :func:`_mixed_w`; the comparison
    replays the recorded array, so a change in the bubble assembly cannot
    silently move the mixed-block vector."""
    from hwave.solver import flex_bond
    s, r = _solver({"flex_second_order": "local", "IterationMax": 1,
                    "Nmat": _NMAT_VEC, **_GATE})
    gi = r.get_param("green")
    s._phase_b_reset(gi)
    s._phase_b_preflight(gi)
    s._calc_epsilon_k({})
    nvol, nd = s.lattice.nvol, s.norb ** 2
    B = s._bond_view.n_channels
    G = s._calc_dressed_green(0.5, 0.1,
                              np.zeros((1, s.nmat, nvol, s.norb, s.norb), complex))
    with flex_bond.BondBlockStore(s.nmat, nvol, B * nd, nd, ("chibar",)) as store:
        s._phase_b_prepare_vertices()
        flex_bond.assemble_bubble(store, G, None, 0.5, s._bond_view, _SHAPE, 1)
        return np.array(store.get_freq_batch("chibar", 0, s.nmat))


def _mixed_w(chibar=None):
    """``(chibar, w_mixed_rows, w_mixed_cols)``: the bond gate's MIXED
    second-order strips ``W[:, :, :nd, nd:]`` and ``W[:, :, nd:, :nd]``,
    built from ``chibar`` through :func:`hwave.solver.flex_bond.
    dress_and_build_w` (the recipe of
    ``tests/test_flex_second_order_bond.py::_bond_w``)."""
    from hwave.solver import flex_bond
    s, r = _solver({"flex_second_order": "local", "IterationMax": 1,
                    "Nmat": _NMAT_VEC, **_GATE})
    gi = r.get_param("green")
    s._phase_b_reset(gi)
    s._phase_b_preflight(gi)
    s._calc_epsilon_k({})
    nvol, nd = s.lattice.nvol, s.norb ** 2
    B = s._bond_view.n_channels
    chibar = _bond_bubble() if chibar is None else np.asarray(chibar)
    with flex_bond.BondBlockStore(s.nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
        s._phase_b_prepare_vertices()
        store.put_freq_batch("chibar", 0, s.nmat, chibar)
        with flex_bond.BondDeviceContext.for_view(np, s._bond_S, s._bond_C, s._bond_S_on,
                                                  s._bond_C_on, s._bond_view, s.norb) as dev:
            flex_bond.dress_and_build_w(store, dev, nb=s.nmat, output_full=False,
                                        nmat=s.nmat, nvol=nvol, nd=nd, spatial_shape=_SHAPE,
                                        factors=s._second_order_factors,
                                        second_order="local")
        W = store.get_freq_batch("W", 0, s.nmat)
        return chibar, np.array(W[:, :, :nd, nd:]), np.array(W[:, :, nd:, :nd])


def _smoke_state(second_order):
    """``(sigma, green, run)`` of a gate-on run: 4x4x1, ``Nmat = 16``, run to
    CONVERGENCE (Anderson mixing, ``EPS = 1e-12``, reached in 16-17
    iterations and 0.4 s; ``IterationMax`` is a ceiling, not the recipe).

    The end-to-end pin: the state the dropped develop identity case used to
    cover (``"takimoto"`` with the bond gate), plus the production default
    (``"local"``). A converged state is pinned rather than a fixed number of
    iterations because a fixed point is defined by the equations alone --
    the mid-trajectory state of a non-converged run carries the mixer's
    history as well, which makes it both noisier under re-association and
    less meaningful to compare. It also keeps the run free of the
    non-convergence and linear-mixing warnings.

    The solver refuses to be reused, so the convergence is asserted here:
    a pin recorded from a run that silently stopped at the ceiling would
    be a different quantity from the one this function documents.

    ``run`` is the third element: the iteration count and final residual
    that convergence check produced. The generator stamps them into the
    ``smoke`` file's provenance -- "converged" is a property of the run, and
    a later regeneration that needed four times the iterations to reach the
    same fixed point is worth seeing."""
    s, r = _solver({"flex_second_order": second_order, "IterationMax": 200,
                    "mixing_scheme": "anderson", "EPS": 1e-12, **_GATE})
    gi = r.get_param("green")
    with tempfile.TemporaryDirectory() as out:
        s.solve(gi, out)
    if not s.scf_converged:
        raise AssertionError(
            "the smoke baseline did not converge in {} iterations (residual "
            "{:.3e}); the vector would pin a mixer trajectory, not a fixed "
            "point".format(s.scf_iterations, s.scf_sigma_residual))
    run = {"iterations_" + second_order: int(s.scf_iterations),
           "residual_" + second_order: float(s.scf_sigma_residual)}
    return np.asarray(gi["sigma"]), np.asarray(gi["green"]), run


# --- the committed files ---------------------------------------------------

def _path(name):
    return os.path.join(_VECTORS, name + ".npz")


def regenerate():
    """Write every vector file from the CURRENT code. Deliberate only: run
    it when an intended change has already been re-gated by the ED/oracle
    modules named in the module docstring."""
    os.makedirs(_VECTORS, exist_ok=True)
    rho_r, sigma_hf = _hf_map_state()
    np.savez_compressed(_path("hf_map"), rho_r=rho_r, sigma_hf=sigma_hf,
                        metadata=_metadata())
    chibar, w2 = _w2_state()
    np.savez_compressed(_path("w2"), chibar=chibar, w2=w2, metadata=_metadata())
    cb, rows, cols = _mixed_w()
    np.savez_compressed(_path("mixed_w"), chibar=cb, w_mixed_rows=rows, w_mixed_cols=cols,
                        metadata=_metadata())
    payload, run = {}, {}
    for so in ("local", "takimoto"):
        sigma, green, info = _smoke_state(so)
        payload["sigma_" + so], payload["green_" + so] = sigma, green
        run.update(info)
    np.savez_compressed(_path("smoke"), metadata=_metadata(**run), **payload)
    return sorted(os.listdir(_VECTORS))


def _regeneration_requested():
    """True when ``HWAVE_REGENERATE_ORIENTATION_VECTORS`` asks for a rewrite.

    Read only by the ``__main__`` entry point below -- never at import
    time. Importing this module with the variable set (which every
    ``pytest``/``unittest`` collection of the whole suite would do once an
    operator exported it) must not overwrite the committed vectors: a
    rewrite is a deliberate, single act, and a regeneration that happened
    as a side effect of collecting an unrelated test would replace the pins
    with whatever the working tree currently produces and then "verify"
    the code against it."""
    return os.environ.get(_REGENERATE_ENV, "").strip() not in ("", "0", "false", "no", "off")


class TestOrientationVectors(unittest.TestCase):
    """The moved paths against their committed fixed state.

    These tests take no skip at all. A MISSING vector file is a broken
    checkout, not a configuration: the files are committed, so the only
    honest verdict for an absent one is a failure."""

    def _load(self, name):
        p = _path(name)
        if not os.path.exists(p):
            self.fail(
                "the committed fixed-state vector {} is missing; it is part of "
                "the repository -- restore it, or regenerate deliberately with "
                "{}=1".format(p, _REGENERATE_ENV))
        return np.load(p)

    def _compare(self, got, want, rtol, what, archive=None):
        """``got`` against the recorded ``want``.

        ``archive`` is the loaded vector file: its provenance stamp goes
        into the failure message, so that a pin that stops holding says in
        the same breath WHERE it was recorded -- which numpy, which
        platform, which revision. It is a diagnostic only; no tolerance
        here reads it."""
        want = np.asarray(want)
        scale = float(np.abs(want).max())
        where = "" if archive is None else " [recorded on: {}]".format(_provenance(archive))
        self.assertGreater(scale, 1e-8,
                           "the recorded vector {} is empty{}".format(what, where))
        self.assertEqual(np.asarray(got).shape, want.shape, what + where)
        np.testing.assert_allclose(got, want, rtol=0, atol=rtol * scale,
                                   err_msg=what + where)

    def test_hf_map_matches_the_committed_vector(self):
        """The FLEX Hartree-Fock map at the recorded density."""
        d = self._load("hf_map")
        _rho, sigma = _hf_map_state(d["rho_r"])
        self._compare(sigma, d["sigma_hf"], _RTOL, "hf_map/sigma_hf", d)

    def test_w2_matches_the_committed_vector(self):
        """The local kernel's second order at the recorded random batch."""
        d = self._load("w2")
        _cb, w2 = _w2_state(d["chibar"])
        self._compare(w2, d["w2"], _RTOL, "w2/w2", d)

    def test_mixed_bond_w_matches_the_committed_vector(self):
        """The bond gate's mixed (channel-0 x bond) blocks of W at the
        recorded bond bubble -- the blocks the pair permutation of spec
        2026-09-16 R3 acts on."""
        d = self._load("mixed_w")
        _cb, rows, cols = _mixed_w(d["chibar"])
        self._compare(rows, d["w_mixed_rows"], _RTOL, "mixed_w/w_mixed_rows", d)
        self._compare(cols, d["w_mixed_cols"], _RTOL, "mixed_w/w_mixed_cols", d)

    def test_every_vector_carries_its_provenance(self):
        """Every committed vector has a ``metadata`` member, it is JSON, and
        it carries :data:`_METADATA_KEYS` at :data:`_METADATA_SCHEMA`; the
        ``smoke`` file also carries the iteration count and final residual
        of each converged run.

        A vector without provenance is still a valid pin -- which is
        exactly why this has to be a test of its own: nothing else in the
        module would notice its absence, and a regeneration that quietly
        stopped stamping the files would leave the next reader of a failure
        with no way to tell which tree produced the numbers."""
        for name in ("hf_map", "w2", "mixed_w", "smoke"):
            with self.subTest(vector=name):
                d = self._load(name)
                self.assertIn("metadata", d.files,
                              "the committed vector {}.npz carries no provenance "
                              "stamp".format(name))
                stamp = json.loads(str(d["metadata"]))
                for key in _METADATA_KEYS:
                    self.assertIn(key, stamp, (name, key))
                    self.assertTrue(str(stamp[key]),
                                    "{}.npz records an empty {}".format(name, key))
                self.assertEqual(stamp["schema"], _METADATA_SCHEMA, name)
                # the versions really are versions and not the fallback of a
                # generator that could not read anything at all
                self.assertNotEqual(stamp["numpy"], "absent", name)
                if name == "smoke":
                    for so in ("local", "takimoto"):
                        self.assertGreater(stamp["iterations_" + so], 0, so)
                        self.assertGreater(stamp["residual_" + so], 0.0, so)
                        self.assertLess(stamp["residual_" + so], 1e-11, so)

    def test_gate_on_smoke_matches_the_committed_vector(self):
        """End to end: sigma and green of a CONVERGED gate-on run under both
        kernels (0.8 s measured for the pair: below the 5 s opt-in rule of
        ``tests/heavy_tests.py``, so it stays in the fast gate).

        The ``"takimoto"`` half is the case that left
        ``tests/test_flex_second_order_compat.py``'s develop identity list
        (the bond gate reads the off-site rows through the Hartree-Fock
        term, so it no longer reproduces develop); the ``"local"`` half is
        the production default."""
        d = self._load("smoke")
        for so in ("local", "takimoto"):
            with self.subTest(second_order=so):
                sigma, green, _run = _smoke_state(so)
                self._compare(sigma, d["sigma_" + so], _RTOL_SMOKE, "smoke/sigma_" + so, d)
                self._compare(green, d["green_" + so], _RTOL_SMOKE, "smoke/green_" + so, d)


if __name__ == "__main__":
    if _regeneration_requested():
        # the ONLY route that rewrites the vectors: this module run as a
        # script, with the variable set. It then exits -- verifying the code
        # against the file it has just written would say nothing.
        print("\n".join(regenerate()))
    else:
        unittest.main()
