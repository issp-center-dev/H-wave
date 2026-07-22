"""``run.sh`` must scrub loader-injection environment variables.

LD_LIBRARY_PATH, LD_PRELOAD, and LD_AUDIT must be removed so
the native solver child processes cannot be tampered with via the
caller's shell environment.

This test does not invoke vmcdry.out / vmc.out / UHF (those need MPI +
built binaries + significant runtime); instead it sources the
LD_LIBRARY_PATH sanitization block from run.sh in a subshell that
exports hostile values for each of the three vars, then dumps the
resulting env. A subshell trace is sufficient because the sanitization
is a top-level assignment that runs BEFORE any child solver is
launched: if the shell's own env is clean after the block, so will
every child process's env be, courtesy of POSIX inheritance rules.
"""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path


def shlex_quote(s):
    return shlex.quote(s)


_REPO = Path(__file__).resolve().parents[1]
_RUN_SH = _REPO / "tests" / "validation" / "uhfk_mvmc_pairproduct" / "run.sh"


def _extract_sanitize_block():
    """Return the top-of-file sanitization block from run.sh as a
    standalone bash snippet, up to and including the run_native
    function definition. The end marker is `run_native() { ... }`, so
    if a future refactor
    moves the sanitize block below WORK setup or below any solver
    invocation, the extractor will still terminate but the parallel
    structural test (see test_sanitize_block_precedes_native_solver
    _invocations) will catch it."""
    text = _RUN_SH.read_text()
    lines = []
    seen_run_native_open = False
    for line in text.splitlines():
        lines.append(line)
        stripped = line.strip()
        # Terminate when the run_native function's closing brace is
        # seen (structurally: on a line beginning with '}' after the
        # opening 'run_native() {').
        if stripped == "run_native() {":
            seen_run_native_open = True
        elif seen_run_native_open and stripped == "}":
            break
    else:
        raise RuntimeError(
            "run.sh did not contain 'run_native() { ... }'; "
            "the sanitization wrapper is missing."
        )
    return "\n".join(lines) + "\n"


_GUARDED_VARS = ("LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT")


def _source_block_and_dump_env(snippet, hostile_env):
    """Write snippet + KEEP/UNSET dump to a tmp script, run it under
    bash with `hostile_env` merged into the environment.

    Assert returncode == 0 and return stdout only after that check. A
    non-zero return
    (e.g., snippet aborts on `set -e` before reaching the unsets)
    used to be silently ignored, making the test spuriously pass on
    a moved-block regression."""
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "sanitize_and_dump.sh"
        script.write_text(
            snippet + "\n"
            "for var in LD_LIBRARY_PATH LD_PRELOAD LD_AUDIT; do\n"
            "  if [[ -n ${!var+x} ]]; then\n"
            "    echo \"KEEP $var=${!var}\"\n"
            "  else\n"
            "    echo \"UNSET $var\"\n"
            "  fi\n"
            "done\n"
        )
        env = os.environ.copy()
        for v in _GUARDED_VARS + ("MVMC_LD_LIBRARY_PATH",):
            env.pop(v, None)
        env.update(hostile_env)
        res = subprocess.run(
            ["bash", str(script)],
            capture_output=True, text=True, env=env, cwd=tmp,
        )
        assert res.returncode == 0, (
            f"sanitize snippet exited with {res.returncode}; "
            f"stderr={res.stderr!r} stdout={res.stdout!r}"
        )
        return res.stdout


def _parse_env_dump(output):
    """Return {var: value_or_None} for every guarded var found in the
    KEEP/UNSET dump. Assert every guarded variable appears in the dump
    so a truncated snippet
    (which does not print the dump lines) cannot silently pass a
    var-is-absent check via `state.get(...) is None`."""
    state = {}
    for line in output.splitlines():
        line = line.strip()
        if line.startswith("KEEP "):
            _, kv = line.split(" ", 1)
            k, v = kv.split("=", 1)
            state[k] = v
        elif line.startswith("UNSET "):
            _, k = line.split(" ", 1)
            state[k] = None
    missing = [v for v in _GUARDED_VARS if v not in state]
    assert not missing, (
        f"dump missing guarded vars {missing}; snippet did not reach "
        f"the dump block. Full stdout: {output!r}"
    )
    return state


def test_ld_preload_is_unset_after_sanitization():
    """A hostile LD_PRELOAD supplied in the caller env MUST be unset
    after the run.sh sanitization block runs."""
    snippet = _extract_sanitize_block()
    dump = _source_block_and_dump_env(
        snippet, hostile_env={"LD_PRELOAD": "/tmp/attacker.so"},
    )
    state = _parse_env_dump(dump)
    assert state.get("LD_PRELOAD") is None, (
        f"LD_PRELOAD survived sanitization: {state.get('LD_PRELOAD')!r}"
    )


def test_ld_audit_is_unset_after_sanitization():
    """A hostile LD_AUDIT supplied in the caller env MUST be unset."""
    snippet = _extract_sanitize_block()
    dump = _source_block_and_dump_env(
        snippet, hostile_env={"LD_AUDIT": "/tmp/attacker_audit.so"},
    )
    state = _parse_env_dump(dump)
    assert state.get("LD_AUDIT") is None, (
        f"LD_AUDIT survived sanitization: {state.get('LD_AUDIT')!r}"
    )


def test_ld_library_path_is_unset_when_no_mvmc_path_available():
    """Without a validated MVMC_LD_LIBRARY_PATH, run.sh must unset
    LD_LIBRARY_PATH entirely (not preserve the caller's dirs)."""
    snippet = _extract_sanitize_block()
    dump = _source_block_and_dump_env(
        snippet, hostile_env={
            "LD_LIBRARY_PATH": "/tmp/attacker_lib_dir",
            # Explicitly no MVMC_LD_LIBRARY_PATH to force the fallback.
            "HOME": "/nonexistent",  # kills the auto-detect fallback
        },
    )
    state = _parse_env_dump(dump)
    assert state.get("LD_LIBRARY_PATH") is None, (
        f"LD_LIBRARY_PATH survived sanitization: "
        f"{state.get('LD_LIBRARY_PATH')!r}"
    )


def test_sanitize_block_precedes_native_solver_invocations():
    """Enforce that the sanitization endpoint appears before solvers.

    The `run_native() { ... }` close brace must appear
    in run.sh BEFORE any `run_native` invocation and BEFORE any native
    solver call (``${VMCDRY}`` / ``${VMC}`` / ``${UHF}``). Otherwise a
    future refactor that moves the block below WORK setup or below any
    solver call would silently break the trust boundary, and the
    per-var pass/unset tests would still pass because the extractor
    would still find the sanitize sequence."""
    lines = _RUN_SH.read_text().splitlines()
    sanitize_end = None
    seen_run_native_open = False
    for idx, line in enumerate(lines, start=1):
        s = line.strip()
        if s == "run_native() {":
            seen_run_native_open = True
        elif seen_run_native_open and s == "}":
            sanitize_end = idx
            break
    assert sanitize_end is not None, (
        "run.sh sanitize endpoint (`run_native` closing brace) not "
        "found; run_native wrapper missing entirely."
    )
    # Any native solver reference or run_native call MUST be at line
    # sanitize_end + 1 or later. Excludes comments (# ...).
    solver_needles = ('"${VMCDRY}"', '"${VMC}"', '"${UHF}"',
                      'run_native ')
    for idx, line in enumerate(lines, start=1):
        if idx <= sanitize_end:
            continue
        # ok, this line is after the sanitize block
        break
    for idx, line in enumerate(lines, start=1):
        if idx > sanitize_end:
            continue
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        for needle in solver_needles:
            assert needle not in line, (
                f"run.sh:{idx}: {needle!r} appears BEFORE sanitize "
                f"endpoint at line {sanitize_end}. A native solver "
                "call in that position runs with the caller's "
                "unsanitized loader env."
            )


def test_run_native_uses_absolute_env_and_resists_shadowing():
    """An unqualified `env` in the run_native wrapper can be shadowed by:
      (a) a caller-exported bash function via BASH_FUNC_env%%, or
      (b) a PATH-poisoned fake env binary earlier in PATH.
    Both defeat the per-command LD sanitization. Fix: run_native calls
    ``/usr/bin/env`` directly via the readonly ``_RUN_NATIVE_ENV``
    variable and issues ``unset -f env`` inside the sanitize block.
    This test exports a hostile ``env`` bash function to bash's
    function environment, then invokes ``run_native /bin/env-checker``
    (a tiny stub script that prints its own LD_PRELOAD / LD_AUDIT env
    state) and asserts the stub sees LD_PRELOAD / LD_AUDIT unset."""
    snippet = _extract_sanitize_block()
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        # Stub env-checker: dumps LD_PRELOAD / LD_AUDIT / LD_LIBRARY_PATH
        # so we can verify the child process really saw them unset.
        checker = tmp_path / "env_checker.sh"
        checker.write_text(
            "#!/bin/bash\n"
            "for v in LD_PRELOAD LD_AUDIT LD_LIBRARY_PATH; do\n"
            "  if [[ -n ${!v+x} ]]; then\n"
            "    echo \"CHILD_KEEP $v=${!v}\"\n"
            "  else\n"
            "    echo \"CHILD_UNSET $v\"\n"
            "  fi\n"
            "done\n"
        )
        checker.chmod(0o755)
        # Hostile env function bound via bash's exported-function env
        # var. If run_native called `env` unqualified, this shadow
        # would win and the LD vars would leak into the child.
        hostile_env_fn = (
            "() { echo \"[HOSTILE ENV SHADOW: LD passthrough]\"; "
            "exec \"$@\"; }"
        )
        script = tmp_path / "harness.sh"
        script.write_text(
            snippet + "\n"
            f"run_native {shlex_quote(str(checker))}\n"
        )
        env = os.environ.copy()
        for v in _GUARDED_VARS + ("MVMC_LD_LIBRARY_PATH",):
            env.pop(v, None)
        # Poison every var we can from the caller side.
        env["LD_PRELOAD"] = "/tmp/attacker.so"
        env["LD_AUDIT"] = "/tmp/attacker_audit.so"
        env["LD_LIBRARY_PATH"] = "/tmp/attacker_libs"
        # Kill auto-detect fallback (which would otherwise set
        # LD_LIBRARY_PATH to the miniconda lib on this host).
        env["HOME"] = "/nonexistent-home-for-loader-env-shadow-test"
        # Bash-4+ exported function encoding.
        env["BASH_FUNC_env%%"] = hostile_env_fn
        res = subprocess.run(
            ["bash", str(script)],
            capture_output=True, text=True, env=env, cwd=tmp,
        )
        assert res.returncode == 0, (
            f"harness exited {res.returncode}; "
            f"stderr={res.stderr!r} stdout={res.stdout!r}"
        )
        assert "HOSTILE ENV SHADOW" not in res.stdout, (
            "hostile env function was invoked; run_native was "
            "shadowed. stdout:\n" + res.stdout
        )
        assert "CHILD_UNSET LD_PRELOAD" in res.stdout, (
            "child process saw LD_PRELOAD kept; shadowing worked. "
            "stdout:\n" + res.stdout
        )
        assert "CHILD_UNSET LD_AUDIT" in res.stdout, (
            "child process saw LD_AUDIT kept. stdout:\n" + res.stdout
        )
        assert "CHILD_UNSET LD_LIBRARY_PATH" in res.stdout, (
            "child process saw LD_LIBRARY_PATH kept. stdout:\n"
            + res.stdout
        )


def test_ld_library_path_is_replaced_by_validated_mvmc_path():
    """When MVMC_LD_LIBRARY_PATH is a validated path, the caller's
    LD_LIBRARY_PATH must be REPLACED (not extended)."""
    snippet = _extract_sanitize_block()
    # Use a real trusted directory owned by the current user.
    tmp_lib = Path(tempfile.mkdtemp())
    try:
        # Make it look like a lib dir with libopenblas.so.0.
        (tmp_lib / "libopenblas.so.0").write_bytes(b"stub")
        os.chmod(tmp_lib, 0o755)
        dump = _source_block_and_dump_env(
            snippet, hostile_env={
                "LD_LIBRARY_PATH": "/tmp/attacker_lib_dir",
                "MVMC_LD_LIBRARY_PATH": str(tmp_lib),
                "HOME": "/nonexistent",
            },
        )
        state = _parse_env_dump(dump)
        # Must be present and equal to the accepted MVMC path (NOT
        # appended to the attacker dir).
        assert state.get("LD_LIBRARY_PATH") == str(tmp_lib), (
            f"expected LD_LIBRARY_PATH == {tmp_lib!s}; "
            f"got {state.get('LD_LIBRARY_PATH')!r}"
        )
    finally:
        shutil.rmtree(tmp_lib, ignore_errors=True)
