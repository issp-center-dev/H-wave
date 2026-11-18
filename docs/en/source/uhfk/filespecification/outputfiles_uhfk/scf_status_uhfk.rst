.. highlight:: none

.. _subsec:scf_status_uhfk.dat:

scf_status
~~~~~~~~~~

The convergence status of the self-consistent (SCF) loop: whether it
converged, the last residual, the convergence threshold and the number of
steps. The file is written by default as ``scf_status.dat`` in the output
directory, before any other result file. Its name is set by the keyword
``scf_status`` in the ``file.output`` section; ``scf_status = ""`` disables
it.

An example of the file is shown as follows.

::

    format_version = 1
    solver = UHFk
    converged = true
    residual = 8.537948717029953e-16
    eps = 1e-12
    iterations = 1
    iteration_max = 1000
    residual_definition = two_frobenius_over_size:uhfk_all_r
    note =

File format
^^^^^^^^^^^

Each line is exactly ``key = value`` followed by a newline, with one space
on each side of ``=`` and no other padding. The keys always appear in the
order below. An empty value leaves the line as ``key =`` followed by one space and the
newline (for example ``note =`` in an ordinary run).

-  ``format_version = 1``

-  ``solver = UHFk``

-  ``converged = [converged]``

-  ``residual = [residual]``

-  ``eps = [eps]``

-  ``iterations = [iterations]``

-  ``iteration_max = [iteration_max]``

-  ``residual_definition = two_frobenius_over_size:uhfk_all_r``

-  ``note = [note]``

Parameters
^^^^^^^^^^

-  ``[converged]``

   **Type :**
   ``true`` or ``false``

   **Description :**
   ``true`` when the last residual is below ``eps``, i.e. the loop stopped
   because it converged.

-  ``[residual]``

   **Type :**
   Float (Python ``repr`` of the value, full precision), or ``nan``

   **Description :**
   The residual of the last completed step,

   .. math::

      R = \frac{2\,\lVert G_{\rm new} - G_{\rm old} \rVert_F}{\mathrm{size}(G)},

   where :math:`G_{\rm old}` is the Green function a step starts from,
   :math:`G_{\rm new}` the one it computes (before mixing),
   :math:`\lVert\cdot\rVert_F` the Frobenius norm and
   :math:`\mathrm{size}(G)` the number of elements of :math:`G`. In UHFk,
   :math:`G` is the array of shape ``(nvol, ns, norb, ns, norb)``: every
   translation vector :math:`r` of the cell and every spin and orbital pair
   (in spin-orbital mode ``ns = 1`` and ``norb`` counts spin-orbitals).
   UHFr uses the same formula on its full :math:`2N \times 2N` matrix.

   Because :math:`R` divides a norm by the number of elements, it does not
   scale with the system size: the same ``EPS`` is a different criterion for
   a larger or a smaller system. Compare residuals only between runs of the
   same size.

   ``nan`` when no step was performed, or when the residual of the last
   step was not finite (see ``[note]``).

-  ``[eps]``

   **Type :**
   Float

   **Description :**
   The convergence threshold :math:`10^{-\mathrm{EPS}}` (``EPS`` in
   ``mode.param``).

-  ``[iterations]``

   **Type :**
   Integer

   **Description :**
   The number of completed SCF steps, the converged step included. It equals
   the number of lines of the ``print_check`` file, when one is configured.

-  ``[iteration_max]``

   **Type :**
   Integer

   **Description :**
   ``IterationMax`` of ``mode.param``.

-  ``[note]``

   **Type :**
   String (may be empty)

   **Description :**
   Empty for an ordinary run. Otherwise one of:

   - ``no SCF step performed``: ``IterationMax = 0``;
   - ``residual not finite after step <n>``: the residual of the last
     completed step ``<n>`` was NaN or infinite; ``residual`` is then
     ``nan`` and ``converged`` is ``false``;
   - ``IterationMax <value> treated as 0``: a negative ``IterationMax``
     (accepted only when the parameter checks are relaxed by
     ``trustme_interaction_range``) runs no step and is written as
     ``iteration_max = 0``.

Invariants
^^^^^^^^^^

-  ``eps`` is finite and positive.

-  :math:`0 \le` ``iterations`` :math:`\le` ``iteration_max``.

-  If ``iterations`` > 0, ``residual`` is finite and non-negative, or it is
   ``nan`` together with the note ``residual not finite after step <n>``.

-  If ``converged = true``, then ``iterations`` > 0 and ``residual`` <
   ``eps``.

-  If ``iterations = 0``, then ``converged = false`` and ``residual = nan``.

Runs without any SCF step
^^^^^^^^^^^^^^^^^^^^^^^^^

With ``IterationMax = 0`` no step is performed. One WARNING is logged, and
``scf_status.dat`` is the only result file written (with ``iterations =
0``, ``converged = false``, ``residual = nan`` and the note ``no SCF step
performed``); ``energy``, ``eigen``, ``green`` and the other outputs are not
written. A ``print_check`` file, if configured, is created empty. With
``scf_status = ""`` no result file is written at all, and the WARNING says
so. Earlier versions stopped with an error in this case.

The ``print_check`` log of ``log`` section (step, residual, energy,
``NCond``, ``Sz`` per step) is unchanged and stays available.

.. note::

   The key name ``scf_status`` is matched exactly, like the other keys of
   ``file.output``: ``SCF_Status`` is not the same key. Its value must be a
   string; any other value (``false``, ``0``) is refused before the
   calculation starts.
