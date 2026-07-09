***********************************
Basic usage
***********************************

- Prerequisite

  H-wave requires the following programs and libraries:

  - python 3.x
  - numpy module
  - scipy module
  - requests module
  - tomli module
  - cupy module (optional; only for GPU execution with ``gpu = true`` --
    RPA/FLEX via ``[mode.param]`` and the dynamic Eliashberg solver via
    ``[eliashberg]``. Install the precompiled wheel matching your CUDA
    version, e.g. ``pip install cupy-cuda12x`` for CUDA 12.x -- see the
    `CuPy installation guide <https://docs.cupy.dev/en/stable/install.html>`_.
    Without CuPy, ``gpu = true`` falls back to the CPU path with a warning.)

  Note that `numpy.fft <https://numpy.org/doc/stable/reference/generated/numpy.fft.fft.html>`_ is used for FFT calculations in H-wave UHFk, RPA, and FLEX modes.
  The spatial FFTs of the RPA/FLEX/dynamic-Eliashberg solvers switch to
  `scipy.fft <https://docs.scipy.org/doc/scipy/reference/fft.html>`_ when the
  opt-in ``fft_workers`` parameter requests parallel FFT worker threads, and to
  cuFFT on the GPU.

- Official Page

    - `GitHub repository <https://github.com/issp-center-dev/H-wave>`_

    - `Sample/Tutorial <https://isspns-gitlab.issp.u-tokyo.ac.jp/hwave-dev/hwave-gallery>`_

- Installation

  - From PyPI:

    H-wave is available from PyPI package repository as follows:

    .. code-block:: bash

       $ pip install hwave

  - From source:

    H-wave source archive can be obtained from the release site:

    https://github.com/issp-center-dev/H-wave/releases

    The latest version is available from the development site using git:

    .. code-block:: bash

       $ git clone https://github.com/issp-center-dev/H-wave.git

    Once the source files are obtained, you can install H-wave by running the
    following command. The required libraries will also be installed at the same time.

    .. code-block:: bash

       $ cd ./H-wave
       $ pip install .

- Directory structure

    ::

      .
      |-- LICENSE
      |-- README.md
      |-- pyproject.toml
      |-- docs/
      |   |-- en/
      |   |-- ja/
      |   |-- tutorial/
      |
      |-- src/
      |   |-- qlms.py
      |   |-- hwave/
      |       |-- __init__.py
      |       |-- qlms.py
      |       |-- dos.py
      |       |-- sc.py
      |       |-- qlmsio/
      |       |   |-- __init__.py
      |       |   |-- read_input.py
      |       |   |-- read_input_k.py
      |       |   |-- wan90.py
      |       |-- solver/
      |           |-- __init__.py
      |           |-- base.py
      |           |-- uhfr.py
      |           |-- uhfk.py
      |           |-- rpa.py
      |           |-- flex.py
      |           |-- perf.py
      |-- tests/
       
- Basic usage

  #. Prepare input files

     First, you need to create input files for H-wave that are an input file that specify calculation conditions, and the definition files for the Hamiltonian.
     To generate the definition files, it will be convenient to use `StdFace library <https://github.com/issp-center-dev/StdFace>`_. 
     A brief description of these files is given in Tutorial section.
     You may consult File format sections for the details.

  #. Run

     Run the H-wave program by typing the following command in the directory where the input files are placed, and the calculation will be launched.

     .. code-block:: bash

        $ hwave input.toml

     or

     .. code-block:: bash

        $ python3 path_to_H-wave/qlms.py input.toml

     When the calculation is completed, the results will be written in the output directory.
     See File format sections for the details of the output files.

     The package also provides the post-processing tools ``hwave_dos`` (density of states)
     and ``hwave_sc`` (Eliashberg / superconducting analysis), which take the same input file.
