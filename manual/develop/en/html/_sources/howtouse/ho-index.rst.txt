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

- Official Page

  `GitHub repository <https://github.com/issp-center-dev/H-wave>`_

- Installation

  - From PyPI:

    H-wave is available from PyPI package repository as follows:

    .. code-block:: bash

       $ pip install h-wave

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
      |       |-- qlmsio/
      |       |   |-- __init__.py
      |       |   |-- read_input.py
      |       |   |-- read_input_k.py
      |       |   |-- wan90.py
      |       |-- solver/
      |           |-- __init__.py
      |           |-- base.py
      |           |-- uhf.py
      |           |-- uhfk.py
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
