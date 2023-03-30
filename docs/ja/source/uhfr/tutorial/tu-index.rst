==================
チュートリアル
==================

UHFrでは、入力ファイルとして

#. 環境設定入力ファイル
#. Hamiltonian作成用ファイル
#. 出力結果指定用ファイル

を用意した後、計算を行います。以下では、 ``docs/tutorial/Hubbard/UHFr`` ディレクトリにあるサンプルを例にチュートリアルを実施します。
なお、相互作用定義ファイルは StdFace ライブラリを用いて生成することもできます。詳細は :ref:`Ch:StdFace` の章をご覧ください。

環境設定入力ファイルの作成
------------------------------------------

環境設定入力ファイルでは、入出力を制御する情報を記載します。
``docs/tutorial/Hubbard/UHFr`` ディレクトリ内に ``input.toml`` というファイルがありますが、これが環境設定入力ファイルになります。
以下、ファイルの中身を記載します。

::

    [log]
    print_level = 1
    print_step = 20
    [mode]
    mode = "UHFr"
    [mode.param]
    Nsite = 8
    2Sz = 0
    Ncond = 8
    IterationMax = 1000
    EPS = 8
    RndSeed = 123456789
    T = 0.0
    [file]
    [file.input]
    path_to_input = ""
    OneBodyG = "greenone.def"
    [file.input.interaction]
    Trans = "trans.def"
    CoulombIntra = "coulombintra.def"
    [file.output]
    path_to_output = "output"
    energy = "energy.dat"
    eigen = "eigen"
    green = "green.dat"

このファイルはtoml形式で記述します。

``log`` セクションに ``print_level`` で標準出力のレベル、 ``print_step`` でログファイルに出力するステップ間隔を指定します。

``mode`` セクションに実行モードおよび基本パラメータを指定します。

``file.input`` セクションに入力ファイルが格納されているディレクトリ ``path_to_input`` 、出力したい一体グリーン関数が定義されたファイル  ``OneBodyG``  及び初期配置 ``Initial`` を指定します。
``OneBodyG`` を指定しない場合にはグリーン関数の出力がされません。また、``Initial`` を指定しない場合には初期配置はランダムな配置が設定されます。

``file.input.interaction`` セクションにHamiltonianを作成するための入力ファイルを指定します。

``file.output`` セクションには出力ファイルを格納するディレクトリ ``path_to_output`` を指定します。
また、エネルギーの値を出力するファイル名 ``energy`` 、ハミルトニアンの固有値を出力するファイル名 ``eigen`` 、一体グリーン関数の出力ファイル名 ``green`` を指定します。これらのキーワードがない場合には情報は出力されません。

詳細については :ref:`ファイルフォーマット<Ch:Config_UHFR>` の章をご覧ください。

Hamiltonianの指定
----------------------------------

基本パラメータを設定した後は、Hamiltonianを構築するためのファイルを作成します。

**Transfer部の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Trans`` でひも付けられるファイル(ここではtrans.def)で電子系のTransferに相当するHamiltonian

.. math::

   \mathcal{H} = -\sum_{ij\sigma_1\sigma_2}
   t_{ij\sigma_1\sigma_2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom\dagger}.
   
を指定します。ファイルの中身は下記の通りです。

::

    ========================
    NTransfer      64
    ========================
    ========i_j_s_tijs======
    ========================
        4     0     0     0         1.000000000000000         0.000000000000000
        0     0     4     0         1.000000000000000        -0.000000000000000
        4     1     0     1         1.000000000000000         0.000000000000000
        0     1     4     1         1.000000000000000        -0.000000000000000
        2     0     0     0         1.000000000000000         0.000000000000000
        0     0     2     0         1.000000000000000        -0.000000000000000
        2     1     0     1         1.000000000000000         0.000000000000000
        0     1     2     1         1.000000000000000        -0.000000000000000
    ...

 
Transファイルの詳細はセクション :ref:`Subsec:Trans` をご覧ください。

**二体相互作用部の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

このチュートリアルの例では ``CoulombIntra`` でひも付けられるファイル(ここではcoulombintra.def)で電子系の二体相互作用部に相当するHamiltonian

.. math::

   \mathcal{H} = \sum_{i} U_i n_{i\uparrow}n_{i\downarrow}.

を指定します。ファイルの中身は下記の通りです。

::

    =============================================
    NCoulombIntra          8
    =============================================
    ================== CoulombIntra ================
    =============================================
        0         8.000000000000000
        1         8.000000000000000
        2         8.000000000000000
        3         8.000000000000000
        4         8.000000000000000
     ...

  
なお、 ``CoulombIntra`` 以外にも、Hamiltonianを簡易的に記載するための各種ファイル形式に対応しています。
詳細はセクション :ref:`Subsec:interall` - :ref:`Subsec:pairlift` をご覧ください。

出力ファイルの指定
-------------------------

一体Green関数の計算する成分を、``OneBodyG`` でひも付けられるファイルで指定します。

**一体Green関数の計算対象の指定**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``OneBodyG`` でひも付けられるファイル(ここではgreenone.def)で計算する一体Green関数  :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2} \rangle` の成分を指定します。ファイルの中身は下記の通りです

::

    ===============================
    NCisAjs         16
    ===============================
    ======== Green functions ======
    ===============================
        0     0     0     0
        0     0     1     0
        0     0     2     0
        0     0     3     0
        0     0     4     0
     ...

一体Green関数計算対象成分の指定に関するファイル入力形式の詳細はセクション :ref:`Subsec:onebodyg` をご覧ください。

計算の実行
--------------------------

全ての入力ファイルが準備できた後、計算実行します。
環境設定入力ファイル(ここでは ``input.toml`` )を引数とし、ターミナルからH-waveを実行します。

.. code-block:: bash

    $ hwave input.toml

計算が開始されると以下のようなログが出力されます。

::

    2022-12-01 09:37:30,114 INFO qlms: Read def files
    2022-12-01 09:37:30,116 INFO qlms: Get Hamiltonian information
    2022-12-01 09:37:30,116 INFO qlms: Get Green function information
    2022-12-01 09:37:30,116 INFO qlms.uhfr: Show input parameters
      Nsite               : 8
      Ncond               : 8
      2Sz                 : 0
      Mix                 : 0.5
      EPS                 : 1e-08
      IterationMax        : 1000
      RndSeed             : 123456789
      T                   : 0.0
      ene_cutoff          : 100.0
      threshold           : 1e-12
    2022-12-01 09:37:30,117 INFO qlms: Start UHF calculation
    2022-12-01 09:37:30,117 INFO qlms.uhfr: Set Initial Green's functions
    2022-12-01 09:37:30,117 INFO qlms.uhfr: Initialize green function by random numbers
    2022-12-01 09:37:30,117 INFO qlms.uhfr: Start UHFr calculations
    2022-12-01 09:37:30,117 INFO qlms.uhfr: step, rest, energy, NCond, Sz
    2022-12-01 09:37:30,119 INFO qlms.uhfr: 0, 0.022144468, -27.16081+0j, 8, -7.425e-16
    2022-12-01 09:37:30,134 INFO qlms.uhfr: 20, 1.2083848e-05, -3.399532+0j, 8, -1.055e-15
    2022-12-01 09:37:30,145 INFO qlms.uhfr: UHFr calculation is succeeded: rest=5.7552848630056134e-09, eps=1e-08.
    2022-12-01 09:37:30,145 INFO qlms: Save calculation results.
    2022-12-01 09:37:30,146 INFO qlms: All procedures are finished.
    --------------------------------------------------------------------------------
    Statistics
      function                         :  total elapsed  : average elapsed : ncalls
    --------------------------------------------------------------------------------
      hwave.solver.uhfr.__init__       :      0.357 msec :      0.357 msec :      1
      hwave.solver.uhfr._initial_G     :      0.090 msec :      0.090 msec :      1
      hwave.solver.uhfr._makeham_const :      0.839 msec :      0.839 msec :      1
      hwave.solver.uhfr._makeham_mat   :      0.309 msec :      0.309 msec :      1
      hwave.solver.uhfr._makeham       :      6.001 msec :      0.176 msec :     34
      hwave.solver.uhfr._diag          :      2.468 msec :      0.073 msec :     34
      hwave.solver.uhfr._green         :      3.107 msec :      0.091 msec :     34
      hwave.solver.uhfr._calc_energy   :      1.990 msec :      0.059 msec :     34
      hwave.solver.uhfr._calc_phys     :     12.929 msec :      0.380 msec :     34
      hwave.solver.uhfr.solve          :     28.290 msec :     28.290 msec :      1
      hwave.solver.uhfr.save_results   :      0.852 msec :      0.852 msec :      1
    --------------------------------------------------------------------------------
		
入力ファイル読み込みに関するログが出力されたあと、UHF計算の計算過程に関する情報が出力されます。
出力ファイルは ``input.toml`` の ``file.output`` セクションでの設定にしたがい、
``output`` ディレクトリに 固有値が記載された ``energy.dat`` ,
固有ベクトルが記載された ``spin-down_eigen.npz``, ``spin-up_eigen.npz``,
一体グリーン関数の値が記載された ``green.dat`` ファイルが出力されます。
出力ファイルの詳細については :ref:`ファイルフォーマット<Sec:outputfile>` の章をご覧ください。

