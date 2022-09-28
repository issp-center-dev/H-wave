波数空間版UHFの使いかた
================================

H-waveを波数空間モードで実行するには、入力ファイルとして

1. 入力パラメータファイル
2. 相互作用定義ファイル

を用意した後、プログラムを実行します。2.は RESPACK 等の外部プログラムの出力を利用する他、
外部ツールで生成できるようにする予定です。

入力パラメータファイルの作成
--------------------------------

入力パラメータファイルには、基本パラメータの指定と入出力を制御する情報を記述します。
ディレクトリ内に ``input.toml`` というファイルがありますが、これが入力パラメータファイルになります。
以下、ファイルの内容を記載します。

::

   [log]
     print_level = 1
   [mode]
     mode = "UHFk"
   [mode.param]
     # 2Sz = 0
     Ncond = 1280
     IterationMax = 20
     EPS = 8
     Mix = 0.5
     RndSeed = 123456789
     # ene_cutoff = 1.0e+2
     T = 1.0
     CellShape = [ 8, 8, 8 ]
   [file]
   [file.input]
     path_to_input = ""
     # initial = "green_init.dat.npz"
   [file.input.interaction]
     path_to_input = "../dir-model"
     Geometry = "zvo_geom.dat"
     Transfer = "zvo_hr.dat"
     Hund = "zvo_jr.dat"
     # Coulomb = "zvo_ur.dat"
     # CoulombIntra = "coulomb_intra.dat"
     # CoulombInter = "coulomb_inter.dat"
     # Ising = "ising.dat"
     # Exchange = "exchange.dat"
     # PairLift = "pairlift.dat"
     # PairHop = "pairhop.dat"
   [file.output]
     path_to_output = "output"
     energy = "energy.dat"
     eigen = "eigen.dat"
     green = "green.dat"

このファイルはTOML形式で記述され、内容ごとにセクションに分類されています。

``[log]`` セクション
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ログ出力に関する設定を行います。
``print_level`` で標準出力のレベル、``print_step`` でログ出力を行う繰り返し間隔を指定します。

``[mode]`` セクション
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

実行モードに関する設定および基本パラメータの指定を行います。
``mode`` で実空間版(``UHF``)または波数空間版(``UHFk``)を選択します。
``[mode.param]`` サブセクションには計算用のパラメータを指定します。

``[file]`` セクション
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``[file.input]`` サブセクションでは、入力ファイルを格納するディレクトリ ``path_to_input`` および
初期配位データファイルのファイル名 ``initial`` を指定します。指定がない場合は乱数を用いて初期化されます。
``[file.input.interaction]`` サブセクションには、幾何情報および相互作用定義を格納するファイルのファイル名を相互作用のタイプごとに列挙します。

``[file.output]`` サブセクションには、エネルギーなどの物理量を出力するファイル名 ``energy`` 、
ハミルトニアンの固有値・固有ベクトルを出力するファイル名 ``eigen`` 、一体グリーン関数を書き出す出力ファイル名 ``green`` を指定します。
これらのキーワードがない場合にはその項目は出力されません。

詳細についてはファイルフォーマットの章をご覧ください。


相互作用定義ファイルの作成
----------------------------------------

Hamiltonianを構築するための格子の幾何情報および相互作用係数を格納したデータファイルを作成します。
項目とファイル名の対応付けは、入力パラメータファイルの ``[file.input.interaction]`` セクションで行います。

``Geometry``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

格子の幾何情報を記述します。ファイル例を以下に示します。

::

   3.7599302871   0.0000000000   0.0000000000
   0.0000000000   3.7599302871   0.0000000000
   0.0000000000   0.0000000000   5.4822004186
        10
   -7.179835091886330E-003 -3.812050198019962E-002  1.639284152926924E-003
    1.346463812592166E-002  6.709778405878775E-003 -6.812442303544219E-003
    0.495705070884200      -0.457955704941170      -4.077818544354700E-003
   -1.577970702078565E-004 -2.999005205319096E-004 -1.190284144276225E-004
   -1.302397074478660E-003 -5.021621895411691E-003 -3.514564279609852E-004
    0.504124376959700       0.457760356450585      -2.634809811615298E-003
    0.499384075989520      -0.494227365093439       6.927730957590197E-003
   -5.164444920392309E-003  3.667887236852975E-002  4.972296517752579E-003
    0.500170586121734       0.499747448247510       2.760670734661295E-003
    0.500734036298328       0.494793997305026      -2.212377045150314E-003

基本ベクトル(1〜3行目)、軌道の数(4行目)、各軌道のWannier center(5行目以降)を記載します。

``Transfer``, ``Coulomb``, ``Hund``, etc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Transferに指定するファイルは、電子系のTransferに相当するHamiltonianの係数を格納します。
また、二体相互作用の係数は相互作用のタイプごとに係数を格納するファイルを指定します。

相互作用のタイプは、実空間版UHFの入力であるExpertModeに合わせて、
Coulomb, Hund, Ising, Exchange, PairLift, PairHop が定義されています。
RESPACK との接続を考慮したセットを用意することも検討中です。

これらのファイルはWannier90形式で記述されます。以下に例を示します。
::

   wannier90 format for vmcdry.out or HPhi -sdry
       10
      245
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
   ...(略)
    1    1    1    1    1
   -3   -3   -2    1    1  -0.0000269645  -0.0000000000
   -3   -3   -2    1    2  -0.0000071722  -0.0000018600
   -3   -3   -2    1    3  -0.0000083990   0.0000010972
   -3   -3   -2    1    4  -0.0000000990   0.0000000427
   -3   -3   -2    1    5  -0.0000018628  -0.0000003609
   -3   -3   -2    1    6  -0.0000129504  -0.0000014047
   -3   -3   -2    1    7  -0.0000189169   0.0000024697
   -3   -3   -2    1    8   0.0000238115   0.0000014316
   -3   -3   -2    1    9   0.0000036708  -0.0000003266
   -3   -3   -2    1   10   0.0000361752   0.0000003247
   -3   -3   -2    2    1  -0.0000071722   0.0000018600
   -3   -3   -2    2    2   0.0000105028  -0.0000000000
   ...(略)

コメント行(1行目)、軌道の数(2行目)、並進ベクトルの数 ``nrpts`` (3行目)、
縮重度 ( ``nrpts`` 個を1行あたり15個ずつ)、係数行列の要素を記載します。
行列要素の各行は、並進ベクトル :math:`r_x, r_y, r_z`、軌道のインデックス :math:`\alpha, \beta`、
係数の値の実部・虚部です。

   
計算の実行
----------------------------------------

全ての入力ファイルが準備できた後、プログラムを実行して計算を行います。
入力パラメータファイル(ここでは ``input.toml`` )を引数とし、ターミナルからH-waveを実行します。

.. code-block:: bash

    $ python3 path_to_H-wave/qlms.py input.toml

計算が開始されると以下のようなログが出力されます。

::
   
   2022-08-02 10:16:57,041 INFO qlms: Read definitions from files
   2022-08-02 10:16:57,041 INFO qlms.read_input: >>> QMLSkInput init
   2022-08-02 10:16:57,041 INFO qlms.read_input: QMLSkInput: read Gemoetry from dir-model/zvo_geom.dat
   2022-08-02 10:16:57,041 INFO qlms.read_input: QMLSkInput: read interaction Transfer from dir-model/zvo_hr.dat
   2022-08-02 10:16:57,068 INFO qlms.read_input: QMLSkInput: read interaction Hund from dir-model/zvo_jr.dat
   2022-08-02 10:16:57,091 INFO qlms: Get Hamiltonian information
   2022-08-02 10:16:57,091 INFO qlms: Get output information
   2022-08-02 10:16:57,091 ERROR qlms.read_input: Get_param: key must be mod or ham or output.
   2022-08-02 10:16:57,091 INFO qlms.uhfk: Show parameters
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     Cell Shape  = (8, 8, 8)
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     Cell volume = 512
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     Num orbit   = 10
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     nspin       = 2
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     nd          = 20
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     Ncond       = 1280
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     T           = 1.0
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     Mix         = 0.5
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     RndSeed     = 123456789
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     IterationMax= 20
   2022-08-02 10:16:57,091 INFO qlms.uhfk:     EPS         = 1e-08
   2022-08-02 10:16:57,091 INFO qlms: Start UHF calculation
   2022-08-02 10:16:57,091 INFO qlms.uhfk: Start UHFk calculations
   2022-08-02 10:16:57,091 INFO qlms.uhfk: step, rest, energy, NCond, Sz
   2022-08-02 10:16:57,091 INFO qlms.uhfk: >>> _make_ham_trans
   2022-08-02 10:16:57,097 INFO qlms.uhfk: >>> _make_ham_inter
   2022-08-02 10:16:57,100 INFO qlms.uhfk: >>> _initial_green
   2022-08-02 10:16:57,100 INFO qlms.uhfk: initialize green function with random numbers
   2022-08-02 10:16:57,102 INFO qlms.uhfk: >>> _make_ham
   2022-08-02 10:16:57,102 INFO qlms.uhfk: Transfer
   2022-08-02 10:16:57,103 INFO qlms.uhfk: Hund
   2022-08-02 10:16:57,132 INFO qlms.uhfk: >>> _diag
   2022-08-02 10:16:57,211 INFO qlms.uhfk: >>> _green
   2022-08-02 10:16:57,578 INFO qlms.uhfk: mu = -11.44737114523863
   2022-08-02 10:16:57,592 INFO qlms.uhfk: >>> _calc_energy
   2022-08-02 10:16:57,602 INFO qlms.uhfk: energy: Band = (-16112.29460022662-1.1407667914167632e-15j)
   2022-08-02 10:16:57,606 INFO qlms.uhfk: energy: Hund = (-1.0089285782652733-3.4322714015332136e-19j)
   2022-08-02 10:16:57,606 INFO qlms.uhfk: >>> _calc_phys
   2022-08-02 10:16:57,606 INFO qlms.uhfk: ncond = (3.6639486880793397-8.920951000482291e-18j)
   2022-08-02 10:16:57,606 INFO qlms.uhfk: sz = (-1.2515495809273247e-05+2.6166965323053473e-17j)
   2022-08-02 10:16:57,613 INFO qlms.uhfk: rest = 1.50790573236065
   2022-08-02 10:16:57,617 INFO qlms.uhfk: 0, 7.362821e-06, -16113.304, 3.664, -6.258e-06 
   ...

入力ファイル読み込みに関するログが出力されたあと、波数空間UHF計算の計算過程に関する情報が出力されます。
出力ファイルは ``input.toml`` の ``[file.output]`` セクションの指定に従い、
``output`` ディレクトリに ``energy.dat`` , ``eigen.dat.npz``, ``green.dat.npz`` ファイルが出力されます。

出力ファイルの詳細についてはファイルフォーマットの章をご覧ください。
   
