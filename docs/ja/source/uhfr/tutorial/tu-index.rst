==================
チュートリアル
==================

H-waveでは、入力ファイルとして

#. 環境設定入力ファイル
#. 入力ファイルリスト
#. Hamiltonian作成用ファイル
#. 出力結果指定用ファイル

を用意した後、計算を行います。2-3については、StdFaceライブラリを使用することで簡単に作成することができます。
以下ではStdFaceライブラリを使用した入力ファイルの作成を行った後、1から4までのファイル形式と計算手順を説明します。

StdFaceライブラリのコンパイルと実行
------------------------------------------

ここでは、StdFaceライブラリを用いた入力ファイルの作成について紹介します。
StdFaceライブラリでは、あらかじめ決められた格子模型やWannier90形式で記述された有効模型の情報をもとに入力ファイルを生成します。
厳密対角化ソルバー :math:`{\mathcal H}\Phi` や多変数変分モンテカルロソルバーmVMCでも入力ファイルの生成に用いられています。
ここでは具体例に沿ってUHF用の入力ファイル生成を説明します。

最初に、StdFaceライブラリをダウンロードします。gitをお持ちの方は、以下のコマンドでダウンロードできます。

.. code-block:: bash

    $ git clone https://github.com/issp-center-dev/StdFace.git

ダウンロード終了後、以下のコマンドを入力し、コンパイルします。

.. code-block:: bash

    $ cd StdFace
    $ mkdir build && cd build
    $ cmake -DUHF=ON ../
    $ make

コンパイルに成功すると、 ``src`` ディレクトリに実行ファイル ``uhf_dry.out`` ができます。

以下、 ``sample/Hubbard`` ディレクトリにあるサンプルを例にチュートリアルを実施します。
ディレクトリ内に ``stan.in`` というファイルがありますが、これが ``uhf_dry.out`` の入力ファイルになります。
ファイルの内容は以下のとおりです。

::

    model = "Hubbard"
    lattice = "square"
    a0W = 2
    a0L = 2
    a1W = -2
    a1L = 2
    t = 1.0
    U = 8.0
    ncond = 8
    2Sz = 0

``model`` は対象となる模型を指定するキーワードです。現状では電子数を固定したHubbard模型 ``Hubbard`` のみに対応しています。

``lattice`` は結晶構造を指定するキーワードです。 ここでは正方格子 ``square`` を選択しています。
``a0W, a0L`` はx軸を ``(a0W, a0L)`` ベクトルとして、``a1W, a1L`` はy軸を ``(a1W, a1L)`` ベクトルとして与えるパラメータです。

``t`` はホッピング、 ``U`` はオンサイトクーロン相互作用、 ``ncond`` は全電子数、
``2Sz`` は全 ``Sz`` の値の二倍を与えます。

入力ファイルの詳細については、例えば :math:`{\mathcal H}\Phi` のマニュアルなどを参照してください。

上記のファイルを入力ファイルとして、 ``uhf_dry.out`` を以下のコマンドで実行します。

.. code-block:: bash

    $ cd path_to_Hwave/sample/Hubbard
    $ ln -s path_to_Stdface/build/src/uhf_dry.out .
    $ ./uhf_dry.out stan.in

実行終了後、実行ディレクトリに入力ファイルリストファイル、基本パラメータ用ファイル、Hamiltonian作成用ファイルが生成されます。

環境設定入力ファイルの作成
------------------------------------------

環境設定入力ファイルでは、入出力を制御する情報を記載します。
ディレクトリ内に ``input.toml`` というファイルがありますが、これが環境設定入力ファイルになります。
以下、ファイルの中身を記載します。

::

    [log]
    print_level = 1
    print_step = 20
    [mode]
    mode = "UHF"
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
    namelist = "namelist.def"
    [file.output]
    path_to_output = "output"
    energy = "energy.dat"
    eigen = "eigen.dat"
    green = "green.dat"


このファイルはtoml形式で記述します。

``log`` セクションに ``print_level`` で標準出力のレベル、 ``print_step`` でログファイルに出力するステップ間隔を指定します。

``mode`` セクションに実行モードおよび基本パラメータを指定します。

``file.input`` セクションに入力ファイルが格納されているディレクトリ ``path_to_input`` 、入力ファイルリストファイルの名前  ``namelist``  を指定します。

``file.output`` セクションには出力ファイルを格納するディレクトリ ``path_to_output`` を指定します。
また、エネルギーの値を出力するファイル名 ``energy`` 、ハミルトニアンの固有値を出力するファイル名 ``eigen`` 、一体グリーン関数の出力ファイル名 ``green`` を指定します。これらのキーワードがない場合には情報は出力されません。

詳細についてはファイルフォーマットの章をご覧ください。

入力ファイルリストファイル
------------------------------------------

入力ファイルの種類と名前を指定するファイルnamelist.defには、下記の内容が記載されています。
入力ファイルリストファイルでは、行毎にKeywordで指定するデータの種類と、そのデータを格納するファイル名を記述します。
詳細はセクション :ref:`Subsec:InputFileList` をご覧ください。 ::

         ModPara  modpara.def
           Trans  trans.def
    CoulombIntra  coulombintra.def
        OneBodyG  greenone.def

基本パラメータの指定
--------------------------

基本パラメータは環境設定ファイルinput.tomlの ``mode.param`` セクションに指定します。
``ModPara`` にひも付けられるファイル(ここではmodpara.def)は使用しません。

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

  
なお、 ``CoulombIntra`` 以外にも、Hamiltonianを簡易的に記載するためのファイル形式に対応しています。
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

    $ ln -s path_to_Hwave/python/qlms.py .
    $ python3 qlms.py input.toml

計算が開始されると以下のようなログが出力されます。

::

    2022-05-26 16:27:17,584 INFO qlms :Read def files
    2022-05-26 16:27:17,585 INFO qlms :Get Parameters information
    {'modpara': {'CDataFileHead': ['zvo'], 'CParaFileHead': ['zqp'], '--------------------': [], 'Nsite': ['8'], '2Sz': ['0'], 'Ncond': ['8'], 'IterationMax': ['1000'], 'EPS': ['8'], 'Mix': ['0.5000000000'], 'RndSeed': ['123456789']}}
    2022-05-26 16:27:17,585 INFO qlms :Get Hamiltonian information
    2022-05-26 16:27:17,585 INFO qlms :Get Output information
    2022-05-26 16:27:17,585 INFO qlms :Start UHF calculation
    2022-05-26 16:27:17,586 INFO qlms.uhf :Set Initial Green's functions
    2022-05-26 16:27:17,586 INFO qlms.uhf :Start UHF calculations
    2022-05-26 16:27:17,586 INFO qlms.uhf :step, rest, energy, NCond, Sz
    2022-05-26 16:27:17,589 INFO qlms.uhf :0, 0.022120078, -36.085839+0j, 8, 0.5424
    2022-05-26 16:27:17,628 INFO qlms.uhf :20, 0.0005230403, -5.6054878+0j, 8, 0.2641
    ...
		
入力ファイル読み込みに関するログが出力されたあと、UHF計算の計算過程に関する情報が出力されます。
出力ファイルは ``input.toml`` の ``file.output`` セクションでの設定にしたがい、
``output`` ディレクトリに ``energy.dat`` , ``eigen.dat``, ``green.dat`` ファイルが出力されます。
出力ファイルの詳細についてはファイルフォーマットの章をご覧ください。

