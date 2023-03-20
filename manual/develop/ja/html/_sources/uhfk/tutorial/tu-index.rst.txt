==================
チュートリアル
==================

H-waveを波数空間モード(UHFk)で実行するには、入力ファイルとして

1. 環境設定入力ファイル
2. 相互作用定義ファイル

を用意した後、プログラムを実行します。2.は RESPACK 等の外部プログラムの出力を利用する他、StdFaceライブラリを使って生成することもできます。

以下では、 ``docs/tutorial/Hubbard/UHFk`` ディレクトリにあるサンプルを例にチュートリアルを実施します。
相互作用定義ファイルは StdFace ライブラリを用いて生成します。

環境設定入力ファイルの作成
--------------------------------

環境設定入力ファイルには、基本パラメータの指定と入出力を制御する情報を記述します。
``docs/tutorial/Hubbard/UHFk`` ディレクトリ内に ``input.toml`` というファイルがありますが、これが入力パラメータファイルになります。
以下、ファイルの内容を記載します。

.. literalinclude:: ../sample/input.toml

このファイルはTOML形式で記述され、内容ごとにセクションに分類されています。

``[log]`` セクション
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ログ出力に関する設定を行います。
``print_level`` で標準出力のレベル、 ``print_step`` でログ出力を行う繰り返し間隔を指定します。

``[mode]`` セクション
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

実行モードに関する設定および基本パラメータの指定を行います。
``mode`` で実空間版(``UHF``)または波数空間版(``UHFk``)を選択します。
``[mode.param]`` サブセクションには計算実行時のパラメータを指定します。

``[file]`` セクション
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``[file.input]`` サブセクションでは、入力ファイルを格納するディレクトリ ``path_to_input`` および
初期配位データファイルのファイル名 ``initial`` を指定します。指定がない場合は乱数を用いて初期化されます。
``[file.input.interaction]`` サブセクションには、幾何情報および相互作用定義を格納するファイルのファイル名を相互作用のタイプごとに列挙します。

``[file.output]`` サブセクションには、エネルギーなどの物理量を出力するファイル名 ``energy`` 、
ハミルトニアンの固有値・固有ベクトルを出力するファイル名 ``eigen`` 、一体グリーン関数を書き出す出力ファイル名 ``green`` を指定します。
これらのキーワードがない場合にはその項目は出力されません。

詳細については :ref:`ファイルフォーマット<Ch:Config>` の章をご覧ください。


相互作用定義ファイルの作成
----------------------------------------

Hamiltonianを構築するための格子の幾何情報および相互作用係数を格納したデータファイルを作成します。
項目とファイル名の対応付けは、入力パラメータファイルの ``[file.input.interaction]`` セクションで行います。

``Geometry``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

格子の幾何情報を記述します。ファイル例を以下に示します。

.. literalinclude:: ../sample/geom.dat

基本ベクトル(1〜3行目)、軌道の数(4行目)、各軌道のWannier center(5行目以降)を記載します。

``Transfer``, ``CoulombIntra``, ``CoulombInter``, ``Hund``, etc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Transferに指定するファイルは、電子系のTransferに相当するHamiltonianの係数を格納します。
また、二体相互作用の係数は相互作用のタイプごとに係数を格納するファイルを指定します。

相互作用のタイプは、実空間版UHFの入力ファイル形式と対応して、CoulombItra, CoulombInter, Hund, Ising, Exchange, PairLift, PairHop が定義されています。

これらのファイルはWannier90(-like)形式で記述されます。以下に例を示します。

.. literalinclude:: ../sample/transfer.dat

コメント行(1行目)、軌道の数(2行目)、並進ベクトルをすべて収める直方体内のセルの総数 ``nrpts`` (3行目)、
縮重度 ( ``nrpts`` 個を1行あたり15個ずつ)、係数行列の要素を記載します。

行列要素の各行は、並進ベクトル :math:`r_x, r_y, r_z` 、軌道のインデックス :math:`\alpha, \beta` 、係数の値の実部・虚部です。

   
計算の実行
----------------------------------------

全ての入力ファイルが準備できた後、プログラムを実行して計算を行います。
入力パラメータファイル(ここでは ``input.toml`` )を引数とし、ターミナルからH-waveを実行します。

.. code-block:: bash

    $ hwave input.toml

計算が開始されると以下のようなログが出力されます。

.. literalinclude:: ../sample/run.log

入力ファイル読み込みに関するログが出力されたあと、波数空間UHF計算の計算過程に関する情報が出力されます。
出力ファイルは ``input.toml`` の ``[file.output]`` セクションの指定に従い、
``output`` ディレクトリに ``energy.dat`` , ``eigen.npz``, ``green.npz`` ファイルが出力されます。

出力ファイルの詳細については :ref:`ファイルフォーマット<Sec:outputfile_uhfk>` の章をご覧ください。


状態密度の計算 (``hwave_dos``)
----------------------------------------

ポストツール ``hwave_dos`` を用いることで、状態密度を計算することができます。
ブリルアンゾーン積分を精度良く計算するために `libtetrabz <https://pypi.org/project/libtetrabz/>`_ を利用しています。 ``pip`` を利用してインストールしてください。 ::

    $ python3 -m pip install libtetrabz

``hwave_dos`` は、 ``hwave`` で利用した入力パラメータファイルを引数として受け取ります ::

    $ hwave_dos input.toml

``hwave_dos`` は、 ``hwave`` と同様に ``[file.output]`` セクションの指定に従い、
``output`` ディレクトリに ``dos.dat`` ファイルを出力します。
ファイル名は ``--output`` オプションで変更することができます。 ::

    $ hwave_dos input.toml --output dos.dat

状態密度を計算するエネルギーの範囲は ``--ene-window`` オプションで指定します。
省略した場合は、 ``hwave`` で得られたエネルギーの最小値と最大値を :math:`E_\text{min}`, :math:`E_\text{max}` として、 :math:`[E_\text{min}-0.2, E_\text{max}+0.2]` で計算されます。
エネルギーの点数は ``--ene-num`` オプションで指定します（デフォルトは101） ::

    $ hwave_dos input.toml --ene-window -10.0 5.0 --ene-num 201

``--plot`` オプションを指定すると、状態密度をプロットします。 ``matplotlib`` が必要です。 ::

    $ hwave_dos input.toml --plot dos.png


StdFaceライブラリのコンパイルと実行
----------------------------------------------------------------

相互作用定義ファイルは、StdFaceライブラリを使用することで簡単に作成することができます。
以下ではStdFaceライブラリを使用した入力ファイルの作成を行う方法を記載します。

波数空間UHFの入力ファイル形式に対応した StdFace ライブラリは以下から取得できます。

.. code-block:: bash

    $ git clone https://github.com/issp-center-dev/StdFace.git

ダウンロード終了後、以下のコマンドでコンパイルを行います。

.. code-block:: bash

    $ cd StdFace
    $ mkdir build && cd build
    $ cmake -DUHF=ON ..
    $ make

コンパイルに成功すると、 ``src`` ディレクトリに実行ファイル ``uhf_dry.out`` ができます。

``uhf_dry.out`` の入力ファイルはサンプルディレクトリ内の ``stan.in`` ファイルです。
ファイルの内容は以下のとおりです。

.. literalinclude:: ../sample/stan.in

- ``model`` は対象となる模型を指定するキーワードです。現状では電子数を固定したHubbard模型 ``Hubbard`` のみに対応しています。
- ``lattice`` は結晶構造を指定するキーワードです。 ここでは正方格子 ``square`` を選択しています。 ``W``, ``L`` は格子のサイズです。
- ``t`` はホッピング、 ``V`` は隣接サイトクーロン相互作用のパラメータです。
- ``calcmode = "uhfk"`` は波数空間UHF向けに Wannier90(-like)形式での出力を指定します。 ``exportall = 0`` は出力をコンパクトにするオプションです。

入力ファイルの詳細については、 :ref:`Ch:HowToWannier90` のセクションを参照してください。

上記のファイルを入力ファイルとして、 ``uhf_dry.out`` を以下のコマンドで実行します。

.. code-block:: bash

    $ cd path_to_Hwave/docs/tutorial/Hubbard/UHFk
    $ ln -s path_to_Stdface/build/src/uhf_dry.out .
    $ ./uhf_dry.out stan.in

実行終了後、実行ディレクトリに幾何情報ファイル ``geom.dat`` 、相互作用定義ファイル ``transfer.dat``, ``coulombinter.dat`` が生成されます。

