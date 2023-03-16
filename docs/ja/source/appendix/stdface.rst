.. _Ch:StdFace:

StdFaceを用いた相互作用ファイルの作成
========================================

StdFaceライブラリのコンパイル
------------------------------

H-waveで使用する相互作用定義ファイルは、StdFaceライブラリを使用することで簡単に作成することができます。
以下ではStdFaceライブラリを使用した入力ファイルの作成を行う方法を記載します。

H-waveに対応した StdFace ライブラリは以下から取得できます。

.. code-block:: bash

    $ git clone https://github.com/issp-center-dev/StdFace.git

ダウンロード終了後、以下のコマンドでコンパイルを行います。

.. code-block:: bash

    $ cd StdFace
    $ mkdir build && cd build
    $ cmake -DHWAVE=ON ..
    $ make

コンパイルに成功すると、 ``src`` ディレクトリに実行ファイル ``hwave_dry.out`` ができます。

StdFaceライブラリの使用
------------------------------

``hwave_dry.out`` の入力ファイルはサンプルディレクトリ内の ``stan.in`` ファイルにそれぞれ置いてあります。
ファイルの内容は以下のとおりです。

.. literalinclude:: ../rpa/sample/stan.in

- ``model`` は対象となる模型を指定するキーワードです。現状では電子数を固定したHubbard模型 ``Hubbard`` のみに対応しています。
- ``lattice`` は結晶構造を指定するキーワードです。 ここでは正方格子 ``square`` を選択しています。 ``W``, ``L`` は格子のサイズです。
- ``t`` はホッピング、 ``V`` は隣接サイトクーロン相互作用のパラメータです。
- ``calcmode = "uhfk"``, ``calcmode = "rpa"`` では Wannier90(-like)形式での入力ファイルが出力されます。 ``exportall = 0`` は出力をコンパクトにするオプションです。
    ``calcmode = "uhfr"`` はH-waveの実空間UHFプログラムUHFr向けの入力ファイルを出力します。デフォルトは、 ``calcmode = "uhfk"`` となっています。

入力ファイルの詳細については、:ref:`Ch:HowToExpert` , :ref:`Ch:HowToWannier90` , :ref:`Ch:HowToWannier90_rpa` のセクションを参照してください。

上記のファイルを入力ファイルとして、 ``hwave_dry.out`` を以下のコマンドで実行します。

.. code-block:: bash

    $ cd path_to_Hwave/docs/tutorial/Hubbard/RPA
    $ ln -s path_to_Stdface/build/src/hwave_dry.out .
    $ ./hwave_dry.out stan.in

実行終了後、実行ディレクトリに幾何情報ファイル ``geom.dat`` 、相互作用定義ファイル ``transfer.dat``, ``coulombinter.dat`` が生成されます。
StdFaceライブラリで指定できる格子や相互作用の種類については、 HΦもしくはmVMCのマニュアルにあるスタンダードモードの入力ファイル形式を参照してください。

