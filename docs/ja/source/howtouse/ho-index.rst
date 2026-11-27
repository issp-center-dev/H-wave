***********************************
基本的な使用方法
***********************************

- 必要なライブラリ・環境

    H-waveを利用するには、以下のプログラムとライブラリが必要です。

    - python 3.10 以降
    - numpy モジュール（1.26 以上 3 未満。CI では numpy 1.26 と 2.x を検証）
    - scipy モジュール
    - requests モジュール
    - tomli モジュール
    - cupy モジュール（任意。``gpu = true``\ によるGPU実行にのみ必要です
      — RPA/FLEX は\ ``[mode.param]``\ 、動的 Eliashberg ソルバーは
      ``[eliashberg]``\ で指定。CUDAバージョンに合ったビルド済みwheel
      （CUDA 12.x なら\ ``pip install cupy-cuda12x``\ ）の利用を推奨します。詳細は
      `CuPyインストールガイド <https://docs.cupy.dev/en/stable/install.html>`_
      を参照してください。CuPyが無い場合、``gpu = true``\ は警告を出して
      CPU実行にフォールバックします。）
    - sparse-ir モジュール（任意。動的 Eliashberg ソルバーの IR 基底松原軸
      ``[eliashberg] matsubara_basis = "ir"``\ および FLEX ソルバーの
      ``[mode.param] matsubara_basis = "ir"``\ にのみ必要です。
      sparse-ir 2.0 以降が必要です。\ ``pip install "hwave[ir]"``\ （または\
      ``pip install "sparse-ir>=2"``\ ）でインストールできます。未導入で\ ``"ir"``
      を指定した場合はインストール手順つきの明示エラーになります。）
    - threadpoolctl モジュール（任意。\ ``longitudinal_bond_channels``\ で
      ``fft_workers``\ を 2 以上にしたとき、ボンドチャネルのホスト側バッチ計算を
      スレッド並列化します。\ ``pip install "hwave[threads]"``\ （または\
      ``pip install threadpoolctl``\ ）でインストールできます。未導入の場合は
      警告を出して直列実行します。）

    なお、H-waveのUHFk、RPA、FLEXモードでは\ `numpy.fft <https://numpy.org/doc/stable/reference/generated/numpy.fft.fft.html>`_\ をFFT計算に利用しています。
    RPA・FLEX・動的Eliashbergソルバーの空間FFTは、オプトインの\ ``fft_workers``
    パラメータで並列FFTワーカーを指定した場合は\ `scipy.fft <https://docs.scipy.org/doc/scipy/reference/fft.html>`_\ に、
    GPU実行時は cuFFT に切り替わります。

- Python からの計算結果の読み込み

    ``hwave.results``\ パッケージを使うと、UHFr・UHFk の計算の出力ファイルを、
    その計算の入力 TOML ファイルを起点として Python から読み込めます。例:

    .. code-block:: python

       import hwave.results as hr

       run = hr.open_run("input.toml")
       status = hr.read_scf_status(run)   # 収束判定・残差・反復回数
       energy = hr.read_energy(run)

    各読み込み関数は出力ファイルに記録された規約の情報と計算設定を確認し、出力が存在しない、
    この計算で書かれたものではない、または未知の形式である場合には、推測せずに
    ヒントつきのエラーを送出します。関数と結果の型の説明は docstring
    （例: ``help(hr)``\ ）を参照してください。

- 初期状態の候補（プレビュー）

    ``python -m hwave.initial`` は UHFk の入力に対して初期状態の候補を用意して
    候補ごとに実行ディレクトリを作り、 ``hwave`` の実行後に結果を集めて順位を
    付けます（ ``prepare``, ``check``, ``collect``\ 。 ``python -m hwave.initial
    --help`` を参照）。 ``--tier small`` は、事前の常磁性 Hartree-Fock 計算から
    軌道のパターンと不安定性の行を加えます（ ``--para-check`` はその計算を一様な
    軌道占有からも始めます。 ``--stage0`` を付けるとその計算を別のジョブとして
    実行し、実行されるまで ``prepare`` は終了コード 3 で終わります）。
    完全な説明ができるまではプレビュー版です。詳しくは
    :doc:`initial_states` を参照してください。選ばれる状態は、試した候補の中で
    自由エネルギーが最も低いものであり、基底状態であることの証明ではありません。

- Official Page

    - `GitHubリポジトリ <https://github.com/issp-center-dev/H-wave>`_

    - `サンプル・チュートリアル <https://isspns-gitlab.issp.u-tokyo.ac.jp/hwave-dev/hwave-gallery>`_

- ダウンロード方法

    gitを利用できる場合は、以下のコマンドでH-waveをダウンロードすることができます。

    .. code-block:: bash

       $ git clone https://github.com/issp-center-dev/H-wave.git

- インストール方法

    - PyPI から

      H-wave は PyPI ソフトウェアリポジトリに登録されています。以下のコマンドで H-wave を
      インストールできます。

      .. code-block:: bash

        $ pip install hwave

    - ソースパッケージから

      H-wave のソースパッケージは配布サイトから取得できます。

      https://github.com/issp-center-dev/H-wave/releases

      また、git を用いて最新版を開発サイトからダウンロードできます。

      .. code-block:: bash

        $ git clone https://github.com/issp-center-dev/H-wave.git

      H-waveをダウンロード後、以下のコマンドを実行してインストールします。
      H-waveが利用するライブラリも必要に応じてインストールされます。

      .. code-block:: bash

        $ cd ./H-wave
        $ pip install .

- ディレクトリ構成

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
       
- 基本的な使用方法

  #. 入力ファイルの作成

     最初にH-wave用の入力ファイルを作成します。計算条件や入出力ファイル・ディレクトリなどの指定と、Hamiltonianの定義ファイルなどを作成する必要があります。
     後者は、`StdFaceライブラリ <https://github.com/issp-center-dev/StdFace>`_\ の利用が便利です。
     各ファイルの簡単な紹介はチュートリアルの章に記載されています。
     詳細についてはファイルフォーマットの章を参照してください。

  #. コマンドの実行

     入力ファイルのあるディレクトリで、以下のコマンドを実行することで、計算が行われます。

     .. code-block:: bash

        $ hwave input.toml

     または、

     .. code-block:: bash

        $ python3 path_to_H-wave/qlms.py input.toml

     計算終了後、計算結果が出力ディレクトリに出力されます。
     出力ファイルの詳細については、ファイルフォーマットの章を参照してください。

     また、本パッケージには後処理ツール\ ``hwave_dos``\ （状態密度）および
     ``hwave_sc``\ （Eliashberg方程式・超伝導解析）が含まれており、同じ入力ファイルを利用します。
     また、単独で動く変換ツール\ ``hwave_sigma_split total.npz out.npz ...``\ （位置引数で
     指定し、入力ファイルは使いません）は、FLEX の全自己エネルギーアーカイブを
     ``flex_hartree_fock = true``\ の計算が\ ``sigma_init``\ として受理する分割形式に変換します。

.. toctree::
   :hidden:

   initial_states
