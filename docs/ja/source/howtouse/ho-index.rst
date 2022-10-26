***********************************
Basic usage of H-wave
***********************************

- 必要なライブラリ・環境

    H-waveを利用するには、以下のプログラムとライブラリが必要です。

    - python 3.x
    - numpy モジュール
    - scipy モジュール
    - requests モジュール
    - tomli モジュール

- Official Page

    `GitHubリポジトリ <git@github.com:issp-center-dev/H-wave.git>`_

- ダウンロード方法

    gitを利用できる場合は、以下のコマンドでH-waveをダウンロードすることができます。

    .. code-block:: bash

       $ git clone git@github.com:issp-center-dev/H-wave.git

- ディレクトリ構成

    ::

      .
      |-- COPYING
      |-- README.md
      |-- src/
      |   |-- __init__.py
      |   |-- qlms.py
      |   |-- solver/
      |   |   |-- __init__.py
      |   |   |-- base.py
      |   |   |-- uhf.py
      |   |   |-- uhfk.py
      |   |-- qlmsio/
      |   |   |-- __init__.py
      |   |   |-- read_input.py
      |   |   |-- read_input_k.py
      |   |   |-- wan90.py
      |-- sample/
      |   |-- Hubbard/
      |   |   |-- input.toml
      |   |   |-- std.in
      |   |-- UHFk/
      |   |   |-- input.toml
      |   |   |-- std.in
      |   |-- dir-model/

       
- 基本的な使用方法

  #. 入力ファイルの作成

     最初にH-wave用の入力ファイルを作成します。
     入力ファイルの格納ディレクトリ、出力ファイルの格納ディレクトリなどを指定する入力ファイルの作成と、
     計算条件、Hamiltonianなどを指定する入力ファイルを作成する必要があります。
     後者は、`StdFaceライブラリ <https://github.com/issp-center-dev/StdFace>`_ の利用が便利です。
     各ファイルの簡単な紹介はチュートリアルの章に記載されています。
     詳細についてはファイルフォーマットの章を参照してください。

  #. コマンドの実行

     入力ファイルのあるディレクトリで、以下のコマンドを実行することで、計算が行われます。

     .. code-block:: bash

            $ python3 path_to_H-wave/qlms.py input.toml

     計算終了後、計算結果が出力ディレクトリに出力されます。
     出力ファイルの詳細については、ファイルフォーマットの章を参照してください。

