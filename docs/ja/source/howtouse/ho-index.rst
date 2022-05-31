***********************************
Basic usage of H-wave
***********************************

- 要件

    H-waveを利用するには、以下の環境が必要です。

    - python 3.x
    - numpy モジュール
    - requests モジュール
    - tomli モジュール

- Official Page

    `GitHubリポジトリ <git@github.com:issp-center-dev/H-wave.git>`_

- ダウンロード方法

    gitを持っている場合には、以下のコマンドでH-waveをダウンロードすることができます。

    .. code-block:: bash

       $ git clone git@github.com:issp-center-dev/H-wave.git

- ディレクトリ構成

        T.B.A.

- 基本的な使用方法

    1. 入力ファイルの作成

        最初にH-wave用の入力ファイルを作成します。
        入力ファイルの格納ディレクトリ、出力ファイルの格納ディレクトリなどを指定する入力ファイルの作成と、
        計算条件、Hamiltonianなどを指定する入力ファイルを作成する必要があります。
        後者は、`StdFaceライブラリ <https://github.com/issp-center-dev/StdFace>`_ の利用が便利です。
        各ファイルの簡単な紹介はチュートリアルの章に記載されています。
        それ以上の入力ファイルの詳細については、ファイルフォーマットの章を参照してください。

    2.　コマンドの実行

        入力ファイルのあるディレクトリで、以下のコマンドを実行することで、計算が行われます。

        .. code-block:: bash

            $ python3 path_to_H-wave/qlms.py input.toml

        計算終了後、計算結果が出力ディレクトリに出力されます。
        出力ファイルの詳細については、ファイルフォーマットの章を参照してください。

