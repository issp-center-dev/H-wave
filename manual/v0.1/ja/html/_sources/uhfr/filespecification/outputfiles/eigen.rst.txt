.. highlight:: none

.. _subsec:eigen.dat:

eigen
~~~~~~~~~~

収束したハミルトニアンの固有値、固有ベクトルをnpz形式で出力します。
ファイルは環境設定ファイルの ``file.output`` セクション ``eigen`` で指定された文字列 (以下、 ``eigen_str`` ) を用いて、
``{key}_eigen_str.npz`` という名前で出力されます。
ここで、 ``{key}`` は

- ``mode.param`` セクションで ``Sz`` を指定しない場合: ``sz-free``
- ``mode.param`` セクションで ``Sz`` を指定した場合: ``spin-up`` , ``spin-down``

となります ( ``Sz`` を指定した場合には2つのファイルが書き出されます)。
以下、データを読み込む例となります。

.. code-block:: python

    import numpy as np
    data = np.load("key_eigen_str.npz")
    eigenvalue = data["eigenvalue"]
    eigenvector = data["eigenvector"]

``eigenvalue`` には固有値が低い順に格納されます。 ``N`` を全サイト数とした場合、

- ``mode.param`` セクションで ``Sz`` を指定しない場合: 2 ``N`` 個
- ``mode.param`` セクションで ``Sz`` を指定した場合: ``N`` 個

の固有値が出力されます。

``eigenvector`` には対応する固有ベクトルが格納されます。
第1列目のindexは  ``n_site`` をサイトのindex、 ``n_spin`` をスピンのindex(upの場合に0, downの場合に1)として、

- ``mode.param`` セクションで ``Sz`` を指定しない場合: ``n_site`` + ``n_spin`` * ``N``
- ``mode.param`` セクションで ``Sz`` を指定した場合: ``n_site``

に対応し、第2列目のindexは固有値のindexに対応します。

.. raw:: latex

   \newpage


