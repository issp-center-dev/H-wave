.. highlight:: none

.. _subsec:eigen_uhfk.dat:

eigen
~~~~~~~~~~

収束したハミルトニアンの固有値、固有ベクトルをnpz形式で出力します。
ファイル名は、環境設定ファイルの ``file.output`` セクション ``eigen`` で指定された文字列 (以下、 ``eigen_str`` ) を用いて、 ``eigen_str.npz`` という名前で出力されます。

以下、データを読み込む例となります。

.. code-block:: python

    import numpy as np
    data = np.load("eigen_str.npz")
    eigenvalue = data["eigenvalue"]
    eigenvector = data["eigenvector"]
    wavevec = data["wavevec"]

``eigenvalue`` には波数ごとの固有値 :math:`\lambda_l(\vec{k})` が格納されます。
副格子を指定している場合は、副格子を単位とした値になります。
データ形式は numpy ndarray で、データの並びは ``eigenvalue[k][l]`` です。
``k`` は波数ベクトル :math:`\vec{k}` を一次元化したインデックス、
``l`` はセル内の固有値のインデックスです。
Sz固定の場合は、固有値のインデックスは軌道 ``a`` とスピン ``s`` (up-spin は 0, down-spin は 1)
に対して ``a + Norb * s`` となります。 ``Norb`` はセル内の軌道数です。

``eigenvector`` には対応する固有ベクトルが格納されます。
データ形式は numpy ndarray で、データの並びは ``eigenvector[k][l][j]`` です。
``k``, ``l`` は対応する波数および固有値のインデックス、
``j`` はセル内の軌道・スピンのインデックスです。

ファイルには波数の情報も出力されます。
``wavevec`` には波数と固有値のインデックスに対応する波数ベクトルの大きさが
格納されます。データ形式は numpy ndarray で、データの並びは ``wavevec[k][l]`` です。
``k`` は波数ベクトル :math:`\vec{k}` を一次元化したインデックス、
``l`` はセル内の固有値のインデックスです。

.. raw:: latex
