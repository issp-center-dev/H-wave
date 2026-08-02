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
    
    wavevector_unit = data["wavevector_unit"]
    wavevector_index = data["wavevector_index"]
    twist_offset = data["twist_offset"]


``eigenvalue`` には波数ごとの固有値 :math:`\lambda_l(\vec{k})` が格納されます。
副格子を指定している場合は、副格子を単位とした値になります。
データ形式は numpy ndarray で、データの並びは ``eigenvalue[k][l]`` です。
``k`` は波数ベクトル :math:`\vec{k}` を一次元化したインデックス、
``l`` はセル内の固有値のインデックスです。
Sz固定の場合は、固有値のインデックスは軌道部分 ``l'`` とスピン ``s`` (up-spin は 0, down-spin は 1)
に対して ``l' + Norb * s`` となります。 ``Norb`` はセル内の軌道数です。

``eigenvector`` には対応する固有ベクトルが格納されます。
データ形式は numpy ndarray で、データの並びは ``eigenvector[k][j][l]`` です。
``k``, ``l`` は対応する波数および固有値のインデックス、
``j`` はセル内の軌道・スピンのインデックスです。

ファイルには波数の情報も出力されます。
``wavevector_unit`` には逆格子ベクトル :math:`\vec{b}_i` を用いて
:math:`2\pi\vec{b}_i/N_i` で表される単位波数ベクトルが格納されます。
``wavevector_index`` には波数のインデックスと1次元化したインデックスとの対応が格納されます。
インデックス ``k`` に対応する波数ベクトルは以下で求められます。

.. code-block:: python

    k_vec = np.dot(wavevector_index[k], wavevector_unit)

``twist_offset`` には計算実行時の境界条件に対応する方向ごとのツイストオフセット
:math:`\theta_d / (2\pi)` が格納されます。各成分は ``periodic`` で ``0`` 、
``antiperiodic`` で ``1/2`` です。任意の境界条件下での物理波数ベクトルは
インデックス ``k`` に対して次のように求められます。

.. code-block:: python

    k_vec_phys = np.dot(wavevector_index[k] + twist_offset, wavevector_unit)

全方向周期境界の場合は ``twist_offset`` は全成分ゼロとなり、このフィールドを
無視する既存の読み込みコードは従来通りの整数グリッド解釈となります(後方互換)。


.. raw:: latex
