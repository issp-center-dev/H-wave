.. highlight:: none

.. _Subsec:green_uhfk:

green
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

一体グリーン関数\ :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom{\dagger}}\rangle`\ の計算結果をnpz形式で出力します。
ファイル名は、環境設定ファイルの中の ``file.output`` セクションでキーワード ``green`` を用いて指定された文字列(以下、 ``green_str``)を用いて、 ``green_str.npz`` という名前で出力されます。

データはキー ``green`` にバインドされます。
データ配列は ``ndarray(r, s, a, t, b)`` で、インデックスは以下のとおりです。

-  ``r``: 並進ベクトル :math:`[r_x\ r_y\ r_z]` を1次元化したインデックス
-  ``a``, ``b``: 軌道のインデックス :math:`\alpha, \beta`
-  ``s``, ``t``: スピンのインデックス :math:`\sigma_1, \sigma_2`

出力ファイルは、``file.input`` セクションの ``initial`` で指定するグリーン関数の初期データとして使用できます。

副格子を指定している場合は、上記に加えて、副格子を単位としたグリーン関数の値がキー ``green_sublattice`` にバインドされます。並進ベクトルおよび軌道のインデックスは副格子に読み替えます。
   
以下、データを読み込む例となります。

.. code-block:: python

    import numpy as np
    data = np.load("green.dat.npz")
    green = data["green"]

.. raw:: latex
