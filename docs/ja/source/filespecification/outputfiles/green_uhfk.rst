.. highlight:: none

.. _Subsec:green_uhfk:

green (UHFk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

一体グリーン関数\ :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle`\ の計算結果をnpz形式で出力します。
ファイル名は環境設定ファイルの中の ``file.output`` セクションでキーワード ``green`` を用いて指定されます。

データはキー ``green`` にバインドされます。
データ配列は ``ndarray(r, s, a, t, b)`` で、インデックスは以下のとおりです。

-  ``r``: 並進ベクトル :math:`[r_x\ r_y\ r_z]` を1次元化したインデックス
-  ``a``, ``b``: 軌道のインデックス :math:`\alpha, \beta`
-  ``s``, ``t``: スピンのインデックス :math:`\sigma, \sigma^\prime`

出力ファイルは、``file.input`` セクションの ``initial`` で指定するグリーン関数の初期データとして使用できます。
   
.. raw:: latex

   \newpage
