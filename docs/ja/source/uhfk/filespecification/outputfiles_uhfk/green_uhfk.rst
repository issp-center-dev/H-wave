.. highlight:: none

.. _Subsec:green_uhfk:

green
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

一体グリーン関数\ :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom{\dagger}}\rangle`\ の計算結果をnpz形式で出力します。
ファイル名は、環境設定ファイルの中の ``file.output`` セクションでキーワード ``green`` を用いて指定された文字列(以下、 ``green_str``)を用いて、 ``green_str.npz`` という名前で出力されます。

データはキー ``green`` にバインドされます。
データ配列は ``ndarray(r, s, a, t, b)`` で、インデックスは以下のとおりです。

-  ``r``: 並進ベクトル :math:`[r_x\ r_y\ r_z]` を1次元化したインデックスで、``r`` :math:`= r_z + N_z \cdot (r_y + N_y r_x)`
-  ``a``, ``b``: 軌道のインデックス :math:`\alpha, \beta`
-  ``s``, ``t``: スピンのインデックス :math:`\sigma_1, \sigma_2`

出力ファイルは、``file.input`` セクションの ``initial`` で指定するグリーン関数の初期データとして使用できます。

副格子を指定している場合は、上記に加えて、副格子を単位としたグリーン関数の値がキー ``green_sublattice`` にバインドされます。並進ベクトルおよび軌道のインデックスは副格子に読み替えます。

このとき、キー ``green_convention`` に値 ``green_slot_first`` がバインドされ、``green`` キーの折りたたみ符号の規約を記録します。これにより、(例えば RPA の ``green_init`` として再利用する際に)グリーン関数を曖昧さなく折りたたみ戻すことができます。

``BoundaryCondition`` に ``antiperiodic`` を含む方向がある場合、``green`` (および ``green_sublattice``) に保存されるグリーン関数は SCF が内部的に用いるゲージ変換後の基底でのものであり、反周期境界上の物理的な :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom{\dagger}}\rangle` ではない点に注意してください。適用された twist は同じ npz ファイル内のキー ``boundary_theta`` に記録されます。反周期境界での物理的なグリーン関数は ``onebodyg`` (``greenone.dat``) として出力されます(出力時に (i, j) ペアごとに逆ゲージ位相を掛けて戻しています)。反周期境界条件下で保存したグリーン関数を初期データとして再開する場合は、同じ ``BoundaryCondition`` で実行する必要があります(不一致は読み込み時に警告として出力されます)。

以下、データを読み込む例となります。

.. code-block:: python

    import numpy as np
    data = np.load("green.dat.npz")
    green = data["green"]

.. raw:: latex
