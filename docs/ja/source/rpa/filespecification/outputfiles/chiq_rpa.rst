.. highlight:: none

.. _Subsec:chiq_rpa:

chiq, chi0q
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

感受率行列および既約感受率行列の計算結果を Numpy zip 形式で出力します。
ファイル名は、環境設定ファイルの中の ``file.output`` セクションでキーワード ``chiq`` および ``chi0q`` を用いて指定された文字列(以下 ``chiq_str`` と記述)を用いて、 ``chiq_str.npz`` という名前で出力されます。

データはそれぞれキー ``chiq`` および ``chi0q`` にバインドされます。
データ配列は ``ndarray(l,q,a,ap,b,bp)`` で、インデックスは以下のとおりです。

- ``l`` : 松原振動数のラベル。インデックスとラベルの対応はキー ``imat`` で与えられます。
- ``q`` : 波数ベクトルのインデックス :math:`[q_x\ q_y\ q_z]` を1次元化したインデックスで、
  :math:`q = q_z + N_z \cdot (q_y + N_y \cdot q_x)` となります。
- ``a``, ``ap``, ``b``, ``bp`` : 一般化軌道のインデックスで、それぞれ :math:`\tilde\alpha, \tilde\alpha^\prime, \tilde\beta, \tilde\beta^\prime` に対応します。
  一般化軌道インデックス :math:`\tilde\alpha` は軌道 :math:`\alpha`, スピン :math:`\sigma` に対して :math:`\tilde\alpha = \alpha + N_\text{orb}\cdot\sigma` となります。:math:`N_\text{orb}` は軌道数です。

出力する松原振動数の値または範囲をパラメータ ``matsubara_frequency`` で指定することができます。
  
``chi0q`` の出力ファイルは、計算済み既約感受率データとして ``file.input`` セクションの ``chi0q_init`` に指定して使用できます。

副格子を指定している場合は、出力されるデータは副格子を単位とした感受率の値です。波数ベクトルおよび軌道のインデックスは副格子に読み替えます。
   
以下、データを読み込む例となります。

.. code-block:: python

    import numpy as np
    data = np.load('output/chiq.npz')

    chiq = data['chiq']
    imat = data['imat']    

.. raw:: latex
