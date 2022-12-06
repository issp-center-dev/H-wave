.. highlight:: none

.. _subsec:energy_uhfk.dat:

energy
~~~~~~~~~~

波数空間版UHF法で求めたエネルギー、粒子数、スピンに関する計算結果を出力します。
ファイル名は環境設定ファイルの中の ``file.output`` セクションでキーワード ``energy`` を用いて指定することができます。
以下にファイル例を記載します。

::

    Energy_Total = -5.88984624257707
    Energy_Band = -0.9265413257740396
    Energy_Coulomb = -4.963304916803031
    NCond = 8.000000000000007
    Sz = 3.2822430107160017e-07

ファイル形式
^^^^^^^^^^^^

-  Energy_Total = :math:`[`\ energy_total\ :math:`]`

-  Energy_Band = :math:`[`\ energy_band\ :math:`]`

-  Energy_{type} = :math:`[`\ energy_type\ :math:`]`

-  NCond = :math:`[`\ ncond\ :math:`]`

-  Sz = :math:`[`\ sz\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ energy_total\ :math:`]`

   **形式 :** float型

   **説明 :**
   UHF法で求めた固有ベクトルを用い計算した全エネルギー。

-  :math:`[`\ energy_band\ :math:`]`

   **形式 :** float型

   **説明 :** UHF法で求めたハミルトニアン行列の固有値のみ考慮した場合のエネルギー。


-  :math:`[`\ energy_type\ :math:`]`

   **形式 :** float型

   **説明 :** 相互作用分のエネルギー。相互作用のタイプごとに出力される。

-  :math:`[`\ ncond\ :math:`]`

   **形式 :** float型

   **説明 :** 全粒子数の期待値。
    :math:`\sum_{i}\langle n_{i}\rangle`

-  :math:`[`\ sz\ :math:`]`

   **形式 :** float型

   **説明 :** 全スピンの :math:`z` 成分 :math:`S_z` の期待値。
    :math:`\sum_{i}\langle (n_{i\uparrow}-n_{i\downarrow})\rangle/2`


.. raw:: latex
