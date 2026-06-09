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
    FreeEnergy_Total = -5.91234017896271
    FreeEnergy_Band = -0.9490352621596976
    NCond = 8.000000000000007
    Sz = 3.2822430107160017e-07
    ChemicalPotential = 0.123456789

ファイル形式
^^^^^^^^^^^^

-  Energy_Total = ``[energy_total]``

-  Energy_Band = ``[energy_band]``

-  Energy_{type} = ``[energy_type]``

-  FreeEnergy_Total = ``[free_energy_total]``

-  FreeEnergy_Band = ``[free_energy_band]``

-  NCond = ``[ncond]``

-  Sz = ``[sz]``

-  ChemicalPotential = ``[mu]``

パラメータ
^^^^^^^^^^

-  ``[energy_total]``

   **形式 :** float型

   **説明 :**
   UHF法で求めた固有ベクトルを用い計算した全エネルギー。

-  ``[energy_band]``

   **形式 :** float型

   **説明 :** UHF法で求めたハミルトニアン行列の固有値のみ考慮した場合のエネルギー。

-  ``[energy_type]``

   **形式 :** float型

   **説明 :** 相互作用分のエネルギー。相互作用のタイプごとに出力される。

-  ``[free_energy_total]``

   **形式 :** float型

   **説明 :**
   自由エネルギー :math:`F = E - TS` （Helmholtz自由エネルギー）。
   有限温度の自己無撞着計算で最小化される量。
   :math:`T=0` では内部エネルギー ``[energy_total]`` と一致する。

-  ``[free_energy_band]``

   **形式 :** float型

   **説明 :**
   自由エネルギーのバンド寄与 :math:`\mu N + \Omega_0`
   （:math:`\Omega_0 = -T\sum_n \ln(1+e^{-(\varepsilon_n-\mu)/T})`）。
   :math:`T=0` では ``[energy_band]`` と一致する。

-  ``[ncond]``

   **形式 :** float型

   **説明 :** 全粒子数の期待値。
    :math:`\sum_{i}\langle n_{i}\rangle`

-  ``[sz]``

   **形式 :** float型

   **説明 :** 全スピンの :math:`z` 成分 :math:`S_z` の期待値。
    :math:`\sum_{i}\langle (n_{i\uparrow}-n_{i\downarrow})\rangle/2`

-  ``[mu]``

   **形式 :** float型

   **説明 :**
   化学ポテンシャル（Fermi準位）。
   :math:`T=0` ではHOMO-LUMOギャップの中点（有限温度の :math:`\mu` の
   :math:`T\to 0^+` 極限）を採用する。
   なお絶縁体ではギャップ内で粒子数が :math:`\mu` に依らず一定のため、
   有限温度の :math:`\mu` はギャップ内で一意に定まらない点に注意。
   :math:`\mu` グループが複数ある場合は ``ChemicalPotential_{g}`` として
   グループごとに出力される。


.. raw:: latex
