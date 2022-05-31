#####################################################
H-waveのドキュメントへようこそ！
#####################################################

H-waveとは？
------------------------------------------
H-waveでは遍歴電子系を対象に非制限Hartree-Fock(UHF)近似を行うためのプログラムです。
UHF近似では揺らぎ :math:`\delta A \equiv A-\langle A \rangle` の一次までを考慮することで、二体項を一体項へと近似します。
たとえば、サイト間クーロン相互作用

.. math::

   {\cal H}_V = \sum_{i,j, \sigma, \sigma'}V_{ij} n_ {i\sigma}n_{j\sigma'}

について考えます。簡単化のため、 :math:`i\equiv (i, \sigma)`,
:math:`j\equiv (j, \sigma')` とすると相互作用の項は揺らぎの二次を落とすことで、

.. math::

   \begin{aligned}
   n_ {i}n_{j} &=
   (\langle n_{i} \rangle +\delta n_i) (\langle n_{j} \rangle +\delta n_j)
   - \left[ \langle c_{i}^{\dagger}c_j \rangle +\delta (c_{i}^{\dagger}c_j ) \right]
     \left[ \langle c_{j}^{\dagger}c_i \rangle +\delta (c_{j}^{\dagger}c_i )\right]
   \nonumber\\
   &\sim
   \langle n_{i} \rangle n_j+\langle n_{j} \rangle  n_i
   - \langle c_{i}^{\dagger}c_j \rangle  c_{j}^{\dagger}c_i  -  \langle c_{j}^{\dagger}c_i \rangle c_{i}^{\dagger}c_j 
   -\langle n_{i} \rangle \langle n_j \rangle +  \langle c_{j}^{\dagger}c_i \rangle \langle c_{i}^{\dagger}c_j \rangle
   \end{aligned}

と近似されます。このような形式で、その他の相互作用についても近似を行うことで、一体問題に帰着させることができます。
計算では、上記の各平均値がself-consistentになるまで計算を行います。

ライセンス
--------------
本ソフトウェアのプログラムパッケージおよびソースコード一式はGNU General Public License version 3（GPL v3）に準じて配布されています。



コピーライト
------------------

© *2022- The University of Tokyo. All rights reserved.*

本ソフトウェアは2022年度 東京大学物性研究所 ソフトウェア高度化プロジェクトの支援を受け開発されており、その著作権は東京大学が所持しています。

ダウンロード
------------------


Contents
--------
.. toctree::
   :maxdepth: 3
   :numbered: 3

   introduction_ja
   howtouse/ho-index
   tutorial/tu-index
   filespecification/fi-index
   algorithm/al-index
   acknowledgement_ja
   

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
