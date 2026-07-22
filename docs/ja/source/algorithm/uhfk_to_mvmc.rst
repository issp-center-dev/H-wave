.. highlight:: none

UHFk → mVMC PairProduct 変換
********************************

``tools/uhfk_to_mvmc.py`` は H-wave UHFk の自己無撞着解を、mVMC の
PairProduct 波動関数の初期値に変換する。本章では変換の規約と構成法を
述べる。利用方法は :doc:`../uhfk/tools/uhfk_to_mvmc` を参照。

変換は 4 つの段階からなる。占有 Slater 列の取り出し、副格子折り畳みの
ゲージ合成、pair 振幅 :math:`F` の構成、そして SOC を含む場合の
``trans.def`` 出力である。

Bloch 振幅の符号規約
====================

ブリッジは **正号 Bloch 振幅規約** を採る。すなわち実空間振幅を

.. math::

   c_R = \frac{1}{\sqrt{N_\mathrm{folded}}} \sum_k c_k e^{+i k R}

で定義する。この符号は H-wave の内部規約から決まるものであって、
選択の余地があるものではない。H-wave は反周期境界条件を負号ゲージ

.. math::

   \tilde{c}_r = e^{-i \theta r / L_\mathrm{phys}} c_r

で扱い、逆 Fourier 変換に ``ifftn(..., norm="forward")`` を用いる。
この 2 つを合成すると、tilde 側の Bloch 変換は正号となる。

``SubShape = [1, 1, 1]`` では :math:`(k, -k)` の時間反転対称和により
平面波因子が対称化されるため、密度は符号の選び方に依存しない。
``SubShape > [1, 1, 1]`` では折り畳み固有ベクトルの副格子 envelope が
非自明になるため両者は区別され、H-wave の ``greenone.dat`` と
element-wise に一致するのは正号規約のみである。

副格子折り畳みとゲージ合成
==========================

``SubShape`` が各方向で ``CellShape`` を割り切るとき、格子は超格子に
折り畳まれる。折り畳み後の軌道索引は、超格子内のセル位置
``sub_offset`` と元の軌道索引から符号化される。物理サイトはこの索引を
分解して再構成する。

折り畳み下で Slater 行列に乗る位相は 2 つあり、これらを合成する。

.. math::

   \exp\left[-i k_\mathrm{folded} \cdot (\mathrm{folded\_cell} + \mathrm{sub\_offset})\right]
   \times
   \exp\left[-i \theta \cdot r_\mathrm{phys} / L_\mathrm{phys}\right]

第 1 因子が副格子内オフセットのゲージ、第 2 因子が反周期境界条件の
twist である。``build_slater_orbitals`` は両者を合成した位相を
出力する Slater 行列に乗せる。

検証側の ``gauge_lift`` は同じ合成変換を用いて、H-wave の折り畳み
Bloch 基底の ``green_sublattice`` を物理基底へ持ち上げる。両者が同じ
変換を共有することにより、密度チェックは出力する Slater 行列そのものを
検査対象にできる。

pair 振幅 F の構成
==================

一体密度行列の規約は

.. math::

   G_{ij} = \langle c^\dagger_i c_j \rangle = \left(\overline{A} A^{T}\right)_{ij}

である。:math:`A` は占有 Slater 列を並べた行列を表す。

canonical pair
--------------

pair は :math:`\{k, \mathrm{partner}(k)\}` の非順序対ごとに、1 つの
canonical 行から発行する。canonical 行は次で定める。

- 自己対 (:math:`\mathrm{partner}(k) = k`) のとき、:math:`k` 自身を
  canonical とする。
- 非自己対のとき、波数索引の組が辞書式に小さい方を canonical とする
  (同順の場合は行索引の小さい方)。

:math:`F` は :math:`(2N_\mathrm{site}, 2N_\mathrm{site})` の複素反対称
行列であり、下三角は :math:`F_{ji} = -F_{ij}` で埋める。

partner balance
---------------

canonical 対の発行が矛盾なく閉じるためには、対をなす 2 行の占有数が
一致していなければならない。ブリッジは次のいずれかに該当する占有を
変換前に拒否する。

- :math:`n_\mathrm{occ}(k) \neq n_\mathrm{occ}(\mathrm{partner}(k))` である対が存在する
- canonical な自己対 :math:`k` の :math:`n_\mathrm{occ}(k)` が奇数である

この検査は ``validate_general_prerequisites`` が行い、条件を満たさない
占有はエラーとして報告する。黙って近似した出力を返すことはしない。

充填の選定にあたっては、スペクトルの gap だけでは不十分である。gap は
canonical 対の 2 行に占有がどう分配されるかについて何も述べないため、
partner balance の不変条件を別途確認する必要がある。

密度への射影
------------

:math:`F` から一体密度を得るには skew-SVD による射影を用いる。
特異値分解の左特異ベクトルのうち上位 :math:`2 N_\mathrm{pairs}` 本が
占有部分空間を張る。この射影は、出力した :math:`F` が意図した Slater
状態を表しているかを検証する際に用いる。

Transfer の写像規約
===================

スピン軌道結合を含む場合、mVMC の ``vmcdry.out`` が生成する
``trans.def`` はスピン対角成分のみを保持し、:math:`s \neq t` の項を
落とす。このためブリッジが H-wave の ``Transfer.dat`` から
``trans.def`` を直接出力する。

写像規約は次のとおりである。mVMC の規約は
:math:`H = -\sum \mathrm{trans}\, c^\dagger c` である。ComplexUHF は
裸のハミルトニアンを :math:`K = -\mathrm{trans}` として構成する。
H-wave の負号 Bloch 物理基底において、``Transfer.dat`` のエントリ
:math:`(R, s, t, v)` は変位 :math:`+R` で次のように写る。

.. math::

   K[i, t;\, i+R, s] &= \overline{v} \\
   \mathrm{trans}[i, t;\, i+R, s] &= -\overline{v}

サイト端点は :math:`i \to i+R` のまま変わらず、**スピン端点が入れ替わり**、
係数は共役かつ符号反転する。実数のスピン対角成分では入れ替えが no-op
となり、規則は :math:`\mathrm{trans} = -v` に帰着する。

この一般規則は導出されたものである。H-wave の
``ifftn(norm="forward")`` による :math:`e^{+ikR}` 変換、Hermitian
transpose、ComplexUHF の :math:`K = -\mathrm{trans}` という経路を
たどって得られる。数値による裏付けとして、保存された全固有対から
H-wave の裸の :math:`K` を再構成した比較では、本規則の最大絶対偏差は
:math:`1.1 \times 10^{-12}` であり、:math:`10^{-10}` を超えるエントリは
存在しない。

スピン端点を入れ替えない旧規則 (:math:`s = t` で :math:`-v`、
:math:`s \neq t` で :math:`+v`) は導出されたものではなく、数値が合う
ように固定された特殊ケースであった。固定 :math:`R` において
:math:`v[t,s] = -\overline{v[s,t]}` を満たす行列に限り両規則は同じ
行列を出力する。面内 Rashba のみの系はこの条件を満たすため差が現れない
が、スピン対称で実部と虚部をともに持つホッピングでは満たされない。
同じ再構成比較において旧規則の最大絶対偏差は :math:`6.0 \times 10^{-1}`
であった。

境界の wrap 位相は共役を取った **後** に適用する。反周期方向で境界を
跨ぐ行に対し、正の :math:`R` 越えで :math:`e^{+i\theta_d}`、負の
:math:`R` 越えで :math:`e^{-i\theta_d}` を掛ける。

.. note::

   検証済みの範囲は :math:`\theta` の各成分が :math:`\{0, \pi\}` の
   場合に限られる。この範囲では wrap 位相は実数である。一般の複素
   twist は対象外であり、検証もされていない。
