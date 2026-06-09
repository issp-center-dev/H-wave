==========================================
チュートリアル: FLEXソルバー
==========================================

本チュートリアルでは、H-waveのFLEX
(Fluctuation Exchange Approximation; ゆらぎ交換近似) ソルバーの使い方を
説明します。FLEXはRPAを拡張し、裸のグリーン関数の代わりに
ドレスド（自己無撞着な）グリーン関数を用いることで、
相関電子系のより正確な記述を提供します。

サンプルファイルは ``docs/ja/source/flex/sample`` ディレクトリにあります。


概要
----------------------------

FLEX近似 [1]_ は遍歴電子系のための自己無撞着ダイアグラム法です。
裸のグリーン関数 :math:`G_0` を用いるRPAと異なり、
FLEXは以下の自己無撞着ループを収束まで繰り返します:

1. Dyson方程式からドレスドグリーン関数 :math:`G(\mathbf{k}, i\omega_n)` を計算
2. ドレスド :math:`G` から裸感受率 :math:`\chi_0(\mathbf{q}, i\nu_m)` を計算
3. 相互作用をスピンチャネルと電荷チャネルに分解
4. スピン/電荷感受率 :math:`\chi_s`, :math:`\chi_c` のRPA方程式を解く
5. 有効相互作用 :math:`V_{\mathrm{eff}}` を構成
6. FFT畳み込みにより自己エネルギー :math:`\Sigma(\mathbf{k}, i\omega_n)` を計算
7. 収束判定; 未収束なら1に戻る


理論
----------------------------

ドレスドグリーン関数
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ドレスドグリーン関数はDyson方程式から得られます:

.. math::

   G(\mathbf{k}, i\omega_n)
   = \left[ G_0^{-1}(\mathbf{k}, i\omega_n) - \Sigma(\mathbf{k}, i\omega_n) \right]^{-1}

ここで :math:`G_0^{-1}(\mathbf{k}, i\omega_n) = i\omega_n + \mu - H_0(\mathbf{k})`
は裸のグリーン関数の逆です。

スピン・電荷感受率
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

裸感受率はドレスドグリーン関数から計算されます:

.. math::

   \chi_0(\mathbf{q}, i\nu_m) = -\frac{T}{N_k} \sum_{\mathbf{k}, n}
   G(\mathbf{k}+\mathbf{q}, i\omega_n + i\nu_m)\, G(\mathbf{k}, i\omega_n)

スピン感受率と電荷感受率は:

.. math::

   \chi_s = \left[ I - \chi_0 \, U_s \right]^{-1} \chi_0

.. math::

   \chi_c = \left[ I + \chi_0 \, U_c \right]^{-1} \chi_0

ここで :math:`U_s`, :math:`U_c` は全相互作用ハミルトニアンから
分解されたスピンおよび電荷相互作用頂点です。
単一バンドHubbardモデルでは :math:`U_s = U_c = U` となります。

有効相互作用
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FLEX有効相互作用はスピンゆらぎと電荷ゆらぎを結合します [1]_:

.. math::

   V_{\mathrm{eff}}(\mathbf{q}, i\nu_m)
   = W \left[ \frac{3}{2}\chi_s
            + \frac{1}{2}\chi_c - \chi_0 \right] W

ここで :math:`W` は裸の相互作用頂点です。
:math:`\chi_0` は **1回だけ** 減算します。最低次では :math:`\chi_s = \chi_c = \chi_0`
となるため括弧内は :math:`\chi_0` に帰着し、:math:`V_{\mathrm{eff}} = W \chi_0 W`
（2次の :math:`U^2` バブル）を与えます。この1回の減算により、:math:`\chi_s` と
:math:`\chi_c` に含まれる二重計上された2次ダイアグラムが除去されます。

自己エネルギー
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

自己エネルギーは実空間・虚時間での畳み込みにより計算されます:

.. math::

   \Sigma(\mathbf{r}, \tau) = V_{\mathrm{eff}}(\mathbf{r}, \tau) \cdot G(\mathbf{r}, \tau)

この要素ごとの（Hadamard）積はFFTを用いて効率的に評価されます。

.. [1] N. E. Bickers and D. J. Scalapino,
   Ann. Phys. (N.Y.) **193**, 206 (1989).


サンプル 1: 1軌道Hubbardモデル
-----------------------------------------

モデル
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

最初のサンプルは、二次元正方格子上のハーフフィリング (:math:`n = 1`)
における **1軌道Hubbardモデル** です。

.. math::

   H = -\sum_{\langle i,j \rangle, \sigma} t_{ij}\,
       c^\dagger_{i\sigma} c_{j\sigma}
     + U \sum_i n_{i\uparrow} n_{i\downarrow}

最近接ホッピング :math:`t = 1.0`、
次近接ホッピング :math:`t' = 0.5`、
オンサイトクーロン斥力 :math:`U = 4.0`、
温度 :math:`T = 0.5` を使用します。

このモデルは :math:`\mathbf{Q} = (\pi, \pi)` にピークを持つ
強い反強磁性 (AF) スピンゆらぎを示すことが知られています。

入力ファイルの準備
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

サンプルファイルは ``docs/ja/source/flex/sample/1orb/`` にあります。

**パラメータファイル** (``input.toml``):

.. literalinclude:: ../sample/1orb/input.toml

主要パラメータ:

- ``mode = "FLEX"``: FLEXソルバーを選択
- ``T = 0.5``: 温度
- ``CellShape = [8, 8, 1]``: 2D系の 8 x 8 k点メッシュ
- ``Nmat = 64``: 松原周波数の数
- ``filling = 0.5``: ハーフフィリング
- ``IterationMax = 100``: SCF反復の最大回数
- ``Mix = 0.2``: 自己エネルギー更新の混合パラメータ
  (:math:`\Sigma_{\mathrm{new}} = (1 - \alpha)\Sigma_{\mathrm{old}} + \alpha\Sigma_{\mathrm{calc}}`)
- ``EPS = 6``: 収束判定基準 :math:`10^{-6}`

**格子情報** (``geom.dat``):

.. literalinclude:: ../sample/1orb/geom.dat

原点に1つの軌道。

**トランスファー積分** (``transfer.dat``):

.. literalinclude:: ../sample/1orb/transfer.dat

正方格子上の最近接 (:math:`t = 1.0`) および次近接 (:math:`t' = 0.5`) ホッピング。

**オンサイト相互作用** (``coulombintra.dat``):

.. literalinclude:: ../sample/1orb/coulombintra.dat

オンサイトクーロン斥力 :math:`U = 4.0`。


計算の実行
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ cd docs/ja/source/flex/sample/1orb
    $ hwave input.toml

出力ログにSCFの収束過程が表示されます:

.. code-block:: text

    FLEX iteration 1/100
      convergence: |dSigma|/|Sigma| = 1.000e+00
    FLEX iteration 2/100
      convergence: |dSigma|/|Sigma| = 9.827e-01
    ...
    FLEX iteration 72/100
      convergence: |dSigma|/|Sigma| = 1.008e-06
    FLEX iteration 73/100
      convergence: |dSigma|/|Sigma| = 8.292e-07
    FLEX converged after 73 iterations


計算結果
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

収束後、ソルバーは ``output`` ディレクトリに以下の出力ファイルを生成します:

- ``chi0q.npz``: 裸感受率 :math:`\chi_0(\mathbf{q}, i\nu_m)`
- ``chiq_s.npz``: スピン感受率 :math:`\chi_s(\mathbf{q}, i\nu_m)`
- ``chiq_c.npz``: 電荷感受率 :math:`\chi_c(\mathbf{q}, i\nu_m)`
- ``chiq.npz``: 結合感受率ファイル
- ``sigma.npz``: 自己エネルギー :math:`\Sigma(\mathbf{k}, i\omega_n)`
- ``green.npz``: ドレスドグリーン関数 :math:`G(\mathbf{k}, i\omega_n)`

**スピン感受率** :math:`\chi_s(\mathbf{q}, i\nu_0)`:

.. figure:: ../sample/1orb/chi_s.png
   :width: 60%
   :align: center

   ハーフフィリングにおける1軌道Hubbardモデルの静的スピン感受率
   :math:`\chi_s(\mathbf{q})`。
   :math:`\mathbf{Q} = (\pi, \pi)` のピークはフェルミ面の
   ネスティングに起因する強い反強磁性スピンゆらぎを示しています。

**自己エネルギー** :math:`\mathrm{Im}\,\Sigma(\mathbf{k}, i\omega_0)`:

.. figure:: ../sample/1orb/sigma_kspace.png
   :width: 60%
   :align: center

   最低松原周波数における自己エネルギーの虚部。
   k依存性はスピンゆらぎによる準粒子の散乱を反映しており、
   反強磁性ホットスポット付近でより強いダンピングを示します。

**自己エネルギーの周波数依存性**:

.. figure:: ../sample/1orb/sigma_matsubara.png
   :width: 80%
   :align: center

   選択されたk点での自己エネルギーの周波数依存性。
   虚部 :math:`\mathrm{Im}\,\Sigma(i\omega_n) < 0` は
   準粒子のダンピングを示し、高周波数での
   :math:`1/\omega_n` テイルはフェルミ液体的振る舞いを示します。


サンプル 2: 2軌道Hubbardモデル
-----------------------------------------

モデル
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

2番目のサンプルは、軌道間クーロン相互作用とフント結合を含む
**2軌道Hubbardモデル** です:

.. math::

   H = \sum_{\mathbf{k},\alpha,\beta,\sigma}
       \varepsilon_{\alpha\beta}(\mathbf{k})\,
       c^\dagger_{\mathbf{k}\alpha\sigma} c_{\mathbf{k}\beta\sigma}
     + U \sum_{i,\alpha} n_{i\alpha\uparrow} n_{i\alpha\downarrow}
     + V \sum_{i,\alpha\neq\beta} n_{i\alpha} n_{i\beta}
     - 2J \sum_{i,\alpha\neq\beta}
       \mathbf{S}_{i\alpha} \cdot \mathbf{S}_{i\beta}

:math:`U = 4.0`, :math:`V = 1.0`, :math:`J = 0.5`,
温度 :math:`T = 1.0` を使用します。

入力ファイルの準備
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

サンプルファイルは ``docs/ja/source/flex/sample/2orb/`` にあります。

**パラメータファイル** (``input.toml``):

.. literalinclude:: ../sample/2orb/input.toml

1軌道サンプルとの主な違い:

- ``T = 1.0``: 安定性のための高い温度
- ``IterationMax = 200``: 多軌道系の収束のためのより多い反復回数
- 追加の相互作用ファイル: ``CoulombInter`` と ``Hund``

**格子情報** (``geom.dat``):

.. literalinclude:: ../sample/2orb/geom.dat

単位胞あたり2つの軌道。

**トランスファー積分** (``transfer.dat``):

.. literalinclude:: ../sample/2orb/transfer.dat

軌道内ホッピング (:math:`t = 1.0`) と軌道間混成 (:math:`t' = 0.5`)。

**オンサイト相互作用** (``coulombintra.dat``):

.. literalinclude:: ../sample/2orb/coulombintra.dat

両軌道の軌道内クーロン :math:`U = 4.0`。

**軌道間クーロン** (``coulombinter.dat``):

.. literalinclude:: ../sample/2orb/coulombinter.dat

軌道間クーロン :math:`V = 1.0`。

**フント結合** (``hund.dat``):

.. literalinclude:: ../sample/2orb/hund.dat

フント結合 :math:`J = 0.5`。


計算の実行
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ cd docs/ja/source/flex/sample/2orb
    $ hwave input.toml

.. code-block:: text

    FLEX iteration 1/200
      convergence: |dSigma|/|Sigma| = 1.000e+00
    FLEX iteration 2/200
      convergence: |dSigma|/|Sigma| = 3.587e-01
    ...
    FLEX iteration 58/200
      convergence: |dSigma|/|Sigma| = 1.188e-06
    FLEX iteration 59/200
      convergence: |dSigma|/|Sigma| = 9.684e-07
    FLEX converged after 59 iterations


計算結果
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**スピン感受率** :math:`\chi_s(\mathbf{q}, i\nu_0)`:

.. figure:: ../sample/2orb/chi_s.png
   :width: 60%
   :align: center

   2軌道モデルの静的スピン感受率。
   :math:`\mathbf{Q} = (\pi, \pi)` のピークはフント結合により増強されており、
   各サイト内での強磁性的な配列を促進しつつ、
   サイト間の反強磁性相関を許容しています。

**自己エネルギー** :math:`\mathrm{Im}\,\Sigma(\mathbf{k}, i\omega_0)`:

.. figure:: ../sample/2orb/sigma_kspace.png
   :width: 60%
   :align: center

   2軌道モデルの自己エネルギー虚部。
   軌道依存のk構造は多軌道系における異なる散乱チャネルを反映しています。

**自己エネルギーの周波数依存性**:

.. figure:: ../sample/2orb/sigma_matsubara.png
   :width: 80%
   :align: center

   2軌道モデルの自己エネルギーの周波数依存性。
   1軌道の場合と比較して大きな振幅は、
   軌道間相互作用による増強された相関を反映しています。


プロット
----------------------------

上記の図は以下のプロットスクリプトで再現できます:

.. code-block:: bash

    $ cd docs/ja/source/flex/sample
    $ python plot_results.py

または、個々のサンプルディレクトリ内で:

.. code-block:: bash

    $ cd docs/ja/source/flex/sample/1orb
    $ python ../plot_results.py


出力ファイル形式
----------------------------

FLEXソルバーは以下の内容を持つNumPy ``.npz`` ファイルを生成します:

``chi0q.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``chi0q``: 裸感受率 :math:`\chi_0(\mathbf{q}, i\nu_m)`,
  形状 ``(nmat, nvol, nd, nd)``
- ``freq_index``: 松原周波数インデックス
- ``wavevector_unit``: k点ベクトル
- ``wavevector_index``: 波数テーブル

``chiq_s.npz``, ``chiq_c.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``chiq_s`` / ``chiq_c``: スピン / 電荷感受率、
  ``chi0q`` と同じ形状

``sigma.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``sigma``: 自己エネルギー :math:`\Sigma(\mathbf{k}, i\omega_n)`,
  形状 ``(nblock, nmat, nvol, nd_block, nd_block)``
  (``nblock`` はスピンブロック数、spin-freeモードでは1)

``green.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``green``: ドレスドグリーン関数 :math:`G(\mathbf{k}, i\omega_n)`,
  ``sigma`` と同じ形状

これらの出力ファイルはEliashberg方程式ソルバー (``hwave_sc``) の
入力としても使用できます。詳細は :doc:`/rpa/tutorial/sc-index` を参照してください。


FLEX固有のパラメータ
----------------------------

FLEXソルバーは ``[mode.param]`` セクションで以下のパラメータを受け付けます:

.. list-table::
   :header-rows: 1
   :widths: 20 10 10 60

   * - パラメータ
     - 型
     - デフォルト
     - 説明
   * - ``IterationMax``
     - int
     - 100
     - SCF反復の最大回数
   * - ``Mix``
     - float
     - 0.2
     - 自己エネルギー更新の混合パラメータ :math:`\alpha`。
       小さい値はより安定な収束を与えますが、進行が遅くなります。
   * - ``EPS``
     - int/float
     - 6
     - 収束判定基準。整数 :math:`n` の場合、閾値は :math:`10^{-n}`。
       1未満の浮動小数点数の場合、直接閾値として使用。

その他のパラメータ (``T``, ``CellShape``, ``Nmat``, ``filling`` 等)
はRPAソルバーと共通です。詳細は :ref:`Ch:Config_rpa` を参照してください。

.. note::

   FLEXソルバーは縮約形の感受率を利用し、相互作用を密度-密度成分に
   縮約します。そのため ``calc_scheme = "reduced"`` または
   ``calc_scheme = "squashed"`` が必要であり、
   ``calc_type = "ring+ladder"`` （ ``"general"`` スキームを強制します）
   には対応していません。これらを満たさない場合、ソルバーは
   ``ValueError`` を送出します。


サンプル 3: 鉄系超伝導体2軌道モデル
-----------------------------------------

モデル
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

3番目のサンプルは、Raghu et al. [2]_ が提案した
**鉄系超伝導体の2軌道最小モデル** です。
Fe-As面を正方格子（1-Fe単位胞）上の :math:`d_{xz}` と
:math:`d_{yz}` 軌道で記述します。

.. math::

   H_0(\mathbf{k}) = \begin{pmatrix}
   \varepsilon_x(\mathbf{k}) & \varepsilon_{xy}(\mathbf{k}) \\
   \varepsilon_{xy}(\mathbf{k}) & \varepsilon_y(\mathbf{k})
   \end{pmatrix}

ここで

.. math::

   \varepsilon_x(\mathbf{k}) &= -2t_1 \cos k_x - 2t_2 \cos k_y - 4t_3 \cos k_x \cos k_y \\
   \varepsilon_y(\mathbf{k}) &= -2t_2 \cos k_x - 2t_1 \cos k_y - 4t_3 \cos k_x \cos k_y \\
   \varepsilon_{xy}(\mathbf{k}) &= -4t_4 \sin k_x \sin k_y

パラメータ: :math:`t_1 = -1.0`, :math:`t_2 = 1.3`,
:math:`t_3 = t_4 = -0.85`。

相互作用はKanamoriパラメータ化に従います:

.. math::

   H_{\mathrm{int}} = U \sum_{i,\alpha} n_{i\alpha\uparrow} n_{i\alpha\downarrow}
   + U' \sum_{i,\alpha\neq\beta} n_{i\alpha} n_{i\beta}
   - 2J \sum_{i,\alpha\neq\beta} \mathbf{S}_{i\alpha} \cdot \mathbf{S}_{i\beta}
   + J' \sum_{i,\alpha\neq\beta} c^\dagger_{i\alpha\uparrow} c^\dagger_{i\alpha\downarrow}
     c_{i\beta\downarrow} c_{i\beta\uparrow}

:math:`U = 1.5`, :math:`J = J' = 0.25`, :math:`U' = U - 2J = 1.0`,
温度 :math:`T = 0.1`、ハーフフィリング (:math:`n = 2`)。

フェルミ面は :math:`\Gamma` 点のホールポケットと
:math:`M = (\pi, 0)` / :math:`(0, \pi)` の電子ポケットから成ります。
これらのポケット間のネスティングが
:math:`\mathbf{Q} = (\pi, 0)` における強いスピンゆらぎを駆動します。
これが鉄系超伝導体の特徴的な物理です。

.. [2] S. Raghu, X.-L. Qi, C.-X. Liu, D. J. Scalapino, and S.-C. Zhang,
   Phys. Rev. B **77**, 220503(R) (2008).

入力ファイルの準備
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

サンプルファイルは ``docs/ja/source/flex/sample/iron_2orb/`` にあります。

**パラメータファイル** (``input.toml``):

.. literalinclude:: ../sample/iron_2orb/input.toml

**格子情報** (``geom.dat``):

.. literalinclude:: ../sample/iron_2orb/geom.dat

同一サイトに2つの軌道 (:math:`d_{xz}` と :math:`d_{yz}`)。

**トランスファー積分** (``transfer.dat``):

.. literalinclude:: ../sample/iron_2orb/transfer.dat

特徴的な2ポケットフェルミ面を生成するホッピングパラメータ。

**相互作用**:

.. literalinclude:: ../sample/iron_2orb/coulombintra.dat

.. literalinclude:: ../sample/iron_2orb/coulombinter.dat

.. literalinclude:: ../sample/iron_2orb/hund.dat

.. literalinclude:: ../sample/iron_2orb/exchange.dat


計算の実行
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ cd docs/ja/source/flex/sample/iron_2orb
    $ hwave input.toml

.. code-block:: text

    FLEX iteration 1/200
      convergence: |dSigma|/|Sigma| = 1.000e+00
    FLEX iteration 2/200
      convergence: |dSigma|/|Sigma| = 7.139e-01
    ...
    FLEX iteration 62/200
      convergence: |dSigma|/|Sigma| = 1.055e-06
    FLEX iteration 63/200
      convergence: |dSigma|/|Sigma| = 8.419e-07
    FLEX converged after 63 iterations


計算結果
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**スピン・電荷感受率**:

.. figure:: ../sample/iron_2orb/chi_spin_charge.png
   :width: 90%
   :align: center

   静的スピン感受率 :math:`\chi_s(\mathbf{q})`（左）と電荷感受率
   :math:`\chi_c(\mathbf{q})`（右）。
   スピン感受率は :math:`\mathbf{Q} = (\pi, 0)` および :math:`(0, \pi)` に
   ピークを持ち、ホールポケットと電子ポケット間のネスティングを反映しています。
   これは単一バンドHubbardモデル（ :math:`(\pi, \pi)` にピーク）とは
   定性的に異なります。

**軌道分解自己エネルギー**:

.. figure:: ../sample/iron_2orb/sigma_orbital.png
   :width: 90%
   :align: center

   最低松原周波数での軌道分解自己エネルギー虚部。
   :math:`d_{xz}` 軌道は :math:`k_y` 方向でより強い散乱を示し、
   :math:`d_{yz}` 軌道は :math:`k_x` 方向で強い散乱を示します。
   この軌道異方性はフェルミ面の軌道構造に起因します。

**自己エネルギーの周波数依存性**:

.. figure:: ../sample/iron_2orb/sigma_matsubara_orbital.png
   :width: 90%
   :align: center

   高対称k点での軌道分解自己エネルギーの周波数依存性。
   :math:`M = (\pi, 0)` において :math:`d_{yz}` 軌道は
   :math:`d_{xz}` よりも強くダンピングされ、
   軌道選択的な相関を反映しています。

**プロットスクリプト**:

.. code-block:: bash

    $ python plot_results.py


Tips
----------------------------

- **収束しない場合**: SCFループが収束しない場合は、``Mix`` を小さくする
  (例: 0.1 や 0.05)か、温度を上げてみてください。磁気不安定性に近い
  強相関領域では収束が困難になることがあります。

- **松原周波数**: 正確な結果を得るには十分な数の松原周波数 (``Nmat``)
  が必要です。目安として :math:`N_{\mathrm{mat}} \geq 10 / T` で
  低周波構造を捉えられます。

- **k点メッシュ**: メッシュサイズ (``CellShape``) は感受率と自己エネルギーの
  運動量構造を分解するのに十分大きくする必要があります。2D系では
  8x8で定性的な結果が得られ、定量的な計算には32x32以上が推奨されます。

- **計算コスト**: FLEXはSCFループのためRPAよりも計算コストが高くなります。
  コストは :math:`O(N_{\mathrm{iter}} \times N_k \times N_\omega \times N_d^3)`
  (ここで :math:`N_d = N_{\mathrm{orb}} \times N_{\mathrm{spin}}`)
  でスケールします。

- **Eliashberg方程式との連携**: FLEX出力ファイル (``chiq_s.npz``,
  ``chiq_c.npz``) は ``hwave_sc`` の ``[eliashberg]`` セクションで
  ``chi0q_mode = "flex"`` と設定することで使用できます。
  これによりFLEXレベルのスピン・電荷ゆらぎを用いた超伝導不安定性の
  解析が可能になります。


実装の詳細と制限事項
----------------------------

対応する相互作用型
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FLEXソルバーは以下の相互作用型に対応しています:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - 相互作用型
     - 対応
     - 備考
   * - ``CoulombIntra``
     - 対応
     - 軌道内クーロン斥力 :math:`U`
   * - ``CoulombInter``
     - 対応
     - 軌道間クーロン斥力 :math:`V`
   * - ``Hund``
     - 対応
     - フント結合 :math:`J`
   * - ``Exchange``
     - 部分的に対応
     - 交換相互作用 :math:`J'` （密度-密度成分のみ。
       スピンフリップ・対散乱などの非対角頂点は無視されます）
   * - ``Ising``
     - 対応
     - イジング型相互作用
   * - ``PairLift``
     - 部分的に対応
     - ペアリフト相互作用 （密度-密度成分のみ。
       非対角頂点は無視されます）
   * - ``PairHop``
     - 部分的に対応
     - ペアホッピング相互作用 （密度-密度成分のみ。
       非対角頂点は無視されます）
   * - ``InterAll``
     - **非対応**
     - 任意の4体相互作用（UHFrソルバーのみ対応）

.. note::

   ``InterAll`` 形式はk空間ソルバー (RPA/FLEX) では利用できません。
   ``InterAll`` で記述される相互作用は、上記の個別の相互作用型に
   分解して指定してください。


相互作用の運動量依存性（長距離相互作用）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

相互作用はWannier90形式の入力ファイルで指定します。
各行の形式は ``rx ry rz a b Re Im`` であり、
``(rx, ry, rz)`` は実空間の格子ベクトルです。

**オンサイト相互作用** (``rx = ry = rz = 0`` のみ):

FFT変換後、すべてのq点で同一の値を持つ
q非依存な相互作用 :math:`W(\mathbf{q}) = W_0` となります。
現在のサンプルはすべてこのケースです。

**長距離相互作用** (非零の ``(rx, ry, rz)`` を含む場合):

FFT変換により自動的にq依存の相互作用となります:

.. math::

   W(\mathbf{q}) = \sum_{\mathbf{r}} W(\mathbf{r})\, e^{-i\mathbf{q}\cdot\mathbf{r}}

例えば、最近接サイト間のクーロン相互作用を含める場合は、
相互作用ファイルに ``(1,0,0)`` や ``(0,1,0)`` 等の格子ベクトルを持つ
エントリを追加します。

.. note::

   連続的な :math:`1/r` クーロンポテンシャルを自動的に離散化する機能は
   ありません。各格子点での相互作用の値をユーザーが明示的に
   入力ファイルに指定する必要があります。


スピン・電荷チャネル分解の制約
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FLEXソルバーは、自己エネルギーの計算において相互作用テンソルを
スピンチャネルとチャネルに分解します。
この分解では、4体相互作用テンソル
:math:`W_{\alpha\sigma,\beta\sigma',\alpha\sigma,\beta\sigma'}`
を2体の縮約形式 :math:`W_{\alpha\sigma,\beta\sigma'}` に変換します。

具体的には、同一スピン成分と異スピン成分に分離されます:

.. math::

   U_s &= W_{\mathrm{cross}} - W_{\mathrm{same}} \\
   U_c &= W_{\mathrm{cross}} + W_{\mathrm{same}}

ここで :math:`W_{\mathrm{same}}` は同一スピン間、
:math:`W_{\mathrm{cross}}` は異スピン間の相互作用です。

この縮約は **密度-密度型相互作用** に対して正確です。
``CoulombIntra``, ``CoulombInter``, ``Hund``, ``Ising``
はすべて密度-密度型であり正しく扱われます。
一方、 ``Exchange``, ``PairLift``, ``PairHop`` は純粋な密度-密度型では
**なく**、縮約においては密度-密度成分のみが保持され、非対角
（スピンフリップ・対散乱）頂点は無視されます。これらの相互作用が
含まれる場合、ソルバーはこの近似をユーザーに知らせるための警告を
出力します。


スピン自由度の扱い（spin-freeモード）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FLEXソルバーはデフォルトで **spin-freeモード** で動作します。
このモードでは SU(2) スピン対称性を仮定し、
計算量を削減します。

spin-freeモードでは:

- グリーン関数は軌道空間のみで表現されます
  (形状: ``(1, nmat, nvol, norb, norb)``)
- 感受率や有効相互作用は内部的にスピン軌道空間
  (``nd = norb × ns``) に膨張されて計算されます
- 自己エネルギーの計算後、SU(2)対称性により保証される
  :math:`\Sigma_{\uparrow\uparrow} = \Sigma_{\downarrow\downarrow}`,
  :math:`\Sigma_{\uparrow\downarrow} = 0` の性質を利用して
  軌道空間に縮約されます

.. note::

   spin-freeモードは常磁性状態（磁気秩序のない状態）を前提としています。
   磁気秩序相の記述にはスピン自由度を陽に扱う拡張が必要です。
