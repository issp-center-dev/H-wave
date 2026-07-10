=============================================
チュートリアル: Eliashberg方程式ソルバー
=============================================

本チュートリアルでは、H-waveに含まれる線形化Eliashberg方程式ソルバー
``hwave_sc`` の使い方を説明します。
このツールは、RPAソルバーで計算した裸感受率 :math:`\chi_0(\mathbf{q})`
を用いて線形化Eliashberg方程式を解き、超伝導不安定性を解析します。

サンプルファイルは ``docs/ja/source/rpa/sample_sc`` ディレクトリにあります。


計算の流れ
----------------------------

計算は2つのステップで行います:

1. **RPA計算** (``hwave``): 裸感受率 :math:`\chi_0(\mathbf{q})` を計算し、
   ``chi0q.npz`` に保存します。
2. **Eliashberg方程式ソルバー** (``hwave_sc``): ``chi0q.npz`` を読み込み、
   グリーン関数を再構成し、RPA頂点を計算した後、
   線形化Eliashberg方程式を解きます。


モデル
----------------------------

本チュートリアルでは、二次元正方格子上の **2軌道タイトバインディングモデル**
を3/4フィリングで扱います。ハミルトニアンは以下の通りです:

.. math::

   H = \sum_{\mathbf{k},\alpha,\beta,\sigma}
       \varepsilon_{\alpha\beta}(\mathbf{k})\,
       c^\dagger_{\mathbf{k}\alpha\sigma} c_{\mathbf{k}\beta\sigma}
     + U \sum_{i,\alpha} n_{i\alpha\uparrow} n_{i\alpha\downarrow}
     + \sum_{i,\alpha\neq\beta} V_{\alpha\beta}\,
       n_{i\alpha} n_{i\beta}

ここで、オンサイトクーロン斥力 :math:`U = 0.4`、
軌道間クーロン相互作用 :math:`V` を使用します。

このサンプルは有機導体 :math:`\beta`\ -(meso-DMBEDT-TTF)\ :math:`_2`\ PF\ :math:`_6`
の伝導層モデルに基づいており、
トランスファー積分は拡張ヒュッケル法による計算値を使用しています。 [1]_


理論
----------------------------

超伝導感受率
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RPA電荷感受率 :math:`\hat{X}^c` およびスピン感受率 :math:`\hat{X}^s` は
以下のように与えられます:

.. math::

   \hat{X}^c = (\hat{I} + \hat{X}^{(0)} (\hat{U} + 2\hat{V}))^{-1} \hat{X}^{(0)}

.. math::

   \hat{X}^s = (\hat{I} - \hat{X}^{(0)} \hat{U})^{-1} \hat{X}^{(0)}

ここで :math:`\hat{X}^{(0)}` は裸感受率、
:math:`\hat{U}` はオンサイト相互作用行列、
:math:`\hat{V}` はサイト間相互作用行列です。

線形化Eliashberg方程式
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

一重項超伝導の線形化Eliashberg方程式は次のように書けます:

.. math::

   \lambda_S \Sigma^a_{\alpha\sigma;\beta\bar{\sigma}}(\mathbf{k})
   = -\frac{T}{N_L} \sum_{\mathbf{k}',n',\alpha',\beta'}
   P^S_{\alpha\sigma;\beta\bar{\sigma}}(\mathbf{k} - \mathbf{k}')
   G^{(0)}_{\alpha\alpha'}(\mathbf{k}', i\varepsilon_{n'})
   G^{(0)}_{\beta\beta'}(-\mathbf{k}', -i\varepsilon_{n'})
   \Sigma^a_{\alpha'\sigma;\beta'\bar{\sigma}}(\mathbf{k}')

ペアリング相互作用は、一重項の場合:

.. math::

   \hat{P}^S = \hat{U} + \hat{V}
   + \frac{3}{2} \hat{U} \hat{X}^s \hat{U}
   - \frac{1}{2} (\hat{U} + 2\hat{V}) \hat{X}^c (\hat{U} + 2\hat{V})

三重項の場合:

.. math::

   \hat{P}^T = \hat{V}
   - \frac{1}{2} \hat{U} \hat{X}^s \hat{U}
   - \frac{1}{2} (\hat{U} + 2\hat{V}) \hat{X}^c (\hat{U} + 2\hat{V})

:math:`\lambda_S = 1` (:math:`\lambda_T = 1`) が超伝導転移点に対応します。
:math:`\lambda > 1` （正の固有値）のとき常伝導状態は超伝導に対して不安定です。
負の固有値は符号反転ギャップに対応しますが、
自己無撞着条件 :math:`\Delta = K\Delta` を満たさないため
超伝導不安定性を示しません。

Eliashberg方程式の数値解法として、``hwave_sc`` では自己無撞着べき乗反復法
（固有値が最大のモードに収束）とArnoldi法による固有値解析を実装しています。

.. [1] K. Yoshimi, M. Nakamura, and H. Mori,
   J. Phys. Soc. Jpn. **76**, 024706 (2007);
   `arXiv:cond-mat/0608466 <https://arxiv.org/abs/cond-mat/0608466>`_.


入力ファイルの準備
----------------------------

パラメータファイル
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TOML形式のパラメータファイル ``input.toml`` を作成します:

.. literalinclude:: ../sample_sc/input.toml

このファイルは以下のセクションから構成されます。

``[mode.param]`` セクション
""""""""""""""""""""""""""""""""

- ``T``: 温度。
- ``CellShape``: k点メッシュのサイズ（2次元系では 32 x 32 x 1）。
- ``Nmat``: 松原振動数の数（512）。
- ``filling``: 1軌道1スピンあたりの電子充填率（0.75 = 3/4フィリング）。

``[file]`` セクション
""""""""""""""""""""""""""""""""

- ``[file.input.interaction]``: 幾何情報、トランスファー積分、
  相互作用パラメータのファイルを指定します。RPAステップと共通です。
- ``[file.output]``: :math:`\chi_0(\mathbf{q})` および
  :math:`\chi(\mathbf{q})` の出力先ディレクトリとファイル名。

``[eliashberg]`` セクション
""""""""""""""""""""""""""""""""

Eliashberg方程式ソルバーの設定です。主なパラメータ:

- ``solver_mode``: ``"iteration"`` (自己無撞着べき乗法)、
  ``"eigenvalue"`` (Arnoldi固有値解析)、または ``"both"``。
- ``chi0q_mode``: ``"load"`` はRPA出力ファイルから
  :math:`\chi_0(\mathbf{q})` を読み込みます。
  ``"calc"`` は内部で計算します。
  ``"flex"`` はFLEX計算の dressed 感受率を読み込みます
  （``frequency = "dynamic"`` に必須）。
- ``frequency``: ペアリング頂点の振動数の扱い。``"static"`` （デフォルト）は
  ボゾン振動数ゼロでペアリング頂点を評価する静的近似（Nakano--Kuroki 式(9)）で、
  振動数依存性のないギャップを与えます。``"dynamic"`` は、振動数依存のペアリング
  頂点 :math:`V(\mathbf{q}, i\omega_l)` とギャップ :math:`\phi(\mathbf{k}, i\omega_n)`
  を用いた、松原振動数に完全に依存する Eliashberg 方程式を解きます。
  ``chi0q_mode = "flex"`` が必要です
  （下記の :ref:`動的振動数の節 <sc_dynamic_frequency>` を参照）。
- ``pairing_type``: ``"singlet"`` または ``"triplet"``。
- ``init_gap``: 反復法の初期ギャップ対称性。
  ``"cos"`` (:math:`\cos(k_x+k_y+k_z)`)、
  ``"d_x2y2"`` (:math:`\cos k_x - \cos k_y`)、
  ``"random"`` などが利用可能です。
  有効な形状因子の一覧は
  ``"cos"`` 、 ``"s"`` 、 ``"s_ext"`` 、 ``"s_ext_2d"`` 、 ``"d_x2y2"`` 、
  ``"d_xy"`` 、 ``"d_xz"`` 、 ``"d_yz"`` 、 ``"d_z2"`` 、
  ``"p_x"`` 、 ``"p_y"`` 、 ``"p_z"`` 、 ``"random"`` です。
- ``max_iter``: 自己無撞着反復の最大回数。
- ``alpha``: 混合パラメータ（0: 混合なし、1: 古い解を完全保持）。
- ``convergence_tol``: ギャップ関数の収束条件。
- ``num_eigenvalues``: 固有値モードで計算する固有値の数。
- ``eigenvalue_method``: ``"arnoldi"`` （デフォルト）、 ``"subspace"`` 、
  ``"shift-invert-gmres"`` / ``"shift-invert-bicgstab"`` /
  ``"shift-invert-lgmres"``。
- ``gpu``: ``true`` で動的モード（``frequency = "dynamic"``）のカーネル適用を
  GPU（CuPy）で実行します（デフォルト ``false``。下記の
  :ref:`GPU実行の節 <sc_dynamic_gpu>` を参照）。
- ``fft_workers``: 動的モードの空間 FFT のワーカースレッド数
  （デフォルト ``1`` = 従来どおりの直列 numpy。``-1`` で全コア。
  GPU 実行時は無視されます）。
- ``matsubara_basis``: 動的モードの松原軸の表現。``"uniform"``（デフォルト、
  従来どおり）または ``"ir"``（sparse-ir による中間表現基底。下記の
  :ref:`IR基底の節 <sc_dynamic_ir>` を参照）。
- ``ir_tol``: IR基底の打ち切り精度 :math:`\varepsilon`（デフォルト 1e-8）。
- ``ir_wmax``: IR基底の実周波数バンド幅 :math:`\omega_{\max}`（省略時は
  バンド幅と相互作用スケールから自動推定。推定できない場合は明示指定を
  求めるエラーになります）。

相互作用定義ファイル
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

相互作用定義ファイルはWannier90形式で記述し、
RPAソルバーと共通です。詳細は :ref:`Ch:Config_rpa` を参照してください。

``Geometry`` (``geom.dat``):

.. literalinclude:: ../sample_sc/geom.dat

2軌道の単位胞を定義します。

``Transfer`` (``transfer.dat``):

.. literalinclude:: ../sample_sc/transfer.dat

2軌道モデルのホッピング積分を定義します。

``CoulombIntra`` (``coulombintra.dat``):

.. literalinclude:: ../sample_sc/coulombintra.dat

各軌道のオンサイトクーロン斥力 :math:`U = 0.4`。

``CoulombInter`` (``coulombinter.dat``):

.. literalinclude:: ../sample_sc/coulombinter.dat

軌道間・サイト間クーロン相互作用。


ステップ1: RPA計算の実行
----------------------------

まず、RPAソルバーを実行して裸感受率を計算します:

.. code-block:: bash

    $ hwave input.toml

``output/chi0q.npz`` と ``output/chiq.npz`` が生成されます。
32 x 32メッシュの場合、数秒で完了します。


ステップ2: Eliashberg方程式ソルバーの実行
---------------------------------------------

次に、同じ入力ファイルを使ってEliashberg方程式ソルバーを実行します:

.. code-block:: bash

    $ hwave_sc input.toml

ソルバーは以下の処理を行います:

1. ``output/chi0q.npz`` から :math:`\chi_0(\mathbf{q})` を読み込みます。
2. 相互作用ファイルを読み込み、ハミルトニアンを構築します。
3. 非相互作用グリーン関数 :math:`G(\mathbf{k}, i\omega_n)` を構成します。
4. RPA電荷・スピン頂点 :math:`V_c(\mathbf{q})`、
   :math:`V_s(\mathbf{q})` を計算します。
5. 線形化Eliashberg方程式を自己無撞着反復法および/または
   固有値解析で解きます。

実行ログの例:

.. code-block:: text

    hwave_sc: === Self-consistent iteration ===
    hwave_sc: Iteration    0: eigenvalue = 0.924446, diff = 3.544353e-01
    hwave_sc: Iteration    1: eigenvalue = 0.817270, diff = 7.893848e-02
    ...
    hwave_sc: Iteration  192: eigenvalue = 0.959725, diff = 9.900091e-06
    hwave_sc: Converged at iteration 193
    hwave_sc: Iteration result: eigenvalue = 0.959725, converged = True, n_iter = 193

引き続き固有値解析の結果が表示されます:

.. code-block:: text

    hwave_sc: === Eigenvalue analysis ===
    hwave_sc: Leading eigenvalues:
    hwave_sc:     0: 0.959725 (|ev| = 0.959725)
    hwave_sc:     1: 0.836778 (|ev| = 0.836778)
    hwave_sc:     2: 0.810959 (|ev| = 0.810959)
    hwave_sc:     3: -0.887954 (|ev| = 0.887954)
    hwave_sc:     4: -1.071303 (|ev| = 1.071303)
    hwave_sc:     5: -1.349775 (|ev| = 1.349775)
    hwave_sc:     6: 0.976353 (|ev| = 0.976353)  [opposite-parity sector]
    hwave_sc:     7: 0.823774 (|ev| = 0.823774)  [opposite-parity sector]
    ...

正の固有値 :math:`\lambda > 1` は、その温度で超伝導不安定性が
存在することを示します。負の固有値は符号反転ギャップに対応しますが、
超伝導不安定性は示しません。

固有対は **チャネルのパリティ**\ （一重項=偶、三重項=奇）を持つものが
先頭に並べられます。先頭（index 0）が物理解で、自己無撞着反復の結果と
一致します。``[opposite-parity sector]`` のタグが付いた固有値は逆パリティの
セクターに属し、そのスピンチャネルではパウリ原理により禁止され、物理的な
不安定性を表しません（下記のパリティに関する注記を参照）。
``eigenvalue.dat`` の末尾列 ``match`` がこれを記録します（``1`` = チャネルのパリティ
セクター、``0`` = 逆パリティセクター）。


計算結果
----------------------------

ギャップ関数
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

以下の図は、自己無撞着反復から得られた
ギャップ関数 :math:`\Sigma_{\alpha\beta}(\mathbf{k})`
の運動量空間での分布を示しています。

**一重項チャネル** (:math:`\lambda \approx 0.96`):

.. figure:: ../sample_sc/gap_singlet.png
   :width: 90%
   :align: center

   一重項ギャップ関数のk空間分布。
   左: 軌道内成分 :math:`\mathrm{Re}\,\Sigma_{00}(\mathbf{k})`。
   右: 軌道間成分 :math:`\mathrm{Re}\,\Sigma_{01}(\mathbf{k})`。
   軌道間成分が軌道内成分の約5倍大きく、
   軌道間ペアリングが支配的であることを示す。

**三重項チャネル** (:math:`\lambda \approx 0.97`):

.. figure:: ../sample_sc/gap_triplet.png
   :width: 90%
   :align: center

   三重項ギャップ関数のk空間分布。
   左: 軌道内成分 :math:`\mathrm{Re}\,\Sigma_{00}(\mathbf{k})`。
   右: 軌道間成分 :math:`\mathrm{Re}\,\Sigma_{01}(\mathbf{k})`。
   ギャップは :math:`\mathbf{k} \to -\mathbf{k}` （および軌道の入れ替え）に対して
   **奇** であり、スピン三重項ペアリングに要請される対称性を満たす。
   軌道間成分が再び大きい。

固有値スペクトル
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

一重項・三重項チャネルのEliashberg方程式の固有値スペクトルを以下に示します。

.. figure:: ../sample_sc/eigenvalue_spectrum.png
   :width: 70%
   :align: center

   線形化Eliashberg方程式の正の固有値スペクトル :math:`\lambda` 。
   赤破線は :math:`\lambda = 1` （超伝導不安定性の判定基準）を示す。
   塗りつぶしマーカーは **物理的** な固有値（ギャップがチャネルのパリティ
   ＝一重項は偶・三重項は奇を持つ）、白抜きマーカーは **逆パリティの spurious**
   モードである。物理的な固有値はすべて1未満。1を超える2つの白抜きマーカーは
   *三重項* カーネルの偶パリティ解で、スピン三重項ペアリングにはパウリ原理により
   禁止される（下記のパリティに関する注記を参照）。

Arnoldi固有値解析は複数の固有値を検出します。
図には超伝導不安定性の判定基準 :math:`\lambda = 1` に関連する
正の固有値のみを示しています。
一重項チャネルの主要な **物理的**\ （偶パリティ）固有値は
:math:`\lambda_S \approx 0.96 < 1` で、自己無撞着反復法の結果と一致し、
この温度で一重項SC不安定性は存在しません。

三重項チャネルの主要な **物理的**\ （奇パリティ）固有値は
:math:`\lambda_T \approx 0.97 < 1` です。三重項の計算に現れる
:math:`1.58` や :math:`1.09` 付近の大きな値は偶パリティの spurious モードであり、
三重項不安定性では **ありません**\ 。2つの物理的チャネルは非常に近く
（ :math:`\lambda_T \gtrsim \lambda_S` ）、この温度が一重項・三重項の
クロスオーバー近傍にあることを示します。実際の転移（ :math:`\lambda = 1` ）は
より低温で到達します。

描画スクリプト
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

上記の図はサンプルディレクトリに含まれる描画スクリプトで再現できます:

.. code-block:: bash

    $ python plot_results.py


出力ファイル
----------------------------

ソルバーは ``output`` ディレクトリに以下のファイルを出力します。

``gap.dat``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

収束したギャップ関数 :math:`\Delta_{\alpha\beta}(\mathbf{k})`
のk空間表示。各行の形式:

.. code-block:: text

    kx  ky  kz  Re(Δ_00)  Im(Δ_00)  Re(Δ_01)  Im(Δ_01)  Re(Δ_10)  Im(Δ_10)  Re(Δ_11)  Im(Δ_11)

ここで :math:`\alpha, \beta` は軌道インデックスです。

``eigenvalue.dat``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

線形化Eliashberg方程式の固有値:

.. code-block:: text

    # Iteration eigenvalue
    9.59724792e-01
    # Eigenvalue analysis
    # index  Re(eigenvalue)  Im(eigenvalue)  |eigenvalue|  match(1=channel-parity)
       0  9.59725061e-01 -3.63739178e-12  9.59725061e-01 1
       1  8.36778019e-01 -3.64670431e-12  8.36778019e-01 1
       ...

末尾の ``match`` 列は、固有ベクトルがチャネルのパリティを持つ（物理的な
一重項／三重項ギャップ）とき ``1`` 、逆パリティの spurious モードのとき ``0`` です
（上記のパリティに関する注記を参照）。この列を持たない旧形式の出力もそのまま
読み込めます。


物理的解釈
----------------------------

線形化Eliashberg方程式の最大固有値 :math:`\lambda` は
超伝導の発生を判定します:

- :math:`\lambda > 1`: 常伝導状態は超伝導に対して不安定です。
  対応する固有ベクトルがギャップ関数の対称性を与えます。
- :math:`\lambda < 1` （全ての正の固有値について）:
  その温度では常伝導状態が安定です。

負の固有値は、多軌道系における :math:`s_\pm` 波のような
符号反転ペアリング対称性に対応しますが、
自己無撞着条件 :math:`\Delta = K\Delta` は :math:`\lambda = 1`
（ :math:`\lambda = -1` ではない）を要求するため、
負の固有値はその大きさによらず超伝導不安定性を **示しません** 。

温度を変化させて最大の正の固有値が :math:`\lambda = 1` となる点を
求めることで、超伝導転移温度 :math:`T_c` を決定できます。

本チュートリアルでは、:math:`T = 0.1` で両方の物理的チャネルが不安定性の
しきい値のわずか下にあり、:math:`\lambda_S \approx 0.96` （一重項・偶）、
:math:`\lambda_T \approx 0.97` （三重項・奇）です。両者はほぼ縮退しており、
三重項がわずかに先行します。

.. note::

   **パリティとパウリ原理。**
   クーパー対は2電子の交換（スピン交換・軌道交換・ :math:`\mathbf{k} \to -\mathbf{k}`
   の同時操作）に対して反対称でなければなりません。（偶周波数の）ギャップ
   :math:`\Sigma_{\alpha\beta}(\mathbf{k})` では、これが空間パリティ
   :math:`P:\ \Sigma_{\alpha\beta}(\mathbf{k}) \to \Sigma_{\beta\alpha}(-\mathbf{k})`
   を固定します。スピン一重項ギャップは **偶**\ （ :math:`P = +1` ）、スピン三重項
   ギャップは **奇**\ （ :math:`P = -1` ）です。Eliashbergカーネルはギャップ全空間に
   作用し片方のパリティに制限されないため、各チャネルのカーネルは逆パリティの
   固有ベクトルも持ちます。これらはそのスピンチャネルでは物理的意味を持たない
   数学的解です（例えば三重項カーネルの偶パリティ固有値は全対称な対状態に相当し、
   パウリ原理で禁止されます）。そのため ``hwave_sc`` は各反復をチャネルのパリティ
   セクターへ射影し、固有値解析は各モードを物理的（ ``match = 1`` ）か spurious
   （ ``match = 0`` ）かでラベル付けします。三重項の計算に現れる
   :math:`\lambda \approx 1.58,\ 1.09` という大きな値はこのような spurious な
   偶パリティモードであり、三重項不安定性では **ありません** 。

一重項と三重項の比較
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

入力ファイルの ``pairing_type`` を ``"triplet"`` に変更することで、
一重項と三重項のチャネルを比較できます。
三重項チャネルに切り替える際は、``init_gap`` も ``"p_x"`` のような
奇パリティ（三重項）の初期ギャップに設定してください（同梱の
``input_triplet.toml`` はそのように設定済みです）。既定の
``init_gap = "cos"`` は偶パリティ（一重項）の初期ギャップであり、
ソルバーはパリティの異なる初期ギャップをエラーとして拒否するようになりました
（非物理的なセクターへ収束するのを防ぐため）。あるいは ``init_gap`` を
省略すれば、ソルバーがチャネルのパリティを自動選択します。
同じパラメータで :math:`T = 0.1` の場合、主要な物理的固有値は
:math:`\lambda_T \approx 0.97` （三重項）と :math:`\lambda_S \approx 0.96`
（一重項）で、ともに1未満、三重項がわずかに先行します。これは本サンプルが
文献 [1]_ の一重項・三重項クロスオーバー近傍にあることを示します。
同文献では、:math:`T > 0.05` で三重項SC状態が一重項SC状態と競合し、低温
(:math:`T < 0.05`) ではスピンゆらぎの増大により一重項SC転移が支配的になることが
報告されています。実際の転移（ :math:`\lambda = 1` ）は低温で到達します。


.. _sc_dynamic_frequency:

動的（振動数依存）Eliashberg方程式
--------------------------------------------------

``hwave_sc`` は既定では **静的近似** で Eliashberg 方程式を解きます。すなわち
ペアリング頂点をボゾン松原振動数ゼロで評価し（Nakano--Kuroki 式(9) の静的近似）、
ギャップ :math:`\Sigma(\mathbf{k})` は振動数依存性を持ちません。``[eliashberg]``
セクションで ``frequency = "dynamic"`` と設定すると、代わりに **振動数に完全に
依存する** 線形化 Eliashberg 方程式

.. math::

   \lambda\, \phi_{\alpha\beta}(\mathbf{k}, i\omega_n)
   = -\frac{T}{N_L} \sum_{\mathbf{k}', n'}
     V_{\alpha\beta}(\mathbf{k}-\mathbf{k}', i\omega_n - i\omega_{n'})\,
     [G G](\mathbf{k}', i\omega_{n'})\,
     \phi(\mathbf{k}', i\omega_{n'}),

を解きます。ここではギャップ :math:`\phi(\mathbf{k}, i\omega_n)` のフェルミオン
松原軸全体と、振動数依存のペアリング頂点 :math:`V(\mathbf{q}, i\omega_l)` を保持
します。頂点は（静的な係数ではなく）虚時間の積として作用するため、カーネルは
異なる松原振動数どうしを結合します。

FLEXの前提条件
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

動的モードは FLEX ソルバーのみが生成する振動数分解された入力を必要とするため、
``chi0q_mode = "flex"`` で実行する **必要があります**。``hwave_sc`` を呼び出す前に、
Eliashberg ステップが読み込むディレクトリへ以下を書き出す FLEX 計算
（``mode = "FLEX"``）を実行してください:

- ``chiq_s.npz`` と ``chiq_c.npz`` -- **全** ボゾン松原軸（``Nmat`` 個の全振動数）
  上のスピン・電荷感受率、および
- ``green.npz`` -- **dressed** グリーン関数
  :math:`G(\mathbf{k}, i\omega_n)`\ （これからペアバブルを構成します）。

.. note::

   FLEX の感受率ファイルには ``chi_convention`` タグ（reduced/squashed スキームは
   ``"kuroki"``\ 、general full-vertex スキームは ``"myo"``\ ）が付与されており、
   Eliashberg ローダーはこれを用いて軌道レイアウトを解釈します。\ **2軌道系**\
   （\ ``norb = 2``\ ）では reduced のスピン軌道次元と orbital-pair 次元が一致
   （ともに ``4``\ ）するため、両者の区別はこのタグに依存します。本修正以前の H-wave
   は形状のみからレイアウトを推定し、``norb = 2`` の reduced (kuroki) chi を
   orbital-pair と誤認して pairing vertex を壊していました。該当する計算の Eliashberg
   固有値・ギャップ関数は本バージョンで\ **修正されます（したがって変化します）**\ 。
   1軌道系および general (myo) の結果は影響を受けません。

``Nmat`` は偶数で、かつ FLEX 出力と ``[mode.param]`` の値とで一致していなければ
なりません。``chi0q_mode = "flex"`` なしで ``frequency = "dynamic"`` を指定した
場合、``Nmat`` が奇数の場合、あるいは dressed ``green.npz`` が無い場合、ソルバーは
黙って別の処理へフォールバックせず、説明的なエラーで停止します。FLEX 出力
ディレクトリは ``[file.input] path_to_flex_output`` で指定でき（既定は
``[file.output]`` ディレクトリ）、個々のファイル名は ``[eliashberg]`` の
``flex_chi_s`` / ``flex_chi_c`` / ``flex_green`` で上書きできます。

出力
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

動的モードは ``eigenvalue.dat`` （最大固有値 :math:`\lambda`）に加えて、
以下を書き出します:

``gap_dynamic.npz``
   振動数分解されたギャップ全体とそのメタデータ。キー:

   - ``gap``: 形状 ``(norb, norb, Nx, Ny, Nz, Nmat)`` の複素配列 --
     :math:`\phi_{\alpha\beta}(\mathbf{k}, i\omega_n)`。
   - ``iomega``: 中心化されたフェルミオン松原振動数
     :math:`\omega_n = (2n + 1 - N_{\mathrm{mat}})\pi T`。
   - ``T``: 温度。
   - ``pairing_type``: ``"singlet"`` または ``"triplet"``。
   - ``frequency``: ``"dynamic"``。
   - ``eigenvalue``: 最大固有値 :math:`\lambda`。
   - ``axis_order``: ``"(orb1, orb2, kx, ky, kz, iomega)"``。
   - ``normalization``: ゲージ規約 -- ギャップは全成分について L2 規格化され、
     最大絶対値の成分が実正になるよう回転されます。これにより保存されるギャップは
     実行ごと・線形代数バックエンドごとに再現可能です。

``gap.dat``
   最小の正の松原振動数（インデックス ``Nmat//2``）におけるギャップの単一振動数
   スライス。列の並びは静的な ``gap.dat`` と同じ（``kx ky kz`` の後に軌道対ごとの
   ``Re``/``Im``）です。先頭行は ``#`` で始まるヘッダで、``frequency=dynamic`` と
   スライスのインデックス・その :math:`\omega_n` を記録します。

チャネルパリティによる選別
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

静的ソルバーと同様に、動的モードも単に代数的に最大の固有値ではなく、**指定した
ペアリングチャネルの**先頭固有対を報告します。フェルミオンの反対称性により、
ギャップの結合パリティ :math:`\phi_{\alpha\beta}(\mathbf{k}, i\omega_n) \to
\phi_{\beta\alpha}(-\mathbf{k}, -i\omega_n)` が決まります（``singlet`` では偶
――従来型の偶振動数解と奇振動数解の両方を許容――、``triplet`` では奇）。Arnoldi
の固有対はチャネルパリティのモードが先頭に来るよう並べ替えられ、
``eigenvalue.dat`` の固有値表には静的出力と同じ末尾列
``match(1=channel-parity)`` が付きます（指定セクターで ``1``、反対セクターで
``0``）。計算した ``num_eigenvalues`` 個の固有対のいずれも指定セクターに無い
場合は、警告を出して生の先頭固有対にフォールバックします。その際は
``num_eigenvalues`` を増やすか ``pairing_type`` を確認してください。べき乗反復
（``solver_mode = "iteration"``）でも、カーネルがパリティと可換な（中心対称な）
系では各反復ベクトルをチャネルセクターへ射影します。可換でない場合は警告を出して
射影を無効化し、射影なしの反復を用います。

メモリに関する注意
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

動的ソルバーは複数の全振動数テンソル（ペアリング頂点・ペアバブル・ギャップ）を
保持するため、ピークメモリはおおよそ
:math:`\mathcal{O}(N_{\mathrm{orb}}^4\, N_k\, N_{\mathrm{mat}})` で増加し、軌道数・
k メッシュ・松原振動数の数に対して急速に大きくなります。``hwave_sc`` は確保前に
ピーク必要量を見積もり、上限を超える場合は停止します。``[eliashberg] mem_limit_gb``
で上限を明示的に指定でき（``0`` でガードを無効化）、指定しない場合は利用可能な
RAM の一部が上限に使われます。この見積もりは、設定値だけを信頼せず、ディスク上のファイルヘッダ
（``chiq_s.npz`` / ``chiq_c.npz`` / ``green.npz``\ ）から保存された ``Nmat``
（松原振動数軸）を読み取ります。そのため、ファイルに保存された ``Nmat`` が設定と
異なる場合は、ロード途中で OOM で落ちるのではなく、確保前に事前に拒否されます
（k メッシュと軌道数は設定から取得し、その不一致はローダーの reshape エラーとして
現れます）。

性能に関する注意
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

固有値ソルブはカーネルの繰り返し適用が支配的です。2つの最適化でこれを軽く
保ちます：頂点の虚時間変換（カーネルで最も重いステップ）を matvec ごとではなく
一度だけ事前計算し、空間 FFT を ``scipy.fft`` で並列化します。``[eliashberg]
fft_workers`` は FFT のワーカースレッド数を指定します：``1``（デフォルト）は
従来どおりの直列 numpy、``-1`` は全コアを使用します。既存の計算を変えないよう
オプトインです。複数の動的計算を同時に走らせる場合は、CPU の過剰割当てを避ける
ため小さめの値（例：``OMP_NUM_THREADS`` に合わせる）に設定してください。両者を合わせ
て ``norb = 2``、``N_k = 1024``、``N_{mat} = 1024`` でおよそ 4 倍の高速化になります。
GPU 実行（``gpu = true``）では FFT は GPU 上で走るため ``fft_workers`` は無視されます。

.. _sc_dynamic_ir:

IR基底（sparse-ir）による松原軸の圧縮
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``matsubara_basis = "ir"`` を指定すると、動的モードの松原軸を一様グリッド
（``Nmat`` 点）から中間表現（IR）基底のスパースサンプリング点（典型的には
50〜100点、低温ほど圧縮率が向上）へ置き換えます。オプションの
`sparse-ir <https://sparse-ir.readthedocs.io>`_ パッケージが必要です
（``pip install sparse-ir``）。カーネル・固有値反復・パリティ選別はすべて
スパースノード上で実行され、メモリと周波数軸の計算量が
:math:`N_{\mathrm{mat}}/L` 倍（20〜40×）削減されます。なお ``Nmat`` の
役割は従来どおりです：前段の FLEX 計算は引き続き一様 ``Nmat`` グリッドで
出力し（その収束は必要）、下記の出力もそのグリッドへ密評価されます。IR が
圧縮するのは動的ソルバーの**内部**周波数軸で、そのノード数は ``Nmat`` では
なく :math:`\beta\,\omega_{\max}` と ``ir_tol`` で決まります。GPU 実行
（``gpu = true``、CuPy が必要）とも併用でき、``fft_workers`` は従来どおり
CPU の空間 FFT にのみ作用します。

出力（``gap_dynamic.npz`` / ``gap.dat``）は従来どおり一様グリッドへ密評価
して書き出されるため、下流の解析はそのまま動作します（npz には
``matsubara_basis`` などの由来メタデータが追加されます）。

.. note::

   一様グリッドFFTで計算された FLEX 出力（``chiq_s.npz`` 等）には
   :math:`O(\beta/N_{\mathrm{mat}})` の離散化アーティファクト
   （:math:`\delta(\tau)` 由来の定数オフセットとエイリアシング・イメージ）
   が含まれます。IR 読み込みは定数成分を分離・除去し（ログに記録）、
   一様パスとの固有値の差はこの入力データ品質程度（小規模テスト
   フィクスチャの実測で ``Nmat=128`` のとき約1%、``Nmat=512`` で 5×10⁻⁴。
   これらはフィクスチャ固有の値であり一般保証ではありません）になります。
   両者は ``Nmat`` を増やすと同じ連続極限に収束します。本番計算では
   モデルごとに一度、FLEX の ``Nmat`` を上げる（または同じ ``Nmat`` で
   uniform と IR を比較する）ことで先頭固有値の変化が許容範囲か確認して
   ください。

.. _sc_dynamic_gpu:

GPU実行（CuPy）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``[eliashberg]`` セクションで ``gpu = true`` と設定すると、動的モードの
カーネル適用（固有値ソルバーの matvec）を GPU 上で実行します。2つの大きな
不変テンソル（ペアバブル :math:`[GG]` とペアリング頂点）を反復開始前に一度だけ
GPU へ転送して常駐させ、各反復ではギャップベクトルのみを転送します。結果は
CPU 実行と数値的に同一（倍精度の丸め誤差の範囲内）です。

- ``frequency = "dynamic"`` のみに適用されます。静的ソルバーは CPU 専用のため、
  ``frequency = "static"``\ （または省略時。既定は静的）で ``gpu = true`` を設定
  すると、フラグを黙って無視せず ``ValueError`` で即時停止します。
- `CuPy <https://cupy.dev/>`_ がインストールされ CUDA デバイスが利用可能で
  ある必要があります。CuPy が無い・デバイスが見つからない場合は、警告を出して
  自動的に CPU（numpy）実行へフォールバックします（結果は同一で、速度のみ
  低下します）。
- 必要な GPU メモリはおおよそ常駐テンソル 2 本分
  :math:`2 \times 16\, N_{\mathrm{orb}}^4\, N_k\, N_{\mathrm{mat}}` バイト＋
  作業領域です。不足する場合は CuPy が明示的な OutOfMemory エラーで停止します。
- 参考値: :math:`N_{\mathrm{orb}}=2`, :math:`64\times 64` k メッシュ,
  :math:`N_{\mathrm{mat}}=2048` で matvec あたり CPU 比 16 倍程度
  （NVIDIA RTX 6000 Ada、GPU メモリ約 5 GB 使用）。

対応する相互作用
----------------------------

Eliashberg方程式ソルバーは、H-waveで利用可能な
全ての相互作用型に対応しています:

- ``CoulombIntra`` (:math:`U`): 軌道内クーロン斥力
- ``CoulombInter`` (:math:`V`): 軌道間クーロン斥力
- ``Hund`` (:math:`J`): フント結合
- ``Exchange`` (:math:`J'`): ペアホッピング（交換）
- ``Ising`` (:math:`I`): イジング型スピン相互作用
- ``PairHop`` (:math:`P`): ペアホッピング

``Hund``、``Exchange``、``Ising``、``PairHop``
相互作用が存在する場合、ソルバーは自動的に
一般化 :math:`S`/:math:`C` 行列定式化
（Kuroki et al., PRB 79, 224511）を使用し、
4インデックスの頂点構造で計算します。


Tips
----------------------------

- 大規模系では ``chi0q_mode = "calc"`` と設定すると、
  :math:`\chi_0(\mathbf{q})` を内部計算し、
  大きなファイルの読み込みを回避できます。
- ``"arnoldi"`` 固有値法は少数の主要固有値を求めるのに最速です。
  縮退した固有値がある場合は ``"subspace"`` がより堅牢です。
- 反復法では異なる ``init_gap`` 対称性を使用して、
  特定のペアリングチャネルを狙うことができます。
  固有値法は全ての主要対称性を自動的に見つけます。
- ``pairing_type = "triplet"`` オプションで、
  適切な頂点を用いた三重項ペアリング不安定性を解析できます。
