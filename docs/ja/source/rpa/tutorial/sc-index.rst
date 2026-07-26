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

  .. note::

     ``"calc"`` は **この** 入力ファイルの ``mode.param`` の設定で
     :math:`\chi_0(\mathbf{q})` を再計算します。``"load"`` の結果と比較可能に
     するには、chi0q ファイルを生成した計算とこれらの設定を一致させる必要が
     あります。特に ``coeff_tail`` （Matsubara 高周波数裾補正）は中程度の
     ``Nmat`` において :math:`\chi_0(\mathbf{q})` を :math:`O(1)` で変えます。
     デフォルトの ``coeff_tail = 0.0`` （補正なし）は収束が悪い場合があり、
     実用的な RPA/FLEX 計算では通常 ``coeff_tail = 1.0`` を用います。
     chi0q ファイルには生成時の値が記録されており、読み込み時に設定と
     食い違う場合は ``hwave_sc`` が警告を出します。
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
  ``"d_y2z2"`` (:math:`\cos k_y - \cos k_z`) 、
  ``"d_xy"`` 、 ``"d_xz"`` 、 ``"d_yz"`` 、 ``"d_z2"`` 、
  ``"p_x"`` 、 ``"p_y"`` 、 ``"p_z"`` 、
  ``"f_x"`` (:math:`\sin k_x (\cos k_x - \cos k_y)`) 、
  ``"f_y"`` (:math:`\sin k_y (\cos k_y - \cos k_x)`) 、 ``"random"`` です。
  2つの :math:`f` 波的な初期値は奇パリティ（三重項パリティ）の高次調和関数で、
  正方格子では ``"p_x"``/``"p_y"`` と合わせて奇パリティの2次元表現を張るため、
  点群による分離ではなく初期値としての位置づけです。
  準2次元セル ``CellShape = [1, Ny, Nz]`` （:math:`k_x = 0`）では、
  :math:`\sin k_x` を含む形状因子（``"p_x"`` 、 ``"d_xy"`` 、 ``"d_xz"``）は
  恒等的に消えるため無効です。三重項には ``"p_y"``/``"p_z"`` を、面内 :math:`d` 波
  には :math:`(\pi,0)` と :math:`(0,\pi)` でノードを持つ ``"d_yz"`` ではなく、
  同点で逆符号のアンチノードを持つ ``"d_y2z2"`` を用いてください。
- ``max_iter``: 自己無撞着反復の最大回数。
- ``alpha``: 混合パラメータ（0: 混合なし、1: 古い解を完全保持）。
- ``convergence_tol``: ギャップ関数の収束条件。
- ``num_eigenvalues``: 固有値モードで計算する固有値の数。
- ``eigenvalue_method``: ``"arnoldi"`` （デフォルト）、 ``"subspace"`` 、
  ``"shift-invert-gmres"`` / ``"shift-invert-bicgstab"`` /
  ``"shift-invert-lgmres"``。
- ``sigma_shift`` （shift-invert 系の ``eigenvalue_method`` のみ）: shift-invert
  ソルバの実数ターゲット :math:`\sigma` 。 :math:`\sigma` 近傍の固有値が先に
  求まります。素の ``"arnoldi"`` では無視されます（警告あり）。arnoldi では
  代わりに ``spectral_shift`` を使ってください。
- ``spectral_shift`` （ ``eigenvalue_method = "arnoldi"`` および
  ``solver_mode = "iteration"`` で有効）:
  正の数値または ``"auto"`` 。ARPACK の既定の選択 （ ``which='LM'`` ）は
  *絶対値* 最大の固有値を返すため、超伝導不安定性から遠い場合、小さな正
  （引力的）の主固有値が、より大きな負（斥力的）の固有値に隠れて取りこぼされ、
  報告される主固有値が負（非物理的）になることがあります。 ``spectral_shift``
  を指定すると、シフトした演算子 :math:`A + \sigma I` に対して *実部* 最大の
  固有値（ ``which='LR'`` ； :math:`T_c` で :math:`\lambda \to 1` となる物理的な
  SC固有値）を求めます。シフトは内部で差し引かれるので、受け取る／保存される
  固有値はシフト前の正しい値です。 ``"auto"`` はスペクトル半径から
  :math:`\sigma` を自動設定します。明示的に指定する場合は、最も負の固有値の
  *絶対値* より大きい正の :math:`\sigma` を与えます（ :math:`A + \sigma I` の
  スペクトルが全て正の実部になるように）。主固有値が負になる場合や、対形成が
  弱い系（低圧・擬1次元）を走査する場合に推奨します。上記の ``sigma_shift``
  （shift-invert のターゲット）とは別物である点に注意してください。

  ``spectral_shift`` は **べき乗反復** （ ``solver_mode = "iteration"`` 、
  ``"both"`` および動的ソルバの反復パス）にも適用されます。斥力優勢な
  カーネル（ボンド分解パス ``bond_channels = true`` では通常こちら）では
  主固有値が **負** となるため、反復ごとにベクトルの符号が反転し、収束判定量
  :math:`\|\sigma_{n+1}/\|\sigma_{n+1}\| - \sigma_n\|` が :math:`\approx 2`
  に張り付いて、 ``max_iter`` をいくら大きくしても収束しません。シフトを
  指定すると、主固有値が正となるシフト演算子 :math:`K + \sigma I` で反復し、
  最後に :math:`\sigma` を差し引くため、 ``eigenvalue.dat`` に書かれる値は
  元のカーネルの **符号付き** 固有値 :math:`\lambda` （負でも構いません）に
  なります。固有ベクトル（ギャップ）はシフトの影響を受けません。

  ただし、この同一視が成り立つのは、反復ベクトルが :math:`K` の固有ベクトルで
  あり、かつその *シフト後* の固有値が正の実数である場合に限られます。明示的な
  （小さすぎるかもしれない） :math:`\sigma` 、 ``"auto"`` の推定値、
  ``max_iter`` で打ち切られた実行のいずれもこれを保証しません。そのため
  ソルバは結果を検証します。すなわち **シフトしない** カーネルに対する
  符号付きレイリー商
  :math:`\lambda = \langle v|K|v\rangle / \langle v|v\rangle` （行列ベクトル積
  1回分）と残差 :math:`\|Kv-\lambda v\|/\|v\|` を計算し、残差が許容範囲内で
  あればそのレイリー商を固有値として報告し、 ``eigenvalue.dat`` に
  ``VALIDATED`` と明記します。許容範囲外の場合は、その数値が
  :math:`K` の固有値では **なく** シフト反復ノルムの推定値
  :math:`\|(K+\sigma I)v\| - \sigma` である旨をファイルに明記し、警告を
  出します（ ``max_iter`` を増やす、 ``spectral_shift`` を大きくする、または
  ``solver_mode = "eigenvalue"`` を使ってください）。シフトが小さすぎた場合
  （選ばれたモードが :math:`\lambda + \sigma < 0` 、すなわち実部最大ではなく
  絶対値最大のモードに収束した場合）も、ログと出力の注記の双方で個別に
  警告されます。

  arnoldi と
  同じく、選ばれるのは（シードとパリティ射影が定めるセクター内で） *実部最大*
  の固有値、すなわち物理的なSC固有値であり、絶対値最大の斥力モードでは
  ありません。 ``"auto"`` は、初期ギャップと決定論的な乱数ベクトルを
  プローブとして（反復と同じパリティ射影を通して）数回のべき乗反復で
  スペクトル半径 :math:`\rho` を見積もり、
  :math:`\sigma = 1.05\,\rho_\mathrm{est} + 10^{-6}` 、すなわち推定半径のすぐ上
  に設定します。 :math:`\sigma` は条件を満たす範囲でできるだけ小さく取って
  ください。シフトが大きすぎると固有値の *相対* 間隔が詰まり、べき乗反復の
  収束が遅くなります（その場合は ``max_iter`` を増やしてください）。
  採用された :math:`\sigma` と「報告される :math:`\lambda` はシフトを
  戻した値である」旨は INFO レベルでログ出力されます。
- ``bond_channels``: ``true`` / ``false`` （デフォルト ``false`` ）。ボンド分解
  相互作用チャネルを有効にします。サイト間相互作用 :math:`V_{ab}(R)` を
  :math:`\Delta r = 0` （密度）チャネル1本に潰さず、拡大された
  :math:`(Z+1)\times(Z+1)` のボンドチャネル行列として保持するため、
  :math:`V` の Fock（交換）項と裸の対（Cooper）項の構造がペアリング頂点に
  入ります。 ``CoulombInter`` 項が必須です（値が 0 でも宣言されていれば可。
  :math:`V` スキャン中もチャネル構成が変わりません）。本リリースでの制限:
  1軌道（ ``norb = 1`` ）のみ、実かつ反転対称な ``CoulombInter`` のみ、
  ``frequency = "static"`` のみ、 ``chi0q_mode = "flex"`` とは併用不可
  （このモードが読み込むのは FLEX の *chi* であり、ボンド構造を持たないため。
  FLEX の *Green 関数* を使いたい場合は下の ``bond_green`` を参照）。
  いずれも明示的なエラーになります。このパスは Green 関数から自前のボンド分解
  バブルを構成するため、 ``chi0q_mode`` / ``chi0q_tensor`` は使われません
  （chi0q ファイルは読みません）。静的 RPA ラダーによる dressing であり
  保存則を満たす FLEX の結果ではないので、 :math:`\lambda` の絶対値を
  FLEX／動的計算と直接比較することはできません。固有値ファイルの ``#``
  コメント行に近似レベルとチャネル一覧が記録されます。
- ``bond_diagnostics``: ``true`` / ``false`` （デフォルト ``false`` 。
  ``bond_channels = true`` のときのみ有効）。主要状態の **キャラクター解析**
  をオプトインで出力します。出力は固有値ファイルへの ``#`` コメント行の
  追加のみで、数値行は一切変わりません（既存のリーダーに影響なし。フラグを
  切っていれば出力はこれまでとバイト単位で同一です）。記録される内容:

  * ``bond_diagnostics_harmonics`` — 得られた（主要）ギャップを奇パリティの
    調和関数基底
    :math:`\{\sin k_x,\ \sin k_y,\ \sin k_x(\cos k_x - \cos k_y),\
    \sin k_y(\cos k_y - \cos k_x)\}` へ、 :math:`\sqrt{GG}` 対重み計量で
    分解した結果。各値はその調和関数が状態に捕捉される割合（ :math:`[0, 1]` ）
    です。 :math:`V` スキャンの解釈はこの量に依ります。追跡している状態の
    f 成分はボンド機構の物理的な観測量であり、1軌道 Onari モデルでは
    :math:`V = 0` の約 21 % から :math:`V = 1.2` の約 68 % まで増大します。
    なお正方格子では ``p_x`` / ``p_y`` と ``f_x`` / ``f_y`` は同一の奇2次元
    表現に属するため、これは点群による分離ではなくキャラクターの報告です。
  * ``bond_diagnostics_eigenvalue_clusters`` — 計算された固有値の準縮退
    クラスタリング（ :math:`|\lambda_i - \lambda_j| \le` ``deg_tol`` の
    単連結法。既定値 :math:`10^{-3}` も併せて記録）。下の数値行へのインデックス
    のリストで、 :math:`\mathrm{Re}\,\lambda` の降順に並びます。併せて出力される
    ``bond_diagnostics_leading_cluster`` （およびそのサイズ）は **0行目** を
    含むクラスタ、すなわちチャネルパリティによる並べ替えの後に実際に返される
    固有対のクラスタです（ :math:`\lambda` 最大のクラスタとは限りません）。
    報告された状態が正方格子の奇 :math:`E` **2重項** なのか非縮退の枝なのかを
    直接確認できます。実際にスペクトルを計算した場合、すなわち ``solver_mode``
    に ``"eigenvalue"`` を含む場合のみ出力されます（反復パスはベクトル1本のみで
    スペクトルを持たないため、調和関数分解のみを報告します）。
  * 既存の ``lambda_pp`` / ``lambda_fl`` の寄与分解（フラグに関係なく
    全ボンド計算で出力）も同じブロックに並びます。

  スキャンを **またいで** 1つの不変部分空間を追跡する操作は多点計算なので、
  ソルバではなくライブラリ側にあります。スキャンドライバから
  ``hwave.solver.bond_channels.track_subspace`` を各点の固有対
  （固有ベクトルは1本を1 **行** として渡します。
  ``scipy.sparse.linalg.eigs`` は列で返す点に注意）とともに呼んでください。
- ``bond_green`` （ ``bond_channels = true`` のときのみ）: 外部から与える
  Green 関数 ``npz`` へのパス（ H-wave の ``green`` 配列、形状
  ``(nblock, nfreq, nvol, norb, norb)`` = FLEX/UHF が出力する形式）。
  トランスファーから作る裸の Green 関数の **代わりに** これを使います。
  FLEX で自己無撞着に求めた Green 関数をボンドパスに渡す手段であり、
  Onari のマイルストーンの値はこの方法でのみ再現できます（当該パラメータでは
  裸の RPA スピンチャネルは既に不安定領域に入っているため）。
  ``chi0q_mode = "calc"`` / ``"load"`` と併用可能です（どちらもこのパスでは
  使われません）。 ``chi0q_mode = "flex"`` は FLEX の *chi* を読み込む
  モードなので、引き続きエラーです。ファイルの振動数点数は ``Nmat`` より
  優先されます（ボンドバブルはその格子上で構成されるため）。不一致は警告し、
  ファイルが無い・ ``green`` 配列が無い・格子や軌道数が合わない場合は
  エラーになります。実際に使った Green 関数のパスは固有値ファイルに
  ``bond_green`` として、対応する ``approximation`` 記述とともに記録されます。
- ``bond_max_shells`` （ ``bond_channels = true`` のときのみ）: 近接シェル
  ``0..n`` を残します（シェル0 = オンサイト :math:`\Delta r = 0` 、
  シェル1 = 最近接、…）。省略時は宣言された全シェル。同一のフィルタ済み
  相互作用が Hartree・Fock・Cooper のすべてに使われるため、切り詰めの
  不整合は起きません。非ゼロのサイト間 :math:`V` があるのに ``0`` を
  指定した場合はエラーです（本当に局所的な模型なら
  ``bond_channels = false`` を使ってください）。
- ``bond_memory_cap_gb`` （ ``bond_channels = true`` のときのみ）: ボンドパスの
  ピークメモリ上限（GB、デフォルト ``8.0`` ）。事前見積りで
  :math:`N_q \times ND \times ND` 配列（ :math:`ND = n_{\rm orb}^2 B` ）の
  サイズを評価し、超える場合は原因となるチャネル数を示してエラーにします
  （暴走した確保は行いません）。
- ``bond_precondition_atol`` / ``bond_precondition_rtol`` /
  ``bond_precondition_dense_limit`` （ ``bond_channels = true`` のときのみ）:
  実行時エルミート性チェックの許容値と dense／probe 切り替え
  （デフォルト ``1e-8`` / ``1e-8`` / ``1024`` ）。残差は :math:`W K` 上で
  測りますが、判定では残差とその基準スケールの **両方** を厳密な変換
  :math:`\tilde K = W^{-1/2} (W K) W^{-1/2}` で対称化カーネルに換算します。
  したがって ``bond_precondition_rtol`` は :math:`\tilde K` の *相対* 非対称度
  であり、 :math:`\beta` を上げて :math:`\min w` が小さくなっても判定基準が
  勝手に厳しくなることはありません。一方、 :math:`w` の小さい成分に局在した
  破れは（ :math:`W K` 上では :math:`w_i w_j` で抑えられますが
  :math:`\tilde K` 上では :math:`\sqrt{w_i w_j}` にとどまるため）従来どおり
  検出されます。両方の値が固有値ファイルに
  記録されます（ ``kernel_hermiticity_residual`` と ``..._residual_ktilde`` ）。
  ``bond_precondition_dense_limit = 0`` にすると、厳密な dense 残差が高価な
  格子でも確率的 probe 推定を強制できます。

  **dense チェックのコスト。** ``bond_precondition_dense_limit`` 以下では、
  完全なカーネルを dense に構築（``vec_size`` 回のボンドカーネル matvec）する
  だけでなく、``pp``／``fl`` の各パートについても別途 dense に構築する
  （それぞれ ``vec_size`` 回）ため、固有値ソルバーが始まる前に合計
  :math:`3 \times` ``vec_size`` 回のボンドカーネル matvec が実行されます。
  デフォルトの ``bond_precondition_dense_limit = 1024`` では、32×32 の単一
  バンド格子（``vec_size`` :math:`= n_{\rm orb}^2 N_x N_y N_z = 1024`）が
  ちょうどこの上限に一致するため、このプリコンディション チェックだけで
  約 3072 回のボンドカーネル matvec がかかります。各 matvec は安価な dense
  行列ベクトル積ではなく、:math:`O(B)` の FFT 分解演算です。これは（SCF や
  固有値ソルバーの反復ごとではなく）1 回の求解あたりの前払いコストですが、
  この格子サイズでは無視できません。気になる場合は
  ``bond_precondition_dense_limit`` を下げて、より安価な確率的 probe 推定に
  切り替えてください。

  **不安定性ガード。** 拡大 RPA のドレッシングを行う前に、ボンドパスは各
  :math:`q` について 2 つの解行列 :math:`I - \bar\chi S` と
  :math:`I + \bar\chi C` の条件数を、相対基準
  :math:`\sigma_\mathrm{min}/\sigma_\mathrm{max}` （下限 :math:`10^{-3}` ）で
  検査します。磁気・電荷不安定性領域に入った計算は、チャネル（spin か
  charge）、 :math:`q` 点、比の値を明示したエラーで停止します。厳密に特異な
  場合の素の ``LinAlgError`` や、ほぼ特異な場合に巨大で無意味な頂点関数を
  黙って返す挙動を防ぐためです。この検査はボンドパス専用で、
  ``bond_channels = false`` には影響しません。

  ``bond_channels = false`` のときは、これらボンドオプションはすべて
  警告付きで無視されます。

.. note::

   ``bond_channels = true`` では ``solver_mode = "eigenvalue"`` を使うか、
   ``solver_mode = "iteration"`` に ``spectral_shift`` を **併用** して
   ください。シフトなしのべき乗反復が返すのは :math:`\|A\sigma\|` という
   **符号を持たない** ノルムであり、斥力優勢なボンドカーネルでは最大固有値が
   **負** になるため、反復ベクトルが毎ステップ符号反転して ``max_iter`` を
   いくら増やしても収束しません。このパスでは警告を出し、固有値ファイルにも
   「未収束の符号なしノルム」であることを明記します。符号を持つ物理量は
   provenance の ``lambda_rayleigh`` です。 ``spectral_shift = "auto"`` を
   加えるとこの問題は解消し、 :math:`K + \sigma I` で反復して収束し、元の
   カーネルの符号付き :math:`\lambda` （ ``solver_mode = "eigenvalue"`` や
   ``lambda_rayleigh`` と一致）が報告されます。ただしこれはシフト実行の
   レイリー商による検証を通過した場合に限られます。通過しなかった場合は
   固有値ファイルにその旨が明記され、その数値を :math:`\lambda` として
   引用してはいけません（上記 ``spectral_shift`` 参照）。
- ``gpu``: ``true`` でカーネル適用（matvec/matmat、FFT 畳み込み）を GPU（CuPy）で
  実行します。``frequency = "dynamic"`` と ``frequency = "static"`` の **両方** に
  対応します（デフォルト ``false``）。固有値ソルバー本体（ARPACK Arnoldi・べき乗
  反復・部分空間反復・shift-invert）は常にホストで実行されます（CuPy に一般の
  非エルミート固有値ソルバーが無いため）。反復ごとにデバイスへ渡るのはギャップ
  ベクトルのみです。使用可能な CuPy/CUDA デバイスが無い場合は警告を出して CPU に
  フォールバックします（下記の :ref:`GPU実行の節 <sc_dynamic_gpu>` を参照）。GPU
  実行時は ``fft_workers`` は無視されます。``bond_channels = true`` の場合、
  このフラグは無条件で無視されます: ボンド分解カーネルは v1 では CPU 専用
  なので、ソルバーは警告（"gpu=true is ignored for bond_channels=true"）を
  出して CPU で実行します —— エラーには **なりません**。
- ``gpu_required``: ``true`` にすると ``gpu = true`` を厳格化し、CuPy/CUDA が
  使えない場合に静かに CPU へフォールバックせずエラーで停止します
  （デフォルト ``false``）。動的 Eliashberg ソルバー（``[eliashberg]`` に設定）
  および FLEX・RPA ソルバー（``gpu`` フラグと同じく ``[mode.param]`` に設定）で
  有効です。
- ``fft_workers``: 動的モードの空間 FFT のワーカースレッド数
  （デフォルト ``1`` = 従来どおりの直列 numpy。``-1`` で全コア。
  GPU 実行時は無視されます）。
- ``matsubara_basis``: 動的モードの松原軸の表現。デフォルトは ``"uniform"``
  （従来どおり）。``"ir"`` を指定すると sparse-ir による中間表現基底を
  使用します（下記の :ref:`IR基底の節 <sc_dynamic_ir>` を参照）。
- ``ir_tol``: IR基底の打ち切り精度 :math:`\varepsilon` 。デフォルトは 1e-8。
- ``ir_wmax``: IR基底の実周波数バンド幅 :math:`\omega_{\max}` 。省略時は
  分散のスペクトル範囲 :math:`\max|\varepsilon_k-\mu|` と相互作用スケールから
  自動推定します（推定できない場合は明示指定を求めるエラーになります）。
- ``ir_keep_static_chi``: ``true`` / ``false`` 。デフォルトは ``false``。
  スピン・電荷感受率が静的支配的な場合（サンプリング窓内でほぼ振動数に依存せず
  大きい、すなわち臨界近傍）、IR圧縮が :math:`O(\beta/N_\mathrm{mat})` の
  :math:`\delta(\tau)` アーティファクトとして破棄する振動数非依存成分が
  物理的な重みを持ち、破棄すると最大固有値が誤ります。破棄成分がデータスケールを
  超える場合はソルバーが停止します。``true`` にすると停止せず静的成分を保持します
  （あるいは ``ir_wmax`` を下げる、FLEX の ``Nmat`` を増やす）。

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

ボンドチャネルの provenance（``bond_channels = true``）
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``[eliashberg] bond_channels = true`` のとき、``eigenvalue.dat`` の数値表の
上に追加の ``# key = value`` コメント行が付きます（純粋に追加のみで、
``np.loadtxt`` など既存の読み込み側は ``#`` 行を読み飛ばすため影響を受けません）。
記録される内容:

- ``bond_channels = True``、``bond_n_channels``（オンサイトを含む解決済み
  ボンドチャネル数 :math:`B`）、``bond_delta_r``（各チャネルの
  :math:`\Delta r_m` シェルベクトルのリスト、``bond_n_channels`` と
  同じ順序）。
- ``bond_max_shells``（設定した打ち切り、または ``"all"``）と
  ``bond_memory_cap_gb`` （設定した上限）: 再現性のために記録されるだけで、
  この時点ではもう計算結果に影響しません（シェルのフィルタリングと
  リソース事前見積りはオペレータ構築の前にすでに完了しています）。
- ``approximation`` — 実際に何を計算したかを表す自由記述: トランスファー
  から作った裸の Green 関数、または（``bond_green`` を設定した場合は）
  そこで指定した外部 Green 関数のどちらかの上で行った静的 RPA ラダー
  dressing であること、そしてどちらも保存則を満たす FLEX の結果ではない
  ため :math:`\lambda` の絶対値を近似レベルの異なる計算間で比較できない
  ことを明記します。
- ``bond_green`` — ``[eliashberg] bond_green`` を設定したときのみ出現し、
  実際に読み込んだ外部 Green 関数のパスを記録します。
- ``collapsed_to_pure_hubbard = True`` — ``bond_n_channels == 1``
  のとき（例えばオンサイト以外の :math:`V` を宣言せずに
  ``bond_max_shells = 0`` としたときなど、オンサイト以外のチャネルが
  1本も解決されなかった場合）のみ出現し、``bond_channels = true`` でも
  実質的には純粋 Hubbard（潰れたボンド）極限であることを示します。
- ``kernel_hermiticity_residual`` / ``kernel_hermiticity_residual_ktilde`` /
  ``kernel_hermiticity_relative_ktilde`` / ``kernel_hermiticity_method`` —
  実行時エルミート性事前条件チェックの診断値（上記の
  ``bond_precondition_atol`` / ``bond_precondition_rtol`` /
  ``bond_precondition_dense_limit`` を参照）: 生の残差と対称化カーネル
  :math:`\tilde K` の残差、その *相対* :math:`\tilde K` 残差、
  どちらのチェック（``"dense"`` または ``"randomized-probe"``）で
  得られたか。
- ``pair_weight_min_eigenvalue`` — 事前条件チェック中に見つかった
  :math:`\sqrt{GG}` 対重み計量の最小固有値。下記の Rayleigh 商による
  寄与分解が well-defined であるためには正である必要があります。

**「** :math:`\lambda = \lambda^{\rm pp} + \lambda^{\rm fl}` **」の寄与分解。**
これはボンドパスの目玉となる診断量であり、``bond_diagnostics`` の設定に
関係なく全てのボンド計算で出力されます:

- ``lambda_rayleigh`` — 完全な対称化カーネル :math:`\tilde K` の
  **符号付き** Rayleigh 商を、実行が返すギャップで評価した値
  （そのギャップが収束した固有ベクトルである場合に限りソルバー自身の
  固有値と一致します。下記の ``lambda_rayleigh_solver_mode`` /
  ``lambda_rayleigh_converged`` を参照）。
- ``lambda_pp`` — カーネルの **裸の particle-particle（瞬時）部分**
  の Rayleigh 商。ゆらぎによる dressing を受けない直接の Cooper 頂点
  :math:`\tfrac12 V^{\rm pp}` です。
- ``lambda_fl`` — **ゆらぎ媒介部分** の Rayleigh 商。ペアリング頂点への
  particle-hole ラダー（スピン・電荷感受率）の寄与です。
- ``lambda_attribution`` — :math:`\lambda = \lambda^{\rm pp} +
  \lambda^{\rm fl}` が :math:`\tilde K` の実 Rayleigh 商であり、両方とも
  *返された* ギャップで評価されていることを述べる自由記述です。
- ``lambda_rayleigh_solver_mode`` — Rayleigh 商を評価したギャップを
  どのソルバー（``"iteration"``、``"eigenvalue"``、``"both"``）が
  生成したか。
- ``lambda_rayleigh_converged`` — 反復パスでは、そのギャップが収束した
  かどうか（``"True"``/``"False"``）。固有値モードでは ``"n/a"`` です。

**寄与分解の読み方。** サイト間相互作用が斥力的（:math:`V > 0`）の場合、
裸の particle-particle 頂点自体がペアリングチャネルで斥力的であるため、
:math:`\lambda^{\rm pp} \le 0` になります。:math:`V` を増やすにつれて
:math:`\lambda` が上昇し（最終的には引力的にすらなり得る）場合、それは
**完全にゆらぎ駆動** の効果であり、:math:`\lambda^{\rm fl}` によって
運ばれています——ボンド機構は「:math:`V` が直接電子対を作る」のではなく
「:math:`V` が電子対を作るスピン・電荷ゆらぎを駆動する」というものです。
（例えば :math:`V` スキャンにわたって）``lambda_pp`` / ``lambda_fl`` を
別々に見ることで、和だけを見るのではなくこの2つの機構を区別できます。

**オプトインのキャラクター解析（** ``bond_diagnostics = true`` **）。**
有効にすると、さらに3グループのキーが出現します（各項目の物理の詳細は
上記の ``bond_diagnostics`` パラメータの説明を参照）:

- ``bond_diagnostics_harmonics`` と ``bond_diagnostics_harmonics_note`` —
  返されたギャップを
  :math:`\{\sin k_x,\ \sin k_y,\ \sin k_x(\cos k_x-\cos k_y),\
  \sin k_y(\cos k_y-\cos k_x)\}` 基底へ分解した結果で、それぞれ
  :math:`[0, 1]` の割合です。
- ``bond_diagnostics_eigenvalue_clusters``、
  ``bond_diagnostics_leading_cluster``、
  ``bond_diagnostics_leading_cluster_size``、``bond_diagnostics_deg_tol`` —
  計算された固有値スペクトルの準縮退グルーピングと、返された状態を含む
  グループのサイズ（``solver_mode`` に ``"eigenvalue"`` を含む場合のみ
  出現し、反復パスは調和関数分解のみを報告します）。
- ``bond_diagnostics_clusters_note`` — クラスタリングの規約
  （:math:`|\lambda_i-\lambda_j|\le` ``deg_tol`` による単連結法）を
  説明する自由記述です。


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

   このバージョンから、IRネイティブな感受率ファイル（``write_densified =
   false``）はこの ``chi_convention`` タグを必ず持つ必要があり、また
   ローダーはテンソルの形状から推定されるレイアウトと矛盾する
   ``chi_convention`` を検出すると、誤った pairing vertex を黙って構築する
   代わりに明示的なメッセージで即座に停止するようになりました。現行の
   H-wave が書き出すファイルは常に整合したタグを持つため、これは手編集または
   外部生成の感受率ファイルにのみ影響します。

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

静的ソルバーと同様に、動的モードも単に代数的に最大の固有値ではなく、
**指定したペアリングチャネルの** 先頭固有対を報告します。フェルミオンの反対称性により、
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

.. _sc_channel_decomposition:

チャネル分解（診断用）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ペアリング頂点は、スピンゆらぎ項・電荷ゆらぎ項・瞬間的な裸の項の和です。
一重項・三重項チャネルではそれぞれ

.. math::

   V^{\mathrm s} &= \tfrac{3}{2}\, \hat S\, \chi_{\mathrm s}\, \hat S
       - \tfrac{1}{2}\, \hat C\, \chi_{\mathrm c}\, \hat C
       + \tfrac{1}{2}(\hat S + \hat C), \\
   V^{\mathrm t} &= -\tfrac{1}{2}\, \hat S\, \chi_{\mathrm s}\, \hat S
       - \tfrac{1}{2}\, \hat C\, \chi_{\mathrm c}\, \hat C
       + \tfrac{1}{2}(\hat C - \hat S),

であり、いずれもスピン（:math:`\chi_{\mathrm s}`）・電荷
（:math:`\chi_{\mathrm c}`）感受率について線形です（係数と裸の項が異なるだけ）。
したがってこの分解はどちらの ``pairing_type`` にも適用できます。
``[eliashberg]`` セクションの2つの任意のブール値により、頂点を構成する前に
一方のチャネルを 0 にでき、ペアリング強度をスピン／電荷ゆらぎに切り分けられます。

- ``zero_chi_c``\ （既定 ``false``）：:math:`\chi_{\mathrm c}` を 0 にする。
  頂点はスピンゆらぎ項＋裸の項を保持（スピンチャネル）。
- ``zero_chi_s``\ （既定 ``false``）：:math:`\chi_{\mathrm s}` を 0 にする。
  頂点は電荷ゆらぎ項＋裸の項を保持（電荷チャネル）。

.. note::

   これらのフラグは動的頂点（``frequency = "dynamic"``）にのみ適用されます。
   静的ソルバー（``frequency = "static"``、既定）はカーネルを RPA/FLEX の
   静的感受率から構成し、カーネル構成時にこれらを使用しません。
   ``frequency = "static"`` で指定しても効果はなく、警告が出ます。

例えばスピンチャネルの診断は次のように指定します。

.. code-block:: toml

   [eliashberg]
   zero_chi_c = true

0 にするのは動的感受率関数だけであり、裸の項に含まれる静的相互作用行列
:math:`\hat S` と :math:`\hat C` は変更しません。
両方を ``true`` にすると瞬間的な裸の頂点のみが残ります。いずれも既定で off の
ため通常の計算には影響せず、チャネルを 0 にすると警告が出ます。これらは
**診断用** です：瞬間的な裸の項は常に保持され、線形化ギャップ固有値問題は
頂点について非線形なので、別々に 0 にした計算の固有値は\ **加法的ではありません**\
（一般に
:math:`\lambda_{\mathrm s} + \lambda_{\mathrm c} \neq \lambda_{\mathrm{full}}`\ ）。
有効なフラグは ``gap_dynamic.npz``, ``gap.dat``, ``eigenvalue.dat`` に記録
されます。通常の本計算では両方を off にしてください。

固有ベクトル継続（``seed_eigenvector``）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

通常、報告される先頭固有対は「代数的に最大」のものです。振動数依存（非エルミート）
カーネルでは、これは *exceptional point*（2つの実固有値が衝突して複素共役対に分裂
する点）近傍で脆弱になり、FLEX自己エネルギーが滑らかに変化していても、隣接温度で
「先頭」ブランチが不連続にジャンプすることがあります。1つの物理ブランチを追うには、
``[eliashberg] seed_eigenvector`` に近傍の計算（例：1つ上の温度）が出力した
``gap_dynamic.npz`` を指定します。その gap が ARPACK の開始ベクトルとして使われ、
**かつ** 最大に重なる固有対を「最大の固有値」でなく選択します。温度を段階的に下げ、
各計算に直前の ``gap_dynamic.npz`` を与えれば、同じ対称性（例：d波）を連続的に追跡
できます。種は計算と同じ ``CellShape`` と ``Nmat`` でなければならず（不一致は即エラー）、
continuation スイープでは ``Nmat`` を固定してください。IR経路では種の gap は IR ノード
へ自動で再フィットされます。``[eliashberg] sigma_shift`` は shift-invert の狙い値を明示
指定します（未指定なら予備 Arnoldi から推定）。``sigma_shift`` をブランチ近傍に置き
``seed_eigenvector`` と併用するのが、隠れた／複素化する固有値を解決する最も堅牢な方法です。

温度継続計算（``hwave_tsweep``）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

上記のように ``sigma_init`` と ``seed_eigenvector`` を温度掃引で手動でつなぐには、
温度ごとに ``hwave``/``hwave_sc`` を実行し直し、各ステップの出力を次のステップの
入力へ手作業で配線する必要があります。``hwave_tsweep`` コマンド（``hwave``・
``hwave_sc`` と同様にパッケージに同梱されています）はこれを自動化します。1つの
ベース TOML（単一の FLEX+Eliashberg 計算で使う ``[mode]``/``[mode.param]``/
``[file]``/``[eliashberg]`` の設定そのもの）に ``[continuation]`` セクションを
加えて渡すと、降順の温度ラダーに沿って FLEX（および、無効化していなければ
Eliashberg ソルバー）を実行します。各段では、直前の段の収束した自己エネルギーを
``sigma_init`` としてこの段の FLEX へ（ウォームスタート）、直前の段の動的ギャップを
``seed_eigenvector`` としてこの段の Eliashberg 計算へ（固有ベクトル継続）与えます。
ウォームスタートの連鎖全体を自動化することで、各温度点をコールドスタートして
（毎回異なる準安定解に落ち着く恐れを抱えながら）計算するのではなく、1つの物理的
ブランチを低温まで滑らかに追跡できます。

段の間で変化させるのは ``mode.param.T`` のみです。``CellShape``、``Nmat``、
その他の形状を決定するフィールドはラダー全体で固定されます。これにより、各段の
``sigma_init``／``seed_eigenvector`` ファイルが次の段と形状互換になります。

``[continuation]`` セクション
""""""""""""""""""""""""""""""

.. code-block:: toml

    [continuation]
      temperatures   = [0.02, 0.015, 0.01, 0.008, 0.006]  # 明示的なラダー
      # または、`temperatures` が無い場合に生成するラダー:
      #   T_start = 0.02
      #   T_stop  = 0.006
      #   num     = 5
      #   spacing = "linear"          # "linear"（デフォルト）または "log"
      output_dir     = "tsweep"       # デフォルト
      run_eliashberg = true           # デフォルト
      warm_start     = true           # デフォルト
      seed_gap       = true           # デフォルト
      resume         = false          # デフォルト（または --resume）
      summary_file   = "lambda_vs_T.dat"  # デフォルト

- ``temperatures``: 明示的な温度のリスト。与えた順に実行されます。存在する場合、
  ``T_start``/``T_stop``/``num`` より優先されます。
- ``T_start`` / ``T_stop`` / ``num`` / ``spacing``: ``temperatures`` が無い場合に
  ラダーを生成するために使われます。``T_start`` と ``T_stop`` の間に ``num`` 個の
  点を生成し、``spacing`` は ``"linear"``（デフォルト）または ``"log"`` です。
  ``temperatures`` と ``T_start``/``T_stop``/``num`` の三つ組のどちらも与えない場合は
  pre-flight エラーになります。
- ``output_dir``（デフォルト ``"tsweep"``）: 掃引全体の親ディレクトリ。温度 ``T`` の
  段 ``idx`` は ``<output_dir>/<idx>_T<T>/output/`` に出力されます（``idx`` は3桁
  ゼロ埋め、``T`` は ``%g`` 形式）。
- ``run_eliashberg``（デフォルト ``true``）: 各段で Eliashberg ソルバーも実行します。
  これにはベース TOML に ``[eliashberg]`` セクションが必要で、無い場合は pre-flight
  が欠けているセクション名を挙げてエラーになります。``false`` にすると
  ``sigma_init`` のみを連鎖させる FLEX 専用の掃引になります。
- ``warm_start``（デフォルト ``true``）: 各段の収束した自己エネルギーを次の段の
  ``sigma_init`` として連鎖させます。
- ``seed_gap``（デフォルト ``true``）: 各段のギャップを次の段の
  ``seed_eigenvector`` として連鎖させます。これは動的 Eliashberg ソルバー
  （``[eliashberg] frequency = "dynamic"``）でのみ有効です。``seed_eigenvector``
  自体が動的モード専用のため、静的ラダーでは効果を持ちません。
- ``summary_file``（デフォルト ``"lambda_vs_T.dat"``）: ``<output_dir>/<summary_file>``
  に書き出されるサマリー表のファイル名。

掃引の実行
""""""""""""""""""""""""""""""

.. code-block:: bash

    $ hwave_tsweep input.toml

実行を制御する3つのフラグがあります:

- ``--dry-run``: 温度ラダー、各段の出力ディレクトリ、配線される
  ``sigma_init``/``seed_eigenvector`` のパスを解決して表示するだけで、どちらの
  ソルバーも呼び出しません。長い掃引を実行する前に ``[continuation]`` の設定を
  検証するのに使います。
- ``--keep-going``: デフォルトでは、ソルバーがエラーを送出した段で掃引は停止します
  （壊れた段はそれ以降のすべての種を汚染するため。部分的なサマリーは書き出されます）。
  ``--keep-going`` を指定すると、代わりに次の段をコールドスタートし、それが成功すれば
  以降の段への種として再び使われます。これは *1プロセス内でのエラー継続* であって、
  プロセス自体が中断された後の再開ではありません。
- ``--resume``（または ``[continuation] resume = true``）: *ジョブレベルの再開*。
  resume 付きで掃引を再実行すると、``hwave_tsweep`` は既に完了して種として使える段の
  連続する先頭部分を読み飛ばし、最初の未完了段から――直前の有効な段の ``sigma`` と
  動的ギャップを種として――あたかも掃引が止まらなかったかのように再開します。
  ウォールクロック／スケジューラによる kill、クラッシュ、手動中断の後に使います。

  段が「完了」と見なされるのは、記録されたサマリー行が error でなく、**かつ**
  ディスク上の出力が実際に存在して解析可能な場合だけです（途中まで書かれた／壊れた
  ``eigenvalue.dat`` は検出され、その段とそれ以降の段は再計算されます）。resume は
  小さなマニフェスト（初回実行時に書かれる ``tsweep_manifest.json``。解決済みラダーと、
  形状・物理設定=``CellShape``/``SubShape``/``Nmat``/``filling``/``Ncond``/相互作用
  ファイル/``[eliashberg]`` frequency・pairing のフィンガープリントを記録）で保護されます。
  異なるラダーや設定に対して resume すると、非互換な結果を混ぜずに **即時停止** します。
  サマリーとマニフェストは各段の後にアトミックに書かれるため、中断でチェックポイントが
  途中で切れることはありません。``--resume`` 無しでの再実行は新規実行として既存の掃引を
  段ごとに上書きします（既存の掃引を検出すると警告を出します）。

3つは別物です。**ウォームスタート**（``warm_start``/``seed_gap``）は1回の実行内で
ある段の結果を *次の* 段の種に繋ぎ、**--keep-going** は1回の実行内で *段がエラーした後*
の挙動を決め、**--resume** は *実行全体* が再起動されたときの挙動を決めます。

サマリーファイル
""""""""""""""""""""""""""""""

各実行は ``<output_dir>/<summary_file>``（デフォルト ``tsweep/lambda_vs_T.dat``）
に段ごとに1行を書き出します:

.. code-block:: text

    # idx  T  status  error_stage  Re_lambda  Im_lambda  parity_match  flex_converged  flex_iter
    0 0.02   ok    none 0.845000 0.000000 1 1 18
    1 0.015  ok    none 0.902000 0.000000 1 1 22
    2 0.01   error flex nan      nan      -1 0 -1
    ...

``status`` は以下のいずれかです:

- ``ok`` -- FLEX が収束し、``run_eliashberg`` が有効な場合はこの段の
  ``eigenvalue.dat`` から先頭固有対が読み取れた。
- ``not_converged`` -- FLEX が ``EPS`` を満たさないまま ``IterationMax`` に
  達したが、使用可能な自己エネルギー（Eliashberg を実行した場合はギャップも）は
  書き出された。このような段も次の段の種として利用できます。
- ``error`` -- ソルバーが例外を送出した、または（``run_eliashberg`` の場合）
  ``eigenvalue.dat`` が存在しないか解析できなかった場合。``error_stage`` に
  どちらのソルバーで失敗したか（``flex`` または ``eliashberg``）が記録されます。
- ``dry`` -- ``--dry-run`` によって生成された行。ソルバーは呼び出されていません。

浮動小数点値が欠けている場合（Eliashberg を実行しなかった、または失敗した場合の
``Re_lambda``/``Im_lambda``）は ``nan`` と表示され、整数フィールドが欠けている場合
（``parity_match``、``flex_converged``、``flex_iter``）は ``-1`` と表示されます。
``error_stage`` は ``status = error`` でない限り ``none`` です。

設定例
""""""""""""""""""""""""""""""

.. code-block:: toml

    [mode]
      mode = "FLEX"

    [mode.param]
      T         = 0.02
      CellShape = [32, 32, 1]
      Nmat      = 512
      filling   = 0.75

    [file]
    [file.input]
      path_to_input = "."

    [file.input.interaction]
      path_to_input = "."
      Geometry      = "geom.dat"
      Transfer      = "transfer.dat"
      CoulombIntra  = "coulombintra.dat"
      CoulombInter  = "coulombinter.dat"

    [file.output]
      path_to_output = "output"

    [eliashberg]
      frequency     = "dynamic"
      chi0q_mode    = "flex"
      pairing_type  = "singlet"
      solver_mode   = "eigenvalue"   # hwave_tsweep の pre-flight で必須

    [continuation]
      T_start        = 0.02
      T_stop         = 0.005
      num            = 6
      spacing        = "log"
      run_eliashberg = true
      warm_start     = true
      seed_gap       = true

この例では :math:`T = 0.02` から :math:`T = 0.005` まで、対数間隔の6段を
降順に計算します。各段で FLEX と動的 Eliashberg を実行し、``sigma_init`` と
``seed_eigenvector`` の両方を連鎖させ、``tsweep/lambda_vs_T.dat`` --
:math:`\lambda(T)` の表を書き出します。先頭の物理固有値が1を横切る点から
:math:`T_c` を見積もることができます。

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

FLEX 側が ``[mode.param] matsubara_basis = "ir"`` かつ
``write_densified = false``（IRネイティブ出力、``.npz`` 内の
``frequency_grid = "sparse_ir_nodes"`` キーで識別可能）の場合、動的ソルバーは
``matsubara_basis = "ir"`` でそのファイルを直接読み込みます。ノード集合が
一致すれば保存値をそのまま使用（通常のケース）、異なる場合は本実行の基底へ
リフィットします（残差をログ出力）。温度は FLEX 実行と一致している必要が
あります（感受率は物理入力のため、不一致はエラーです）。
``matsubara_basis = "uniform"`` のままIRネイティブ入力を渡すと明示エラーに
なります。本ソルバーを ``"ir"`` に切り替えるか、FLEX を
``write_densified = true`` で再実行してください。

.. note::

   一様グリッドFFTで計算された FLEX 出力（``chiq_s.npz`` 等）には
   :math:`O(\beta/N_{\mathrm{mat}})` の離散化アーティファクト
   （:math:`\delta(\tau)` 由来の定数オフセットとエイリアシング・イメージ）
   が含まれます。IR 読み込みは定数成分が小さい場合に分離・除去します
   （ログに記録）。フィットされた定数がデータスケールに匹敵する場合、
   それは離散化アーティファクトではあり得ないため、黙って結果を壊す
   代わりにエラーで停止します（対処法はエラーメッセージにも表示:
   ``ir_wmax`` を自動推定に戻す/``3×(バンド幅+最大相互作用)`` 程度に
   設定する、FLEX の ``Nmat`` を増やす、真に静的な成分なら
   ``ir_keep_static_chi = true`` で保持する、または
   ``matsubara_basis = "uniform"`` に戻す）。一様パスとの固有値の差は
   この入力データ品質程度（分散ベースの自動 ``ir_wmax`` での小規模
   フィクスチャ実測: ``Nmat=128`` で約1.5×10⁻²、``Nmat=512`` で約
   4×10⁻³。フィクスチャ固有の値であり一般保証ではありません）に
   なります。両者は ``Nmat`` を増やすと同じ連続極限に収束します。
   本番計算ではモデルごとに一度、FLEX の ``Nmat`` を上げる（または同じ
   ``Nmat`` で uniform と IR を比較する）ことで先頭固有値の変化が
   許容範囲か確認してください。

.. warning::

   Issue #57 修正より前のバージョンで計算した動的IRの結果は、ペアリング
   頂点にν非依存成分を持つモデル——特にオフサイト ``CoulombInter``
   を含む系——では誤りです（オンサイト ``CoulombIntra`` のみの模型は
   裸頂点項が厳密に相殺するため影響なし）。該当する計算は再実行して
   ください。λの大きな変化は想定どおりです（変化分がバグであり、物理の
   変更ではありません）。また自動 ``ir_wmax`` 推定は分散ベースに変わり、
   現実的な多ホッピング模型では大幅に小さい（正しい）値になります。

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
- ``gpu_required = true``（デフォルト ``false``）にすると、この暗黙の CPU
  フォールバックをハードエラーに変えます。``gpu = true`` を指定したのに利用可能
  な CuPy/CUDA バックエンドが無い場合、静かに CPU 実行へ移らずエラーで停止する
  ため、短い GPU ジョブが長大な CPU ジョブに化けるのをスケジューラ段階で防げます。
  同じフラグは FLEX・RPA ソルバーでも有効です。
- FLEX・RPA の GPU 経路は、大きなデバイス確保の前に advisory な VRAM プリフライト
  を実行します。常駐テンソルの見積りが空きデバイスメモリを超える場合、ソルバー名と
  見積り／空き容量を示す警告を出します（実際の確保時には CuPy が明示的な
  OutOfMemory エラーで停止します）。
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
  これはボンド分解パス（``bond_channels = true``）には当てはまり
  **ません**: そのカーネルは斥力優勢なので、単純なべき乗反復はそもそも
  収束できません（上記の ``spectral_shift`` の注記を参照）。まず
  ``spectral_shift`` を設定してください。また、それを設定していても
  ``init_gap`` のシードが依然として反復が見つける枝／パリティセクターを
  決めます——固有値法がスペクトル全体を走査するのとは異なり、シードは
  自動的にスキャンされないためです。
- ``pairing_type = "triplet"`` オプションで、
  適切な頂点を用いた三重項ペアリング不安定性を解析できます。
