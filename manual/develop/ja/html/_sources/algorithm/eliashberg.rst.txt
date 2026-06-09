.. highlight:: none

.. _algorithm_eliashberg:

線形化Eliashberg方程式
================================

概要
*****************************

線形化Eliashberg方程式ソルバー（ ``hwave_sc`` ）は、
RPA の枠組みで線形化ギャップ方程式の固有値問題を解くことにより、
超伝導不安定性を解析します。
超伝導転移温度 :math:`T_c` は、
最大固有値が :math:`\lambda = 1` に達する条件から決定されます。

アルゴリズムは以下のステップで構成されます:

1. 裸感受率 :math:`\hat{X}^{(0)}(\mathbf{q})` を計算または読み込む。
2. 非相互作用グリーン関数 :math:`G(\mathbf{k}, i\omega_n)` を構成する。
3. RPAペアリング頂点 :math:`V(\mathbf{q})` を構築する。
4. 線形化Eliashberg方程式を解いて最大固有値を求める。


グリーン関数
*****************************

k空間における非相互作用一体ハミルトニアンは、
ホッピング積分のフーリエ変換により得られます:

.. math::

   \varepsilon_{\alpha\beta}(\mathbf{k})
   = \sum_{\mathbf{R}} t^{\alpha\beta}_{\mathbf{R}}\,
     e^{i\mathbf{k}\cdot\mathbf{R}}

各k点でハミルトニアンを対角化し:

.. math::

   \sum_\beta \varepsilon_{\alpha\beta}(\mathbf{k})\, u_{\beta m}(\mathbf{k})
   = \xi_m(\mathbf{k})\, u_{\alpha m}(\mathbf{k})

非相互作用グリーン関数は以下のように与えられます:

.. math::

   G_{\alpha\beta}(\mathbf{k}, i\omega_n)
   = \sum_{m} \frac{u_{\alpha m}(\mathbf{k})\, u^*_{\beta m}(\mathbf{k})}
                    {i\omega_n - (\xi_m(\mathbf{k}) - \mu)}

ここで :math:`\omega_n = \pi(2n+1)/\beta` はフェルミオン松原振動数、
:math:`\mu` はフィリング条件から決定される化学ポテンシャル、
:math:`\beta = 1/T` は逆温度です。


ペアリング頂点
*****************************

簡易モード
-----------------------------

``CoulombIntra`` （ :math:`U` ）と ``CoulombInter`` （ :math:`V` ）のみが
存在する場合、ペアリング頂点はスピンチャネル（ :math:`W_s` ）と
電荷チャネル（ :math:`W_c` ）を用いて計算されます:

.. math::

   W_s = -U, \qquad W_c = U + 2V

RPA感受率は

.. math::

   \hat{X}^s(\mathbf{q}) = \left[\hat{I} - \hat{X}^{(0)}(\mathbf{q})\, \hat{W}_s\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

.. math::

   \hat{X}^c(\mathbf{q}) = \left[\hat{I} + \hat{X}^{(0)}(\mathbf{q})\, \hat{W}_c\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

一重項ペアリング頂点は

.. math::

   V^S_{\alpha\beta}(\mathbf{q})
   = \frac{1}{2}(W_c + W_s)_{\alpha\beta}
     + \frac{3}{2} (W_s\, X^s\, W_s)_{\alpha\beta}
     - \frac{1}{2} (W_c\, X^c\, W_c)_{\alpha\beta}

三重項ペアリング頂点は

.. math::

   V^T_{\alpha\beta}(\mathbf{q})
   = \frac{1}{2}(W_c - W_s)_{\alpha\beta}
     - \frac{1}{2} (W_s\, X^s\, W_s)_{\alpha\beta}
     - \frac{1}{2} (W_c\, X^c\, W_c)_{\alpha\beta}


一般モード（S/C行列定式化）
-----------------------------------------

``Hund`` （ :math:`J` ）、 ``Exchange`` （ :math:`J'` ）、
``Ising`` （ :math:`I` ）、または ``PairHop`` （ :math:`P` ）相互作用が
存在する場合、ソルバーは一般化 :math:`S` / :math:`C` 行列
定式化 [1]_ を使用します。

:math:`S` 行列と :math:`C` 行列は、軌道インデックス :math:`l_1, l_2` の
複合インデックス空間で定義されます。
:math:`n_{\rm orb}` 軌道の系では、行列の次元は
:math:`n_{\rm orb}^2 \times n_{\rm orb}^2` です。
行列要素は以下の通りです:

.. list-table::
   :header-rows: 1
   :widths: 30 20 25 25

   * - インデックス条件
     - 種類
     - :math:`S` の値
     - :math:`C` の値
   * - :math:`l_1 = l_2 = l_3 = l_4`
     - 軌道内
     - :math:`U`
     - :math:`U`
   * - :math:`l_1 = l_3 \neq l_2 = l_4`
     - クロス
     - :math:`U' - I`
     - :math:`-U' + J - I`
   * - :math:`l_1 = l_2 \neq l_3 = l_4`
     - 密度
     - :math:`J - 2I`
     - :math:`2U' - J`
   * - :math:`l_1 = l_4 \neq l_2 = l_3`
     - 交換
     - :math:`J' + P`
     - :math:`J' + P`

RPA感受率は

.. math::

   \hat{X}^s(\mathbf{q}) = \left[\hat{I} - \hat{X}^{(0)}(\mathbf{q})\, \hat{S}\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

.. math::

   \hat{X}^c(\mathbf{q}) = \left[\hat{I} + \hat{X}^{(0)}(\mathbf{q})\, \hat{C}\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

一重項ペアリング頂点は

.. math::

   \hat{V}^S(\mathbf{q})
   = \frac{3}{2}\, \hat{S}\, \hat{X}^s(\mathbf{q})\, \hat{S}
     - \frac{1}{2}\, \hat{C}\, \hat{X}^c(\mathbf{q})\, \hat{C}
     + \frac{1}{2}(\hat{S} + \hat{C})

三重項ペアリング頂点は

.. math::

   \hat{V}^T(\mathbf{q})
   = -\frac{1}{2}\, \hat{S}\, \hat{X}^s(\mathbf{q})\, \hat{S}
     - \frac{1}{2}\, \hat{C}\, \hat{X}^c(\mathbf{q})\, \hat{C}
     + \frac{1}{2}(\hat{C} - \hat{S})


線形化Eliashberg方程式
***********************************

線形化Eliashberg方程式は固有値問題として定式化されます:

.. math::

   \lambda\, \Sigma_{\alpha\beta}(\mathbf{k})
   = -\frac{T}{N_L} \sum_{\mathbf{k}', n', \alpha', \beta'}
     V_{\alpha\alpha';\beta\beta'}(\mathbf{k} - \mathbf{k}')
     \, G_{\alpha\alpha'}(\mathbf{k}', i\omega_{n'})
     \, G_{\beta\beta'}(-\mathbf{k}', -i\omega_{n'})
     \, \Sigma_{\alpha'\beta'}(\mathbf{k}')

ここで :math:`\Sigma_{\alpha\beta}(\mathbf{k})` は異常自己エネルギー
（ギャップ関数）であり、右辺はEliashbergカーネル :math:`K[\Sigma]` を定義します。

超伝導不安定性は :math:`\lambda = 1` で生じます。
SC転移には正の固有値のみが物理的に関連します。


二粒子グリーン関数
-----------------------------

松原振動数の和を解析的に実行することで、
二粒子グリーン関数が得られます:

.. math::

   G^{(2)}_{\alpha\beta;\gamma\delta}(\mathbf{q})
   = \frac{T}{N_L} \sum_{\mathbf{k}, n}
     G_{\alpha\gamma}(\mathbf{k}, i\omega_n)\,
     G_{\beta\delta}(-\mathbf{k}+\mathbf{q}, -i\omega_n)

これにより、Eliashbergカーネルはk空間での畳み込みに帰着し、
高速フーリエ変換（FFT）を用いて効率的に計算できます。


FFTによるカーネル評価
-----------------------------

カーネル演算 :math:`\Sigma_{\rm new} = K[\Sigma_{\rm old}]` は、
ペアリング頂点 :math:`V(\mathbf{q})` と
:math:`G^{(2)}` ・ギャップ関数の積との畳み込みを含みます。
この畳み込みは以下のように効率的に計算されます:

1. :math:`V(\mathbf{q})` と :math:`G^{(2)}(\mathbf{q}) \cdot \Sigma(\mathbf{q})` の逆FFTにより実空間に変換。
2. 実空間での要素ごとの積。
3. FFTによりk空間に戻す。

これにより計算量が :math:`O(N_k^2)` から :math:`O(N_k \log N_k)` に削減されます。


数値解法
*****************************

自己無撞着べき乗反復法
---------------------------------

自己無撞着べき乗反復法は、最大の正の固有値を持つ固有モードに収束します:

.. math::

   \Sigma^{(i+1)} = K[\Sigma^{(i)}], \qquad
   \lambda^{(i)} = \|\Sigma^{(i+1)}\|, \qquad
   \Sigma^{(i+1)} \leftarrow \Sigma^{(i+1)} / \lambda^{(i)}

収束を安定化するために線形混合が適用されます:

.. math::

   \Sigma^{(i+1)} \leftarrow (1-\alpha)\, \Sigma^{(i+1)}_{\rm norm}
   + \alpha\, \Sigma^{(i)}

ここで :math:`\alpha` は混合パラメータです。
:math:`\|\Sigma^{(i+1)} - \Sigma^{(i)}\| < \epsilon`
のとき収束したと判定します。

初期ギャップ関数は様々な対称性
（s波、 :math:`d_{x^2-y^2}` 、 :math:`\cos k_x + \cos k_y` 、ランダムなど）
を設定でき、特定のペアリングチャネルを狙うことができます。


Arnoldi固有値解析
---------------------------------

Arnoldi法（ARPACKによる陰的再起動法）は、
Eliashbergカーネルを線形演算子として、
その主要固有値を求めます。
この方法はカーネル行列を明示的に構成することなく、
複数の固有値を同時に効率的に計算できます。


部分空間反復法
---------------------------------

部分空間反復法は複数のベクトルを同時に伝播します:

1. 全ベクトルにカーネルを適用: :math:`W = K \cdot V`
2. レイリー商を計算: :math:`H = V^T K V`
3. 小行列 :math:`H` を固有値分解
4. リッツベクトルで部分空間を更新
5. QR分解による再直交化

この方法は縮退した固有値に対してより堅牢です。


シフト逆反復法
---------------------------------

シフト逆反復変換 :math:`(K - \sigma I)^{-1}`
を用いて、目標値 :math:`\sigma` 近傍の固有値を求めます。
線形方程式系はBiCGSTAB、GMRES、またはLGMRESにより
反復的に解かれます。


.. [1] K. Kuroki, S. Onari, R. Arita, H. Usui, Y. Tanaka, H. Kontani,
   and H. Aoki, Phys. Rev. Lett. **101**, 087004 (2008);
   K. Kuroki, H. Usui, S. Onari, R. Arita, and H. Aoki,
   Phys. Rev. B **79**, 224511 (2009).
