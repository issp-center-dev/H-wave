.. highlight:: none

.. _algorithm_sec:

非制限Hartree-Fock法
==========================

Hartree-Fock近似
*************************

Unrestricted Hartree-Fock(UHF)近似では一体の演算子の揺らぎについて一次のみを取り入れ、二体項を一体項へと近似します。
一般的な二体相互作用については以下の近似を行うことに相当します。

.. math::
   \begin{aligned}
   c_{i}^{\dagger}c_{j}^{\dagger}c_{k}c_{l} 
   &\sim \langle c_{i}^{\dagger} c_l\rangle c_{j}^{\dagger} c_k   +  c_{i}^{\dagger} c_l \langle c_{j}^{\dagger} c_k\rangle - \langle c_{i}^{\dagger} c_k\rangle c_{j}^{\dagger} c_l -  c_{i}^{\dagger} c_k \langle c_{j}^{\dagger} c_l\rangle \nonumber\\
   &\qquad -(\langle c_{i}^{\dagger} c_l\rangle \langle c_{j}^{\dagger} c_k\rangle - \langle c_{i}^{\dagger} c_k\rangle \langle c_{j}^{\dagger} c_l\rangle)
   \end{aligned}

H-waveでは以下の形式で二体相互作用を定義しています。
   
.. math::
   \begin{aligned}
   \mathcal{H}_\text{InterAll} &= \sum_{ijkl\alpha\beta\gamma\delta} \sum_{\sigma_1 \sigma_2 \sigma_3 \sigma_4}  I_{ijkl\alpha\beta\gamma\delta} c^\dagger_{i\alpha\sigma_1} c_{j\beta\sigma_2} c^\dagger_{k\gamma\sigma_3} c_{l\delta\sigma_4} \nonumber\\
   &= \sum_{ijkl\alpha\beta\gamma\delta} \sum_{\sigma_1 \sigma_2 \sigma_3 \sigma_4}  I_{ijkl\alpha\beta\gamma\delta} (c^\dagger_{i\alpha\sigma_1} c^\dagger_{k\gamma\sigma_3} c_{l\gamma\sigma_4} c_{j\beta\sigma_2} + c^\dagger_{i\alpha\sigma_1} c_{l\delta\sigma_4}\delta_{j,k}\delta_{\beta,\gamma}\delta_{\sigma_2,\sigma_3})
   \end{aligned}

ここで、上式の第2項のように一体項が存在することに注意が必要です。
さて、一体項で与えられたハミルトニアンは一般的に以下のように書けます。

.. math::
   \begin{aligned}
   \mathcal{H}_\text{UHF} &= \sum_{ij} H_{ij} c^\dagger_{i} c_{j} = \hat{c}^\dagger H \hat{c}
   \end{aligned}

ここで、簡単化のため、 :math:`i\equiv(i, \alpha, \sigma_1), j\equiv(j, \beta, \sigma_2)` 、 :math:`H` は :math:`H_{ij}` を成分に持つ行列、 :math:`\hat{c}` は :math:`c_{i}` を成分にもつ行ベクトルを表します。
このとき、 :math:`H` はエルミート行列なので、 :math:`\hat{\xi}` を :math:`H` の固有値を対角成分に持つ行列、:math:`U` は各固有ベクトルに対応する行列として、:math:`H=U \hat{\xi} U^\dagger` のように変形できることから、:math:`\hat{d} = U^\dagger \hat{c}` とすると、

.. math::
   \begin{aligned}
   \mathcal{H}_\text{UHF} &= \hat{d}^\dagger \hat{\xi} \hat{d} =  \sum_{k} \xi_k d_k^\dagger d_k 
   \end{aligned}
   
のように書き直すことが出来ます。したがって、UHF近似の一体相互作用からくるエネルギーは

.. math::
   \begin{aligned}
   E_\text{UHF} = \langle \mathcal{H}_\text{UHF} \rangle = \sum_{k} \xi_k \langle d_k^\dagger d_k \rangle
   \end{aligned}

として求められます。実際の数値計算では、:math:`H` はUHF近似を通して一体グリーン関数 :math:`\langle c_{i}^\dagger c_{j}\rangle` に依存するため、一体グリーン関数を最初に初期値として与え、

.. math::
   \begin{aligned}
   \langle c_{i}^\dagger c_{j}\rangle = \sum_{l} U_{il}^* U_{jl} \langle d_l^\dagger d_l \rangle = \sum_{l} \frac{U_{il}^* U_{jl}}{1+\exp^{\beta(\xi_l -\mu)}}
   \end{aligned}

の関係から一体グリーン関数を更新し、一体グリーン関数が収束するまで計算を繰り返します。ただし、上式において :math:`\beta` は逆温度 :math:`1/ k_B T` , :math:`\mu` は化学ポテンシャルとしました。 
なお、粒子数を固定するカノニカル計算では、粒子数を :math:`N` とした場合に、

.. math::
   \begin{aligned}
   N = \sum_{i} \langle c_i^{\dagger} c_i \rangle
   \end{aligned}

を満たすように、各ステップで :math:`\mu` を決定・更新しながら計算を進めます。H-waveでは、更新に対してはsimple-mixingを現在行っています。simple-mixingでは :math:`n` 番目のステップの一体グリーン関数を :math:`\langle c_{i}^\dagger c_{j}\rangle^{(n)}`  とした場合に、 :math:`n` 番目の一体グリーン関数と :math:`n+1` 番目のステップで求められた一体グリーン関数を混ぜた上で更新を行います：

.. math::
   \begin{aligned}
   \langle c_{i}^\dagger c_{j}\rangle^{(n+1)} = (1-\alpha) \langle c_{i}^\dagger c_{j}\rangle^{(n)} +  \alpha \langle c_{i}^\dagger c_{j}\rangle^{(n+1)} 
   \end{aligned}

ここで、 :math:`\alpha` は0から1までの定数を表します。現在のH-waveでは実装していませんが、simple-mixing以外の更新方法としては、Anderson-mixingなどもあります。
なお、H-waveの実空間UHFでは、全ての相互作用をInterAll形式にマップした後に解析を行う仕様になっています。
有限温度の自由エネルギーは

.. math::
   \begin{aligned}
   F = \mu N -\frac{1}{\beta}\sum_k \ln \left[ 1+\exp (-\beta(\xi_k - \mu)) \right]
    - \sum_{ijkl} I_{ijkl} (\langle c_{i}^{\dagger} c_j\rangle \langle c_{k}^{\dagger} c_l\rangle - \langle c_{i}^{\dagger} c_l\rangle \langle c_{k}^{\dagger} c_j\rangle)
   \end{aligned}

で与えられます。

波数空間への拡張
*************************

一体項で与えられたハミルトニアンについて、Fourier変換 :math:`c_i = \dfrac{1}{\sqrt{V}} \sum_k e^{ikr_i} c_k` を用いて
波数表示に書き換えます。

.. math::
   \begin{aligned}
   \mathcal{H}_\text{UHF}
   &=
   \sum_{k\alpha\beta\sigma\sigma^\prime}
   h_{\alpha\beta\sigma\sigma^\prime}(k)\,
   c_{k\alpha\sigma}^\dagger c_{k\beta\sigma^\prime}^{\phantom\dagger}
   \end{aligned}

相互作用は並進対称性を仮定し、係数 :math:`J_{ij\alpha\beta}` は空間座標について並進ベクトル :math:`r_{ij}=r_j - r_i` のみに
依存するとします。なお、波数空間版UHFでは一般化された InterAll 形式の相互作用は扱いません。

波数空間表示のハミルトニアンは波数 :math:`k` について対角的なので、固有値・固有ベクトルの計算は、
:math:`N_\text{site}N_\text{orbit} \times N_\text{site}N_\text{orbit}` 行列の対角化から、
:math:`N_\text{site}` 個の :math:`N_\text{orbit} \times N_\text{orbit}` 行列の対角化に計算量を抑えることができます。ここで :math:`N_\text{site}` はサイト数、 :math:`N_\text{orbit}` はスピン自由度も含めた1サイト当たりの軌道縮絨度です。
