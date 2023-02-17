.. highlight:: none

相互作用指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

次の形で表わされるハミルトニアンの一体項および二体相互作用項について、その係数 :math:`T_{\alpha\beta}(r_{ij})` および :math:`J_{\alpha\beta}(r_{ij})`, :math:`V_{\alpha\beta}(r_{ij})`, :math:`U_{\alpha}` を共通のWannier90(-like)形式で記述します。
なお、波数空間版UHFでは一般化二体相互作用 InterAll 形式には対応していません。

    
    **Transfer**:
      :math:`\sum_{ij\alpha\beta\sigma} T_{\alpha\beta}(r_{ij})\,c_{i\alpha\sigma}^{\dagger}c_{j\beta\sigma}^{\phantom{\dagger}}`
    **CoulombIntra**:
      :math:`\sum_{i\alpha} U_\alpha\,n_ {i\alpha\uparrow} n_{i\alpha\downarrow}` (\ :math:`n_{i\alpha\sigma}=c_{i\alpha\sigma}^{\dagger}c_{i\alpha\sigma}^{\phantom{\dagger}}`)
    **CoulombInter**:
      :math:`\sum_{ij\alpha\beta} V_{\alpha\beta}(r_{ij})\,n_{i\alpha} n_{j\beta}` (\ :math:`n_{i\alpha}=n_{i\alpha\uparrow}+n_{i\alpha\downarrow}`)
    **Hund**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Hund}(r_{ij}) \left( n_{i\alpha\uparrow} n_{j\beta\uparrow} + n_{i\alpha\downarrow} n_{j\beta\downarrow} \right)`
    **Ising**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ising}(r_{ij}) S^{z}_{i\alpha} S^{z}_{j\beta}` (\ :math:`S^{z}_{i\alpha}=\frac{1}{2}(n_{i\alpha\uparrow} - n_{i\alpha\downarrow})`)
    **PairHop**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm PH}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{j\beta\uparrow}^{\phantom{\dagger}} c_{i\alpha\downarrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}} + h.c.`
    **Exchange**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ex}(r_{ij})\,c_{i\alpha\uparrow}^\dagger c_{j\beta\uparrow}^{\phantom{\dagger}} c_{j\beta\downarrow}^\dagger c_{i\alpha\downarrow}^{\phantom{\dagger}}`
    **PairLift**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm PairLift}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow}^{\phantom{\dagger}} c_{j\beta\uparrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}}`


以下にファイル例を示します。

::

   wannier90 format for vmcdry.out or HPhi -sdry
       10
      245
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
   ...(略)
    1    1    1    1    1
   -3   -3   -2    1    1  -0.0000269645  -0.0000000000
   -3   -3   -2    1    2  -0.0000071722  -0.0000018600
   -3   -3   -2    1    3  -0.0000083990   0.0000010972
   -3   -3   -2    1    4  -0.0000000990   0.0000000427
   -3   -3   -2    1    5  -0.0000018628  -0.0000003609
   -3   -3   -2    1    6  -0.0000129504  -0.0000014047
   -3   -3   -2    1    7  -0.0000189169   0.0000024697
   -3   -3   -2    1    8   0.0000238115   0.0000014316
   -3   -3   -2    1    9   0.0000036708  -0.0000003266
   -3   -3   -2    1   10   0.0000361752   0.0000003247
   -3   -3   -2    2    1  -0.0000071722   0.0000018600
   -3   -3   -2    2    2   0.0000105028  -0.0000000000
   ...(略)


ファイル形式
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: ``[Norbit]``

-  3行: ``[Npts]``

-  4 - :math:`\lceil N_\text{pts} / 15 \rceil + 3` 行:
   ``[n_1] [n_2] ...``

-  :math:`\lceil N_\text{pts} / 15 \rceil + 4` 行以降:
   ``[r_x] [r_x] [r_x] [alpha] [beta] [J.real] [J.imag]``

パラメータ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``[Norbit]``

   **形式 :** int型

   **説明 :**
   ユニットセル内の軌道の数\ :math:`N_\text{orbit}`\ を指定します。

-  ``[Npts]``

   **形式 :** int型

   **説明 :**
   並進ベクトル全体が入る直方体に含まれるセルの数を指定します。

-  ``[n1]``, ``[n1]``, ...

   **形式 :** int型

   **説明 :**
   各セルの縮重度を指定します(通常は 1)。一行あたり15点を列挙します。

-  ``[r_x]``, ``[r_y]``, ``[r_z]``

   **形式 :** int型

   **説明 :**
   並進ベクトルを指定します。
   
-  ``[alpha]``, ``[beta]``

   **形式 :** int型

   **説明 :**
   軌道のインデックスを指定します。
   ``[alpha]`` が元のセル内の軌道、``[beta]`` が :math:`\vec{r}` 離れたセル内の軌道を指します。

-  ``[J.real]``, ``[J.imag]``

   **形式 :** float型

   **説明 :**
   係数 :math:`J_{\alpha\beta}(\vec{r})` の実部と虚部を指定します。


使用ルール
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行うため、ヘッダの省略はできません。

-  係数行列のうち、省略された要素は 0と仮定します。

-  並進ベクトルは全て ``CellShape`` 内に収まるとします。
   ``r_x``, ``r_y``, ``r_z`` の範囲が ``CellShape`` のx,y,z軸のサイズを超える場合はエラーで終了します。

.. raw:: latex
