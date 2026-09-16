.. highlight:: none

相互作用指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

次の形で表わされるハミルトニアンの一体項および二体相互作用項について、その係数\ :math:`T_{\alpha\beta}(r_{ij})`\ および\ :math:`J_{\alpha\beta}(r_{ij})`, :math:`V_{\alpha\beta}(r_{ij})`, :math:`U_{\alpha}`\ を共通のWannier90(-like)形式で記述します。
なお、波数空間版UHFでは一般化二体相互作用 InterAll 形式には対応していません。

    
    **Transfer**:
      :math:`\sum_{ij\alpha\beta\sigma} T_{\alpha\beta}(r_{ij})\,c_{i\alpha\sigma}^{\dagger}c_{j\beta\sigma}^{\phantom{\dagger}}`
    **CoulombIntra**:
      :math:`\sum_{i\alpha} U_\alpha\,n_ {i\alpha\uparrow} n_{i\alpha\downarrow}` (\ :math:`n_{i\alpha\sigma}=c_{i\alpha\sigma}^{\dagger}c_{i\alpha\sigma}^{\phantom{\dagger}}`)
    **CoulombInter**:
      :math:`\sum_{ij\alpha\beta} V_{\alpha\beta}(r_{ij})\,n_{i\alpha} n_{j\beta}` (\ :math:`n_{i\alpha}=n_{i\alpha\uparrow}+n_{i\alpha\downarrow}`)
    **Coulomb**:
      CoulombIntra と CoulombInter をまとめた形式（RESPACK の\ ``zvo_ur.dat``\ を想定）。:math:`r=0`\ かつ同一軌道（:math:`\alpha=\beta`\ ）の成分を CoulombIntra、それ以外を CoulombInter として読み込みます。CoulombIntra/CoulombInter を個別に指定するのと等価です。
    **Hund**:
      :math:`-\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Hund}(r_{ij}) \left( n_{i\alpha\uparrow} n_{j\beta\uparrow} + n_{i\alpha\downarrow} n_{j\beta\downarrow} \right)`
    **Ising**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ising}(r_{ij}) (n_{i\alpha\uparrow} - n_{i\alpha\downarrow})(n_{j\beta\uparrow} - n_{j\beta\downarrow})`
    **PairHop**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm PH}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{j\beta\uparrow}^{\phantom{\dagger}} c_{i\alpha\downarrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}} + h.c.`
    **Exchange**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ex}(r_{ij})\,c_{i\alpha\uparrow}^\dagger c_{j\beta\uparrow}^{\phantom{\dagger}} c_{j\beta\downarrow}^\dagger c_{i\alpha\downarrow}^{\phantom{\dagger}}`
    **PairLift**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm PairLift}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow}^{\phantom{\dagger}} c_{j\beta\uparrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}} + h.c.`

.. note::

   上の式のうち2つは、ソルバーが以前から実装している内容に合わせて修正した
   ものです。\ **入力ファイルを変更する必要はなく**\ 、数値結果も変わりません。

   - **Hund**: 実装されている規約は上のとおりマイナス符号を伴います。
     したがって強磁性的な（同スピン間で引力的な）Hund 結合は、係数
     :math:`J^{\rm Hund}_{\alpha\beta}(r_{ij})`\ を正の値として宣言します。
   - **PairLift**: 宣言した各行はそのエルミート共役とともにハミルトニアンに
     含まれます（``PairHop``\ と同様）。改訂前の本ページでは\ ``+ h.c.``\ の
     記載が漏れていました。

.. _uhfk_interaction_orientation:

.. note::

   **規約：**\ 上の全ての式で\ :math:`r_{ij} = R_j - R_i`\ です。すなわち
   ``[rx] [ry] [rz] [alpha] [beta] ...``\ という行は、軌道\ ``[alpha]``\ を
   元のセルに、軌道\ ``[beta]``\ を\ :math:`\vec{r} = (r_x, r_y, r_z)`\ だけ
   並進したセルに置きます。これは一体項のファイルでも二体項のファイルでも
   同じです。

   **2.0.0 との互換性：**\ H-wave 2.0.0 までは、平均場ソルバー（``UHFk``\ 、
   および\ ``FLEX``\ の Hartree-Fock 項）はオフサイトの二体行を2つのセルを
   入れ替えて読んでいましたが、RPA・FLEX・Eliashberg の頂点は上の規約に
   従っていました。今回から平均場も上の規約に従います。\ ``UHFk``\ の結果
   （および\ ``flex_hartree_fock = true``\ または
   ``flex_second_order = "local"``\ を指定した\ ``FLEX``\ の結果）が変わるのは、
   異なる2つの軌道をもつオフサイト行、または複素係数をもつオフサイト行が
   ある場合だけです。全てのオンサイト相互作用、実数で軌道対角な全ての
   オフサイト行、実数の単一軌道入力、およびボンド分解チャネルを用いない
   （``longitudinal_bond_channels = false``\ 、既定値）、かつ影響を受けない
   感受率から計算した RPA と Eliashberg の全ての結果は変わりません
   （``chi0q_mode = "flex"``\ の Eliashberg 計算は、感受率を生成した FLEX 計算
   に従って変わります）。同じ相互作用ファイル形式を用いる FLEX/RPA 計算に
   ついては：\ ``calc_type = "ring"``\ と
   ``longitudinal_bond_channels = true``\ を指定した計算も、実数の軌道間
   オフサイトボンドでは結果が変わります（混合した2次のブロックが厳密に
   なりました。issue #192）。Hartree-Fock 項とボンド分解チャネルのいずれも
   使わない\ ``flex_second_order = "takimoto"``\ の計算は変わりません。\ **入力ファイルを変更する必要はありません。**

   **再現手順と複合パイプライン：**\ 2.0.0 の\ ``UHFk``\ の数値を再現するには、
   全てのオフサイト二体行の変位を反転してください
   （``[rx] [ry] [rz]`` → ``[-rx] [-ry] [-rz]``\ 。一体項の\ ``Transfer``\ の行、およびオンサイトの行は
   対象外です。また、軌道の
   添字を入れ替えてはいけません。複素結合に対して誤りになります）。ただし、
   このように変位を反転したファイルは、RPA と Eliashberg のソルバーが計算する
   内容を変えてしまいます。元のファイルに対するこれらの 2.0.0 の結果は、
   もともと正しいものでした。\ ``flex_hartree_fock = true``\ を指定した 2.0.0 の
   ``FLEX``\ 計算は、単一の入力ファイルでは再現できません（2.0.0 は一回の計算
   の中で平均場項と頂点項を異なる向きで評価していたためです）。\ ``UHFk``\ の
   出力を RPA や Eliashberg のソルバーに渡していた 2.0.0 のパイプラインを
   再現するには、\ ``UHFk``\ だけを変位を反転したファイルで、後段は元の
   ファイルで実行する必要があります。2.0.0 では、2つの段階が異なる
   ハミルトニアンを解いていたためです。従来の食い違いに対する回避策を入れて
   いた場合は取り除き、影響を受ける自己無撞着計算は 2.0.0 の初期値からでは
   なく最初からやり直してください。


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

-  4 - :math:`\lceil N_\text{pts} / 15 \rceil + 3`\ 行:
   ``[n_1] [n_2] ...``

-  :math:`\lceil N_\text{pts} / 15 \rceil + 4`\ 行以降:
   ``[r_x] [r_y] [r_z] [alpha] [beta] [J.real] [J.imag]``

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

-  ``[n1]``, ``[n2]``, ...

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
   ``[alpha]``\ が元のセル内の軌道、``[beta]``\ が\ :math:`\vec{r}`\ 離れたセル内の軌道を指します。

-  ``[J.real]``, ``[J.imag]``

   **形式 :** float型

   **説明 :**
   係数\ :math:`J_{\alpha\beta}(\vec{r})`\ の実部と虚部を指定します。


使用ルール
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行うため、ヘッダの省略はできません。

-  係数行列のうち、省略された要素は 0と仮定します。なお、エルミート共役の相手\ :math:`X_{ba}(-R)`\ が省略された宣言済みエントリは読み込み時に拒否されます（バージョン 2.0 以降）。結合の両方向を宣言してください。

-  並進ベクトルは全て\ ``CellShape``\ 内に収まるとします。
   ``r_x``, ``r_y``, ``r_z``\ の範囲が\ ``CellShape``\ のx,y,z軸のサイズを超える場合はエラーで終了します。

-  ``mode.enable_spin_orbital``\ が\ ``true``\ の場合、Transfer項の軌道のインデックスはスピン自由度を含む一般化軌道インデックスと読み替え、1〜\ :math:`2 N_\text{orbital}`\ （ジオメトリファイルの\ ``Norbit``\ 。このモードではスピン軌道の数）の値をとります。スピンを内側に並べたインターリーブ順で、奇数インデックス (1, 3, 5, …) が各軌道の spin up、偶数インデックス (2, 4, 6, …) が spin down に対応します（軌道\ :math:`\alpha` (0 起点) とスピン\ :math:`s` (0: up, 1: down) に対しファイル上のインデックス（1 起点）は\ :math:`2\alpha + s + 1`\ ）。\ ``mode.enable_spin_orbital``\ が\ ``false``\ の場合は、インデックスの範囲が 1〜\ :math:`N_\text{orbital}`\ の行のみ考慮します。

-  スピン軌道モードでは相互作用項（CoulombIntra, CoulombInter, Coulomb, Hund, Ising, Exchange, PairLift, PairHop）も利用でき、仮想スピン分解を介して取り扱われます。倍化された\ ``2α+s+1``\ のインデックス規約は Transfer ファイルのみに適用されます。相互作用定義ファイルでは物理軌道のインデックス（1 〜\ :math:`N_\text{orbital}`\ 、すなわちジオメトリファイルの\ ``Norbit``\ の半分）を用います。

.. note::

   Ising 項の規約の履歴。（1）本ページの旧版はこの項を
   :math:`J S^z S^z`\ （\ :math:`S^z = (n_\uparrow - n_\downarrow)/2`\ ）
   と定義しており、UHFk ソルバーはその定義に従っていました（上記の形式
   に対して実効的に因子 1/4）。（2）RPA/FLEX ソルバーは一貫して上記の
   密度差形式でファイルを読んでおり、これは厳密対角化で頂点内容が判定
   された演算子でもあります。（3）現在は wannier90 形式の k 空間ソルバーとページはすべて密度差
   形式を共有します（UHFr の実空間リーダーは独自の S^z 規約を保持し
   ます）。旧バージョンの UHFk の結果は、同じファイルが現在
   与える結合の 4 分の 1 に対応します。

.. note::

   二体相互作用の宣言ファイル（Coulomb, CoulombIntra, CoulombInter,
   Hund, Exchange, Ising, PairLift, PairHop）はエルミート共役で閉じて
   いる必要があります：
   各エントリ\ :math:`X_{ab}(R)`\ には相手\ :math:`X_{ba}(-R) = X_{ab}(R)^{*}`
   が伴わなければなりません。相手が欠けている、または値が一致しない場合
   は読み込み時に拒否されます。 CoulombIntra はさらにオンサイト同一軌道・有限実数値の
   エントリのみを受け付けます。（これはここで読む wannier90 形式の
   k 空間ファイルに適用されます。別系統の UHFr 実空間リーダーは独自の
   規約を保持します。）
