.. highlight:: none

.. _Ch:HowToWannier90_rpa:

RPA用入力ファイル
--------------------------------

乱雑位相近似で使用する入力ファイルに関して説明します。
入力ファイルの種別は以下の2つに分類され、データフォーマットは Wannier90 形式です。

(1) 幾何情報:
    
    **Geometry**: 格子の幾何情報を設定します。

(2) 相互作用定義:
    
    RPAのハミルトニアンを電子系の表式により指定します。
    具体的には、キーワードで指定する各相互作用のタイプについて、ハミルトニアンの係数をデータとして与えます。

    :math:`{\mathcal H}\Phi` や mVMC のExpert Mode入力に相当する以下のキーワードを指定できます。

    **Transfer**:
      :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`\ で表される一体項を指定します。
    **CoulombIntra**:
      :math:`n_ {i \uparrow}n_{i \downarrow}`\ で表される相互作用を指定します(\ :math:`n_{i \sigma}=c_{i\sigma}^{\dagger}c_{i\sigma}`)。
    **CoulombInter**:
      :math:`n_ {i}n_{j}`\ で表される相互作用を指定します(\ :math:`n_i=n_{i\uparrow}+n_{i\downarrow}`)。
    **Hund**:
      :math:`n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}`\ で表される相互作用を指定します。
    **Ising**:
      :math:`S_i^z S_j^z`\ で表される相互作用を指定します。  
    **Exchange**:
      :math:`S_i^+ S_j^-`\ で表される相互作用を指定します。  
    **PairLift**:
      :math:`c_ {i \uparrow}^{\dagger}c_{i\downarrow}c_{j \uparrow}^{\dagger}c_{j \downarrow}`\ で表される相互作用を指定します。	    
    **PairHop**:
      :math:`c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow}`\ で表される相互作用を指定します。

以下のセクションでデータフォーマットについて記述します。

.. toctree::
  :maxdepth: 1

  geometry
  interaction

