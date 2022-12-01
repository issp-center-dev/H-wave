.. highlight:: none

.. _Ch:HowToExpert:

UHFr用入力ファイル
--------------------------------

H-waveで使用する入力ファイル(\*.def)に関して説明します。
入力ファイルの種別は以下の2つで分類されます。

(1) Hamiltonian:
    
    Hamiltonianを電子系の表式により指定します。
    具体的には以下のファイルで指定されます。
    
    | **Trans**:
      :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}`\ で表される一体項を指定します。
    | **InterAll**:
      :math:`c_ {i \sigma_1}^{\dagger}c_{j\sigma_2}c_{k \sigma_3}^{\dagger}c_{l \sigma_4}`\ で表される一般二体相互作用を指定します。
	    
    | なお、使用頻度の高い相互作用に関しては下記のキーワードで指定することも可能です。
    | **CoulombIntra**:
      :math:`n_ {i \uparrow}n_{i \downarrow}`\ で表される相互作用を指定します(\ :math:`n_{i \sigma}=c_{i\sigma}^{\dagger}c_{i\sigma}`)。
    | **CoulombInter**:
      :math:`n_ {i}n_{j}`\ で表される相互作用を指定します(\ :math:`n_i=n_{i\uparrow}+n_{i\downarrow}`)。
    | **Hund**:
      :math:`n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}`\ で表される相互作用を指定します。
    | **PairHop**:
      :math:`c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow}`\ で表される相互作用を指定します。
    | **Exchange**:
      :math:`S_i^+ S_j^-`\ で表される相互作用を指定します。
    | **Ising**:
      :math:`S_i^z S_j^z`\ で表される相互作用を指定します。
    | **PairLift**:
      :math:`c_ {i \uparrow}^{\dagger}c_{i\downarrow}c_{j \uparrow}^{\dagger}c_{j \downarrow}`\ で表される相互作用を指定します。
	    
(2) Green functions:

    | **Initial** :初期値として入力する一体Green関数を指定します。
      :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}\rangle`\ を入力します。
    | **OneBodyG** :出力する一体Green関数を指定します。
      :math:`\langle c^{\dagger}_{i\sigma_1}c_{j\sigma_2}\rangle`\ が出力されます。

以下、各入力ファイルのファイルフォーマットについて詳細を記載します。

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Trans_file
   InterAll_file
   CoulombIntra_file
   CoulombInter_file
   Hund_file
   PairHop_file
   Exchange_file
   Ising_file
   PairLift_file
   Initial_file
   OneBodyG_file
