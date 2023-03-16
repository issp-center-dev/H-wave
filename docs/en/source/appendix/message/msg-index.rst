エラーメッセージ一覧
=====================

- ``mode is not defined in [mode].``

  **説明 :** パラメータファイルの ``[mode]`` セクションに ``mode`` パラメータが指定されていない。

  **モード :** main

- ``Get_param: key must be mod or ham or output.``

  **説明 :** ``get_param()`` の引数が不正

  **モード :** UHFr (read_input)

- ``duplicate items found in`` *file*

  **説明 :** *file* に重複するエントリがある

  **モード :** UHFr (read_input)

- ``incorrect number of lines in`` *file* ``: expected=`` *N* ``, found=`` *M*

  **説明 :** 入力ファイルの行数とファイル内の行数指定が異なっている

  **モード :** UHFr (read_input)

- ``Unknown keyword`` *keyword*

  **説明 :** ``[file.input.interaction]`` に不明なキーワードがある

  **モード :** UHFk (read_input_k)

- ``initial and initial_uhf can not be specified simultaneously.``

  **説明 :** ``initial`` と ``initial_uhf`` は同時に指定できない

  **モード :** UHFk (read_input_k)

- ``read_input_k: file`` *file* ``not found``

  **説明 :** *file* が存在しない

  **モード :** UHFk (read_input_k)

- ``Get_param: key must be mod or ham or output.``

  **説明 :** ``get_param()`` の引数が不正

  **モード :** UHFk (read_input_k)

- ``read_geom: file`` *file* ``not found``

  **説明 :** ``Geometry`` に指定されているファイル *file* が存在しない

  **モード :** UHFk (wan90)

- ``mode.param.2Sz must be even(odd) when Ncond is even(odd).``

  **説明 :**  パラメータ ``2Sz`` と ``Ncond`` の偶奇が一致していない

  **モード :** solver base

- ``range check for`` *type* ``failed.``

  **説明 :** *type* の値が範囲外

  **モード :** UHFk

- ``_check_cellsize failed. interaction range exceeds cell shape.``

  **説明 :** 相互作用項の並進ベクトルが CellShape 内に収まっていない

  **モード :** UHFk

- ``Hermiticity check failed: |T_ba(-r)^* - T_ab(r)| =`` *val*

  **説明 :** Transfer項が Hermite でない

  **モード :** UHFk

- ``Parameter range check failed for param_mod.``

  **説明 :** [mode.param] のパラメータの値が範囲外

  **モード :** solver base

- ``Parameter check failed for param_mod.``

  **説明 :** [mode.param] のパラメータが不正

  **モード :** solver base

- ``Hermite check failed for`` *type*

  **説明 :** *type* が Hermite でない

  **モード :** UHFr

- ``Parameter check failed for info_mode.``

  **説明 :** [mode] のパラメータが不正

  **モード :** solver base

- ``value not integer``

  **説明 :** パラメータの値が整数ではない

  **モード :** RPA

- ``Lattice initialization failed: 'CellShape' not found.``

  **説明 :** [mode.param] に ``CellShape`` が指定されていない

  **モード :** RPA

- ``Ncond must be greater than zero: Ncond=`` *Ncond*

  **説明 :** ``Ncond`` に0以上の値が指定されていない

  **モード :** RPA

- ``Nmat must be greater than zero: Nmat=`` *Nmat*

  **説明 :** ``Nmat`` に0以上の値が指定されていない

  **モード :** RPA

- ``RPA._find_mu: not converged. abort``

  **説明 :** ``mu`` の計算が収束しなかった

  **モード :** RPA

- ``SubShape is not compatible with CellShape.``

  **説明 :** 副格子の指定が不正

  **モード :** RPA

- ``T must be greater than or equal to zero: T=`` *T*

  **説明 :** ``T`` に 0以上の値が指定されていない

  **モード :** RPA

- ``both mu and Ncond or filling are specified``

  **説明 :** ``mu` と ``Ncond`` または ``filling`` が同時に指定されている

  **モード :** RPA

- ``dimension of CellShape must be one, two, or three.``

  **説明 :** ``CellShape`` の指定が不正

  **モード :** RPA

- ``dimension of SubShape does not match with that of CellShape.``

  **説明 :** 副格子の指定が不正

  **モード :** RPA

- ``invalid CellShape.``

  **説明 :** CellShapeの指定が不正

  **モード :** RPA

- ``invalid SubShape.``

  **説明 :** SubShapeの指定が不正

  **モード :** RPA

- ``none of mu, Ncond, nor filling is specified``

  **説明 :** ``mu`` または ``Ncond``, ``filling`` のいずれも指定されていない

  **モード :** RPA

- ``read_chi0q failed:`` *info*

  **説明 :** ``chi0q`` の読み込みに問題があった

  **モード :** RPA

- ``round_to_int: unknown mode`` *mode*

  **説明 :** 丸めモードの指定が不正

  **モード :** RPA

- ``unexpected data size`` *error*

  **説明 :** データサイズが不正

  **モード :** RPA

- ``mode is not defined in [mode].``

  **説明 :** ``[mode]`` に mode パラメータが指定されていない

  **モード :** RPA

- ``orbital index check failed for`` *type*

  **説明 :** 軌道のインデックスが不正

  **モード :** UHFk

- ``initial green function in coord space requires geometry.dat``

  **説明 :** 実空間でのグリーン関数の読み込みは ``geometry.dat`` を同時に指定する必要がある

  **モード :** UHFk

- ``CellShape is missing. abort``

  **説明 :** ``CellShape`` パラメータが指定されていない

  **モード :** UHFk

- ``Ncond or Nelec is missing. abort``

  **説明 :** ``Ncond`` または ``Nelec`` パラメータが指定されていない

  **モード :** UHFk

- ``SubShape is not compatible with CellShape. abort``

  **説明 :** 副格子の指定が不正

  **モード :** UHFk

- ``_check_orbital_index failed. invalid orbital index found in interaction definitions.``

  **説明 :** 相互作用定義ファイルの軌道インデックスが不正

  **モード :** UHFk

- ``_save_greenone: onebodyg_uhf and geometry_uhf are required``

  **説明 :** ``onebodyg_uhf`` と ``geometry_uhf`` が指定されていない

  **モード :** UHFk

- ``find mu: not converged. abort``

  **説明 :** ``mu`` の計算が収束しなかった

  **モード :** UHFk

- ``range check failed for Initial``

  **説明 :** ``Initial`` の指定が不正

  **モード :** UHFr

- ``OneBodyG is required to output green function.``

  **説明 :** グリーン関数の出力のための ``OneBodyG`` の指定がない

  **モード :** UHFr

- ``hermite check failed for Initial``

  **説明 :** ``Initial`` が Hermite でない

  **モード :** UHFr

- ``Range check failed for Transfer``

  **説明 :** ``Transfer`` 定義ファイルのインデックスが範囲外

  **モード :** UHFr

- ``Range check failed for`` *type*

  **説明 :** *type* 定義ファイルのインデックスが範囲外

  **モード :** UHFr

- ``parameter range check failed.``

  **説明 :** パラメータの値が不正

  **モード :** UHFr

- ``mode is incorrect: mode=`` *mode*

  **説明 :** ``mode`` の指定が不正

  **モード :** UHFr

- ``mode.param.`` *key* ``must be greater than`` *value*

  **説明 :** パラメータ *key* の値が不正

  **モード :** solver base [warning]

- ``"mode.param.`` *key* ``must be smaller than`` *value*

  **説明 :** パラメータ *key* の値が不正

  **モード :** solver base [warning]

- ``mode.param.`` *key* ``is not defined.``

  **説明 :** パラメータ *key* が設定されていない

  **モード :** solver base [warning]

- ``mode.`` *key* ``in mode section is incorrect:`` *values*

  **説明 :** ``[mode]`` セクションの ``mode`` パラメータの値が不正

  **モード :** solver base [warning]

- ``mode.`` *key* ``is not defined.``

  **説明 :** ``[mode]`` セクションに ``mode`` パラメータが指定されていない

  **モード :** solver base [warning]

- ``TRUST-ME mode enabled. parameter checks are relaxed``

  **説明 :** ``TRUST-ME`` モードが有効。パラメータのチェックを行わない

  **モード :** solver base [warning]

- ``value not integer``

  **説明 :** 設定値が整数でない

  **モード :** RPA [warning]

- ``mode is incorrect: mode=`` *mode*

  **説明 :** ``mode`` パラメータの値が不正

  **モード :** RPA [warning]

- ``FATAL: 2Sz=`` *value* ``. 2Sz should be even for calculating fij``

  **説明 :** :math:`f_{ij}` の計算で ``2Sz`` は偶数でなければならない

  **モード :** UHFr [warning]

- ``FATAL: Ne=`` *value* ``. Ne should be even for calculating fij``

  **説明 :** :math:`f_{ij}` の計算で ``Ne`` は偶数でなければならない

  **モード :** UHFr [warning]

- ``NOT IMPLEMENTED: Sz even and Sz != 0: this case will be implemented in near future``

  **説明 :** :math:`f_{ij}` の計算で ``Sz`` が 0以外の偶数の場合は未サポート

  **モード :** UHFr [warning]

- ``key`` *key* ``is wrong!``

  **説明 :** キーワード *key* が不正

  **モード :** UHFr [warning]

- ``UHFr calculation is failed: rest=`` *residue* ``, eps=`` *eps*

  **説明 :** UHFr の計算が収束しなかった

  **モード :** UHFr [warning]

