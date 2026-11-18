.. highlight:: none

.. _subsec:scf_status.dat:

scf_status
~~~~~~~~~~

自己無撞着(SCF)計算の収束状態 (収束したか、最後の残差、収束判定値、ステップ数) を出力します。
このファイルはデフォルトで出力ディレクトリに\ ``scf_status.dat``\ として、他のどの結果ファイルよりも先に出力されます。
ファイル名は\ ``file.output``\ セクションのキーワード\ ``scf_status``\ で指定でき、\ ``scf_status = ""``\ とすると出力されません。

以下にファイル例を記載します。

::

    format_version = 1
    solver = UHFr
    converged = true
    residual = 5.755284863831164e-09
    eps = 1e-08
    iterations = 34
    iteration_max = 1000
    residual_definition = two_frobenius_over_size:uhfr_full
    note =

ファイル形式
^^^^^^^^^^^^

各行はちょうど\ ``key = value``\ と改行からなり、\ ``=``\ の両側に空白を1つずつ置き、それ以外の空白は入りません。
キーは常に以下の順に並びます。値が空のときは\ ``key =``\ の後に空白1つと改行だけが続きます (通常の計算の\ ``note =``\ など)。

-  ``format_version = 1``

-  ``solver = UHFr``

-  ``converged = [converged]``

-  ``residual = [residual]``

-  ``eps = [eps]``

-  ``iterations = [iterations]``

-  ``iteration_max = [iteration_max]``

-  ``residual_definition = two_frobenius_over_size:uhfr_full``

-  ``note = [note]``

パラメータ
^^^^^^^^^^

-  ``[converged]``

   **形式 :** ``true``\ または\ ``false``

   **説明 :**
   最後の残差が\ ``eps``\ 未満、すなわち収束してループが終了した場合に\ ``true``\ となります。

-  ``[residual]``

   **形式 :** float型 (Pythonの\ ``repr``\ による全桁表示) または\ ``nan``

   **説明 :**
   最後に完了したステップの残差

   .. math::

      R = \frac{2\,\lVert G_{\rm new} - G_{\rm old} \rVert_F}{\mathrm{size}(G)}

   です。ここで\ :math:`G_{\rm old}`\ はステップ開始時のグリーン関数、\ :math:`G_{\rm new}`\ はそのステップで計算されたもの (混合前)、\ :math:`\lVert\cdot\rVert_F`\ はFrobeniusノルム、\ :math:`\mathrm{size}(G)`\ は\ :math:`G`\ の要素数です。
   UHFrでは\ :math:`G`\ は\ :math:`N`\ サイトと両スピンからなる\ :math:`2N \times 2N`\ の行列全体で、\ :math:`R = \sqrt{\sum|\Delta G|^2}/(2N^2)`\ となります。
   UHFkは同じ式をすべての並進ベクトルにわたる配列に適用します。

   :math:`R`\ はノルムを要素数で割っているため、系のサイズに対してスケールしません。
   同じ\ ``EPS``\ でも系の大きさによって判定条件の厳しさが変わるので、残差は同じサイズの計算どうしでのみ比較してください。

   ステップを1回も行わなかった場合、または最後のステップの残差が有限でなかった場合 (\ ``[note]``\ を参照) は\ ``nan``\ です。

-  ``[eps]``

   **形式 :** float型

   **説明 :** 収束判定値\ :math:`10^{-\mathrm{EPS}}`\ (\ ``mode.param``\ の\ ``EPS``\ )。

-  ``[iterations]``

   **形式 :** int型

   **説明 :** 完了したSCFステップの数 (収束したステップを含む)。\ ``print_check``\ ファイルを指定した場合はその行数と一致します。

-  ``[iteration_max]``

   **形式 :** int型

   **説明 :** ``mode.param``\ の\ ``IterationMax``\ 。

-  ``[note]``

   **形式 :** 文字列 (空の場合あり)

   **説明 :**
   通常の計算では空です。それ以外では次のいずれかです。

   - ``no SCF step performed``\ : ``IterationMax = 0``\ の場合。
   - ``residual not finite after step <n>``\ : 最後に完了したステップ\ ``<n>``\ の残差がNaNまたは無限大だった場合。このとき\ ``residual``\ は\ ``nan``\ 、\ ``converged``\ は\ ``false``\ です。
   - ``IterationMax <value> treated as 0``\ : 負の\ ``IterationMax``\ (\ ``trustme_interaction_range``\ でパラメータチェックを緩めた場合にのみ受け付けられます) はステップを行わず、\ ``iteration_max = 0``\ として出力されます。

不変条件
^^^^^^^^

-  ``eps``\ は有限の正の値です。

-  :math:`0 \le`\ ``iterations``\ :math:`\le`\ ``iteration_max``\ です。

-  ``iterations``\ > 0 のとき、\ ``residual``\ は有限の非負の値か、または注記\ ``residual not finite after step <n>``\ を伴う\ ``nan``\ です。

-  ``converged = true``\ のとき、\ ``iterations``\ > 0 かつ\ ``residual``\ < ``eps``\ です。

-  ``iterations = 0``\ のとき、\ ``converged = false``\ かつ\ ``residual = nan``\ です。

SCFステップを行わない場合
^^^^^^^^^^^^^^^^^^^^^^^^^^

``IterationMax = 0``\ の場合はステップを行いません。WARNINGが1つ出力され、結果ファイルは\ ``scf_status.dat``\ のみが出力されます (\ ``iterations = 0``\ 、\ ``converged = false``\ 、\ ``residual = nan``\ 、注記\ ``no SCF step performed``\ )。
``energy``\ 、\ ``eigen``\ 、\ ``green``\ などの他の出力は行われません。
``print_check``\ ファイルを指定している場合は空のファイルが作成されます。
``scf_status = ""``\ の場合は結果ファイルが1つも出力されず、WARNINGにもそのように記されます。
以前のバージョンでは、この場合はエラーで停止していました。

UHFrの\ ``file.output``\ の他のキーは指定した場合にのみ出力されますが、\ ``scf_status``\ は\ ``""``\ を指定しない限り出力されます。

``log``\ セクションの\ ``print_check``\ (各ステップのステップ番号、残差、エネルギー、\ ``NCond``\ 、\ ``Sz``\ ) はこれまでどおり利用できます。

.. note::

   キー名\ ``scf_status``\ は、\ ``file.output``\ の他のキーと同様に完全一致で照合されます (\ ``SCF_Status``\ は別のキーとして扱われます)。
   値は文字列でなければならず、それ以外の値 (\ ``false``\ 、\ ``0``\ など) は計算開始前にエラーとなります。
