.. highlight:: none

.. _Subsec:InputFileList:


入力ファイル指定用ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~

| 計算で使用する入力ファイル一式を指定します。ファイル形式に関しては、以下のようなフォーマットをしています。

::

    ModPara  modpara.def
    Trans    ztransfer.def
    InterAll zinterall.def
    OneBodyG zcisajs.def

| 

ファイル形式
^^^^^^^^^^^^

[keyword] [filename]

パラメータ
^^^^^^^^^^

-  :math:`[`\ keyword\ :math:`]`

   **形式 :** string型 (固定)

   **説明 :** キーワードを指定します。

-  :math:`[`\ filename\ :math:`]`

   **形式 :** string型

   **説明 :** キーワードにひも付けられるファイル名を指定します(任意)。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  キーワードを記載後、半角空白を開けた後にファイル名を書きます。ファイル名は自由に設定できます。

-  必ず指定しなければいけないパラメータはModParaです。

-  各キーワードは順不同に記述できます。

-  指定したキーワード、ファイルが存在しない場合はエラー終了します。

-  :math:`\#`\ で始まる行は読み飛ばされます。

-  ファイル読込用キーワードは :numref:`Table 4.2` により指定します。

.. _Table 4.2:
.. csv-table:: 定義ファイル一覧
    :header: "Keywords", "指定ファイルの概要"
    :widths: 4, 20

    "ModPara","計算で用いるパラメータの指定をします。"
    "Trans","一般的一体相互作用に関する設定をします。"
    "InterAll", "一般的二体相互作用に関する設定をします。"
    "CoulombIntra", "内部クーロン相互作用に関する設定をします。"
    "CoulombInter", "サイト間クーロン相互作用に関する設定をします。"
    "Hund", "フント結合に関する設定をします。"
    "PairHop", "ペアホッピングに関する設定をします。"
    "Exchange", "交換相互作用に関する設定をします。"
    "Ising", "イジング相互作用に関する設定をします。"
    "PairLift", "ペアリフト相互作用に関する設定をします。"
    "Initial", "入力する一体グリーン関数の設定をします。"
    "OneBodyG", "出力する一体グリーン関数 \ :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle` に関する設定をします。"

.. raw:: latex

   \newpage
