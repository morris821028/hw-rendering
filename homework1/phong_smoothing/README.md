## 實驗操作 ##

### debug mode ###

* IDE : Microsoft Visual Studio 2012

若不在 debug mode 下，command line 直接 `$ ./pbrt.exe heightfields2_input/landsea-1.pbrt` 有可能會發生兩件事情，很難抓到 bug 所在之處。

* 連 create shape 都失敗，意即進度條都沒有跑出來，根本看不到物件被創建，程序就自動關閉。
* 進度條跑到一半，直接關閉程序，圖片檔案並沒有輸出成功。(通常是計算中發生除零，或者是數值過小，這次實驗則會是正規化向量時，發生零向量的計算比較多。)

偷懶一點，直接再 VS 2012 上面使用命令引數

* `工具列 > 偵錯(D) > pbrt 專案 ... > 組態(C) Debug > 命令引數`
* 打上絕對位置 `c:\Users\morris821028\Desktop\pbrt-v2\bin\heightfields2_input\hftest.pbrt`，預設是在 `pbrt-v2\src\pbrt.vs2012\pbrt.sln` 的路徑，相對路徑就看著辦。
* 產出來的圖片在  `pbrt-v2\src\pbrt.vs2012\xxx.exr`
* 測試都沒問題，才去編譯成 release 的版本進行測試，debug mode 編出來的速度慢上好幾倍。


## 實驗結果 ##

### Machine ###

* Windows 7 64-bits
* Intel(R) Core(TM) i5-3470 CPU @ 3.20GHz 3.20GHz
* RAM 4.00 GB
* Intel(R) HD Graphics (內顯)

### Calcuate Time ###

```
$ time ./pbrt.exe heightfields_input/hftest.pbrt
```

### Default Triangle ###

| \                | Origin(without Phong) | My Implementation |
|------------------|-----------------------|-------------------|
| hftest.pbrt      |                 125 ms|             130 ms|
| landsea-0.pbrt   |                 815 ms|            1050 ms|
| landsea-1.pbrt   |                 850 ms|            1020 ms|
| landsea-2.pbrt   |                 780 ms|             870 ms|
| landsea-big.pbrt |                6200 ms|            2200 ms|
| texture.pbrt     |                 450 ms|             520 ms|

### Phong Smoothing ###

| \                | Origin(without Phong) | My Implementation |
|------------------|-----------------------|-------------------|
| hftest.pbrt      |                 125 ms|             135 ms|
| landsea-0.pbrt   |                 815 ms|            1100 ms|
| landsea-1.pbrt   |                 850 ms|            1050 ms|
| landsea-2.pbrt   |                 780 ms|             900 ms|
| landsea-big.pbrt |                6200 ms|            2300 ms|
| texture.pbrt     |                 450 ms|             530 ms|
