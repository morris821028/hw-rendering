## Homework 2 ##

### Description of implementation approach and comments ###


* Windows 7 64-bits
* Comipler & IDE: Microsoft Visual Studio 2012
* Image View Software: IrfanView
* Intel(R) Core(TM) i5-3470 CPU @ 3.20GHz 3.20GHz
* Intel(R) HD Graphics 
* RAM 4.00 GB

### Ray Sphere Intersection ###

計算射線和球體的交點，可以參照 `/pbrt-v2/shapes/sphere.cpp` 中 `bool Shpere::Intersect()` 的算法。

假設球體圓心座標 $O$，射線單位向量 $I$ 的起點座標 $C$，且最近目標交點座標 $P$，原半徑 $\mathrm{radius}$。射線走單位時間 $t$ 會到達球面上。

* $\overrightarrow{OC} + \overrightarrow{I} \times t = \overrightarrow{OP}$
* $|\overrightarrow{OP}| = \text{radius}$
* $|\overrightarrow{OC} + \overrightarrow{I} \times t| = |\overrightarrow{OP}|$
* $|\overrightarrow{I}|^2 t^2 + 2 \; (\overrightarrow{OC} \cdot \overrightarrow{I}) \; t + |\overrightarrow{OC}|^2 - \text{radius}^2 = 0$

解一元二次方程式之後可以得到 $t$ 得值，並得到交點座標 $P$。

### Snell's Law ###

根據物理法則斯乃爾定律計算折射的方向向量，課程中提供三種做法

* Whitted's Method
* Heckbert's Method
* Other Method

其中以 Heckbert's Method 消耗最少計算數。

