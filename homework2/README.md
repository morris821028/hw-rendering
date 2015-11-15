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

### RasterToCamera ###

這部分處理後得到 `Transform RasterToCamera`。若計算錯誤，會造成一片黑或者圖片顯示的大小問題。座標轉換處理細節可以參考實際的例子，如下圖所示：

![http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays](submission/images/cambasic1A.png)

### Ray Weight ###

作業要求經過 `float GenerateRay()` 回傳 $\mathrm{weight} = \frac{\cos^4 \theta'}{Z^2}$，這麼設置會過暗，離 camera 最遠的透鏡半徑 $r$，則根據實驗需要多呈上一個常數 $r^2 \pi$，其原因猜測在於 ray 經過多次反射後 weight 若小於某一定值後會停止反射，所導致 weight 過小而黑幕情況。在 $Z$ 值部分，採用在最後一個透鏡 (離 camera 最遠) 得表面得到與 camera 的 z 軸距離 $Z$ 以及從最後一個透鏡折射出去的向量與 $\mathrm{Vector}(0, 0, 1)$ 的夾角 $\theta'$ 較為貼近期望所需。