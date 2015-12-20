## Homework 3 ##

### Description of implementation approach and comments ###


* Windows 7 64-bits
* Comipler & IDE: Microsoft Visual Studio 2012
* Image View Software: IrfanView
* Intel(R) Core(TM) i5-3470 CPU @ 3.20GHz 3.20GHz
* Intel(R) HD Graphics 
* RAM 4.00 GB

#### Median Cut Algorithm ####

根據論文 [P. Debevec, A Median Cut Algorithm for Light Probe Sampling, SIGGRAPH 2006 Sketches and Applications.](http://dl.acm.org/citation.cfm?id=1187029) 中，預期要將 `pbrt/lights/infinite.cpp` 中的 `class InfiniteAreaLight` 用數個點光源取代 Infinite Area Light 的寫法，提升均勻取樣的效能，而 Median Cut Algorithm 在預處理非常快，根據用多少量的點光源將影響品質，若在品質不用太好的 rendering 環境下，這是一個不錯的降質提升速度的方案。

算法的步驟如下：

1. 將入射光場影像 (light probe image) 切成好幾個矩形區域，每一個區域將取用一個點光源代替。將入射光場影像轉換成灰階亮度 $Y$，如論文中所提及的方案 $Y = 0.2125 R + 0.7154 G + 0.0721 B$ 這類型的轉換。
2. 對於每一個區域將會根據最長維度切割成兩個子區域。切割成兩個子區域的亮度總和越接近越好。
3. 若切割區域數量不到目標的數量，則重複步驟 2。
4. 最後將每一個區域算出代表的點光源，並加總區域內的亮度和，隨後取樣根據亮度和作為取樣根據 (在 `Spectrum MedianCutEnvironmentLight::Sample_L(const Point&, float, const LightSample&, float, Vector&, float*, VisibilityTester*)` 中使用)，用每一個區域內部的 pixel 位置和亮度計算重心作為代表的點光源。

算法類似於 k-d Tree，但特別的是每次選擇區域維度最長的那一段進行切割，而不是像 k-d Tree 則是採用輪替維度進行。

Median Cut Algorithm 需要 $\mathcal{O}(HW)$ 時間 $\mathcal{O}(HW)$ 空間來預處理亮度資訊。若要切割成 $N$ 個區域，需要 $\mathcal{O}(\log N)$ 次迭代，每一次迭代增加兩倍區域數量。將一個區域切割時，針對維度最長的那一軸二分搜尋，二分搜尋計算其中一個區域的亮度和是否是整個區域亮度的一半，由於預處理完成的表格，可以在 $\mathcal{O}(1)$ 時間內計算任意平行兩軸的矩形內部元素總和。


維護 `sumTable[i][j]` 計算左上角 $(0, 0)$ 右下角 $(i, j)$ 的亮度和，計算任意平行兩軸矩形內部和只需要 $\mathcal{O}(1)$ 時間。

```c
float getEnergy(float sumTable[], int width, int height) {
	float e = sumTable[VERT(ru, rv)];
	if (lu)	e -= sumTable[VERT(lu-1, rv)];
	if (lv)	e -= sumTable[VERT(ru, lv-1)];
	if (lu && lv)	e += sumTable[VERT(lu-1, lv-1)];
	return e;
}
```


> 另一種方案，可以從 pbrt 原生的 `class MIPMap` 取代上述寫法，`MIPMap` 的寫法採用分層式的架構，每一層將原圖長寬各縮小一半。計算某矩形元素總和時，藉由兩層的總和估計進行內插。理論上，這種寫法雖然不夠精準，但提供很好的快取優勢。

### 重心計算 ####



