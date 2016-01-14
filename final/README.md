## Final Project ##

### Description of implementation approach and comments ###


* Windows 7 64-bits
* Comipler & IDE: Microsoft Visual Studio 2012
* Image View Software: IrfanView
* Intel(R) Core(TM) i5-3470 CPU @ 3.20GHz 3.20GHz
* Intel(R) HD Graphics 
* RAM 4.00 GB

從一般 BVH 架構中，一般都是用 full binary tree，子節點要不 0 個要不 2 個。若有 $N$ 個 primitive 物件，則表示有 $N$ 個葉節點放置這些 primitive 物件，而有 $N-1$ 個內部節點紀錄 Bounding Box 的部分。在測試交點和遍歷走訪的使用上，最慘有一半是多餘計算和走訪，而另一半屬於加速結構。

在論文 [Ray Specialized Contraction on Bounding Volume Hierarchies]() 中，企圖想要利用 generic tree 取代一般使用 full binary tree 實作，在不更動 BVH 的效能下，減少運行上較沒用的內部節點，用以提升遍歷走訪效能，以及提升內存快取的效率。

降低內部節點步驟如下所示：

1. 利用原生方法建立 full binary tree 的 BVH (利用各種分割策略完成)
2. 進行坍倒 (flatten)，將二元樹不連續的記憶體分布調整成線性分布，加速遍歷走訪的內存快取效率。
3. 靜態調整成 generic tree 版本，藉由啟發式收縮 (Contract)，利用節點與其子節點的 Bounding Box 表面積比例，評估浪費的走訪節點。
4. 動態調整部分，採用隨機取樣，根據取樣 ray，取樣走訪次數，將比較容易打到的節點盡可能收縮到接近根節點。

### 實作部分 ###

從實作中，在步驟 2. 約略可以減少 $25\%$ 的節點，在記憶體方面的影響量沒有太大影響，節點紀錄資訊也增加 (`sizeof(struct Node)` 相對地大上許多)。

在步驟 3. 中，根據 pbrt-v2 的架構，加速結構能取得的場景資訊並不容易竄改，大部分的類別函數都是 `const function()`，意即無法變動 object member 的值，但針對指向位址的值可以改變。這類寫法，猜想所有加速結構都是靜態決定，在多核心運行時較不會存在同步 overhead 的影響。

在此，換成另外一種修改方案，在 `pbrt-v2/core/scene.h` 的 `bool Scene::Intersect(...)` 函數中加入 `aggregate->tick();`。利用 `aggregate->tick();` 這個函數，大部分呼叫都沒有更動樹狀結構。當呼叫次數 $T$ 達到一定次數時，加速結構會進行大規模的結構調整。根據 pbrt rendering 的步驟，儘管不斷地測試或者是估計 $T$ 要設定的值，都無法接近全局的取樣評估，其中最大的原因是取樣順序和局部取樣調整，從理論上得知不可能會有比較好的結果。

修改檔案如下：

#### 修改檔案清單 ####
```
.
├── accelerators
│   ├── bvhcontract.cpp
│   └── bvhcontract.h
└── core
    ├── api.cpp
    ├── primitive.cpp
    ├── primitive.h
    └── scene.h
```

#### `core/api.cpp` ####

```cpp
Primitive *MakeAccelerator(const string &name,
                           const vector<Reference<Primitive> > &prims,
                           const ParamSet &paramSet) {
	...
	else if (name == "bvhcontract")
		accel = CreateBVHContractAccelerator(prims, paramSet);
	...
}
```

#### `core/primitive.h` ####

```cpp
class Primitive : public ReferenceCounted {
public:
	...
	// MORRIS ADD
	virtual void tick();
	...
protected:
	...
};
```

#### `core/scene.h` ####
```cpp
class Scene {
public:
    // Scene Public Methods
	bool Intersect(const Ray &ray, Intersection *isect) const {
		...
		aggregate->tick();
		...
	}
};
```

### Test sences ###

* [pbrt.org scenes](http://www.pbrt.org/scenes.php)

### Reference paper ###

[Ray Specialized Contraction on Bounding Volume Hierarchies, Yan Gu Yong He Guy E. Blelloch](http://www.cs.cmu.edu/~ygu1/paper/PG15/conference.pdf)