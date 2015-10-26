##  測試暴力法  ##

### 安裝環境 ###

* Windows 7
* IDE Microsoft Visual Studio 2012
* Image View Software: [IrfanView](http://www.irfanview.com/) and extra plgin for `.exr`  [iv_formats.zip](http://www.irfanview.com/plugins.htm)

### 程式碼探索 ###

一開始為了測試，隨便亂來，只要稍微顯示即可，交點或轉換座標無所謂。想說只寫以下函數

```
bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
        DifferentialGeometry *dg) const {
	...
}
```

然後讓交點恆真，接著測試 `hftest.exr`

```
bool Heightfield2::IntersectP(const Ray &r) const {
	return true;
}
```

上述的寫法會造成一片黑，因為對於每一個交點會想辦法連到光源，一旦有一個 object 恆存在交點，會造成沒有一條連到光源，因此產生的所有圖片都是黑的。