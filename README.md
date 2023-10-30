# SimpleContinuation
## 


## A módszerről röviden
Adott a paraméterfüggő 1D folytonos dinamikai rendszer:
```math
\dot{x} = f(x,p),
```
ahol x és p valós skalármennyiségek. 
A cél, hogy automatizáljuk a fixpontok bifurkáció-elemzését. Ha ismerjük egy adott $p_0$ paraméterre a fixpont helyzetét $x^{\*}(p_0)$, azaz
```math
f(x^{*}(p_0),p_0) = 0
```
a $(x^{\*}(p_0),p_0)$ megoldás felhasználható a fixpont helyzetének becslésére a $p_0 + \Delta p$ értékre. Gyakorlatilag "folytatni" szeretnénk a meglévő megoldást, vagyis kiszámítani egy teljes megoldás-görbét (ág,branch). Ezeknek a módszereknek a gyűjtőneve a numerikus kontinuáció (**numerical continuation**).

## Egyszerű közelítés
Numerikusan szeretnénk közelíteni a $x^{\*}(p)$ görbét egy adott $p \in [p_{min},p_{max}]$ intevallumon. 

Itt explicit görbeként kezeljük, ami csak bizonyos esetekben lehetséges (a módszer nem érvényes, ha a görbe "visszafordul").
A fixpont helyzetét megközelíthetjük a $p = p_0 + \Delta p$ paraméterértékre ($|\Delta p| << 1$) ha sorbafejtünk első rendig,
```math
x^{*}(p) \approx x_0^{*}(p_0) + \frac{d x^{*}}{dp} \Delta p 
```
Felhasználjuk a kezdeti megoldást:
```math
f(x^{*}(p_0),p_0) = 0
```
Differenciálva mindkét oldalt:
```math
\frac{\partial f}{\partial x} dx + \frac{\partial f}{\partial p} dp = 0 
```
vagy rövidebben: $f_x dx + f_p dp = 0$, ahonnan
```math
\frac{dx^*}{dp} = -\frac{f_p}{f_x}
```
Ezt visszaírva az $x^{*}(p)$-re kapott közelítésbe:
```math
 x^{*}(p) \approx x_0^{*}(p_0) - \left[ \frac{f_p}{f_x}  \right]_{x^{*}_0,p_0} \cdot \Delta p 
```
Ezt az összefüggést használhatjuk, hogy egy lépéssel folytassuk a görbét. A csomag `continuation` függvénye numerikusan közelíti az $f_x,f_p$ deriváltakat, így ezeket nem szükséges megadni. 

## Használat

