(ns shakhov.snip.utils)

(defn clip
  [{xmin :min xmax :max} x]
  (let [x (if xmin (max xmin x) x)
        x (if xmax (min xmax x) x)]
    x))

(defn table-2d
  [{:keys [xp yp data clip]}]
  (fn [x y]
    (let [i1 (dec (count (take-while #(< % x) xp)))
          i2 (inc i1)
          i1 (if (= i1 -1) 
               (if (:l clip) 0 1)
               i1)
          i2 (if (= i2 (count xp)) 
               (if (:r clip) (- i2 1) (- i2 2)) 
               i2)
          j1 (dec (count (take-while #(< % y) yp)))
          j2 (inc j1)
          j1 (if (= j1 -1) 
               (if (:t clip) 0 1) j1)
          j2 (if (= j2 (count yp)) 
               (if (:b clip) (- j2 1) (- j2 2)) 
               j2)
          x1 (get xp i1)
          x2 (get xp i2)
          y1 (get yp j1)
          y2 (get yp j2)
          dx (- x2 x1)
          dy (- y2 y1)
          dxdy (* dx dy)
          d11 (get-in data [j1 i1])
          d21 (get-in data [j1 i2])
          d12 (get-in data [j2 i1])
          d22 (get-in data [j2 i2])]
      (if (zero? dxdy)
        (cond (and (zero? dx) (zero? dy)) d11
              (zero? dx) (+ (* d11 (/ (- y2 y ) dy))
                            (* d12 (/ (- y  y1) dy)))
              (zero? dy) (+ (* d11 (/ (- x2 x ) dx))
                            (* d21 (/ (- x  x1) dx))))
        (+ (* d11 (/ (* (- x2 x ) (- y2 y )) dxdy))
           (* d21 (/ (* (- x  x1) (- y2 y )) dxdy))
           (* d12 (/ (* (- x2 x ) (- y  y1)) dxdy))
           (* d22 (/ (* (- x  x1) (- y  y1)) dxdy)))))))