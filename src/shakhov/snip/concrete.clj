
(ns shakhov.snip.concrete
  (:refer-clojure :exclude [time force + - * / < > <= >= = zero? pos? neg? sgn abs
                            sin cos tan asin acos atan exp log min max])

  (:use [shakhov.flow.core]
        [shakhov.snip.utils])

  (:use [clojure.algo.generic.arithmetic :only [+ - * /]]
        [clojure.algo.generic.comparison :only [< > <= >= = zero? pos? neg? min max]]
        [clojure.algo.generic.math-functions :only [pow sqrt sgn abs sin cos tan
                                                    asin acos atan exp log]])
  (:require [shakhov.snip.dimensions :as dim]
            [shakhov.snip.units :as si]))

(def xi-flow
  (flow
   {
    :xi<=xi-y
    (fnk
     ""
     [xi xi-y]
     (<= xi xi-y))
    
    :xi
    (fnk
     ""
     [h0 x]
     (/ x h0))
    
    :xi-y
    (fnk
     ""
     [omega sigma-1 sigma-2]
     (/ omega
        (+ 1 (* (/ sigma-1 sigma-2)
                (- 1 (/ omega 1.1))))))
    
    :omega
    (fnk
     ""
     [Rb]
     (- 0.85 (* ((/ si/MPa) 0.008) Rb)))
    
    :sigma-1
    (fnk
     ""
     [Rs Rp sigma-p]
     (if (zero? sigma-p)
       Rs
       (+ Rp (si/MPa 500) (- sigma-p))))
    
    :sigma-2
    (fnk
     ""
     []
     (si/MPa 500))
    
    }))
(def rect-bending-flow
  (flow
   {
    :M-max
    (fnk
     ""
     [Rb b x h0 h01 Asc-ef Rsc asc sigma-pc Apc apc]
     (+ (* Rb b x
           (- h0 (* 0.5 x)))
        (* Rsc Asc-ef
           (- h01 asc))
        (* sigma-pc Apc
           (- h0 apc))))
    
    :x
    (fnk
     ""
     [Rb b As Asc-ef Ap Apc Rs Rp Rsc sigma-pc]
     (/ (+ (* Rp Ap)
           (* Rs As)
           (* -1 Rsc Asc-ef)
           (* -1 sigma-pc Apc))
        (* Rb b)))
    
    :h01
    (fnk
     [h as]
     (- h as))
    
    :h0
    (fnk
     [h h01 ap As Rs Ap Rp]
     (if (zero? Ap)
       h01
       (/ (+ (* As Rs h01)
             (* Ap Rp (- h ap)))
          (+ (* As Rs)
             (* Ap Rp)))))
    
    :sigma-pc
    (fnk
     ""
     [Rpc sigma-pc1]
     (if (<= sigma-pc1 Rpc)
       (si/MPa 0.0)
       (- Rpc sigma-pc1)))
    
    :M-max-sc
    (fnk
     ""
     [Rp Ap Rs As h0 asc]
     (* (+ (* Rp Ap)
           (* Rs As))
        (- h0 asc)))
    }))

(let [lazy-rect (lazy-compile rect-bending-flow)
      lazy-xi   (lazy-compile xi-flow)]
  (def rect-bending
    (fnk
     {:keys [Rb b h Rs As as] :as args}
     (let [input (merge {:Rsc Rs :Asc ((pow si/m 2) 0) :asc (si/m 0)
                         :Ap ((pow si/m 2) 0) :Apc ((pow si/m 2) 0)
                         :Rp  (si/MPa 0) :Rpc (si/MPa 500)
                         :ap (si/m 0) :apc (si/m 0)
                         :sigma-pc1 (si/MPa 0)
                         :sigma-p (si/MPa 0)}
                        args)
           no-Asc  (lazy-rect (assoc input :Asc-ef ((pow si/m 2) 0)))
           all-Asc (lazy-rect (assoc input :Asc-ef (:Asc input)))]
       (lazy-xi (cond
                 (<  (:x no-Asc) (* 2 (:asc input))) (dissoc  no-Asc :M-max-sc)
                 (>= (:x all-Asc)(* 2 (:asc input))) (dissoc all-Asc :M-max-sc)
                 :else (assoc (dissoc all-Asc :M-max)
                              :x (* 2 (:asc input)))))))))
