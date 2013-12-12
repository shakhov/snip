
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
     [Rc]
     (- 0.85 (* ((/ si/MPa) 0.008) Rc)))
    
    :sigma-1
    (fnk
     ""
     [Rr Rp sigma-p]
     (if (zero? sigma-p)
       Rr
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
     [Rc b x h0 h01 Arc-ef Rrc arc sigma-pc Apc apc]
     (+ (* Rc b x
           (- h0 (* 0.5 x)))
        (* Rrc Arc-ef
           (- h01 arc))
        (* sigma-pc Apc
           (- h0 apc))))
    
    :x
    (fnk
     ""
     [Rc b Ar Arc-ef Ap Apc Rr Rp Rrc sigma-pc]
     (/ (+ (* Rp Ap)
           (* Rr Ar)
           (* -1 Rrc Arc-ef)
           (* -1 sigma-pc Apc))
        (* Rc b)))
    
    :h01
    (fnk
     [h ar]
     (- h ar))
    
    :h0
    (fnk
     [h h01 ap Ar Rr Ap Rp]
     (if (zero? Ap)
       h01
       (/ (+ (* Ar Rr h01)
             (* Ap Rp (- h ap)))
          (+ (* Ar Rr)
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
     [Rp Ap Rr Ar h0 arc]
     (* (+ (* Rp Ap)
           (* Rr Ar))
        (- h0 arc)))
    }))
(def rect-crack-width-flow
  (flow
   {
    :a-cr
    (fnk
     [sigma-r Er psi-cr]
     (* psi-cr (/ sigma-r Er)))
    
    :psi-cr
    (fnk
     [R-cr]
     (si/cm (* 1.5 (sqrt (:magnitude (si/cm R-cr))))))
    
    :R-cr
    (fnk
     [A-cr beta-cr nd d]
     (/ A-cr
        beta-cr nd d))
    
    :A-cr
    (fnk
     [b hr]
     (* b hr))
    
    :hr
    (fnk
     [ar-cr h x-el d]
     (min (- h x-el)
          (+ ar-cr (* 6 d))))
    
    :ar-cr
    (fnk
     [ar]
     ar)
    
    :sigma-r
    (fnk
     [M-ser I-red-el Zr n']
     (* n' (/ M-ser I-red-el)
        Zr))
    
    :x-el
    (fnk
     [n' b h ar arc Ar Arc]
     (let [a (* 1/2 b)
           b (* n' (+ Arc Ar))
           c (* n' (- (* Ar ar)
                             (* Ar h)
                             (* Arc arc)))
           D (- (pow b 2) (* 4 a c))]
       (/ (+ (- b) (sqrt D))
          2 a)))
    
    :Zr
    (fnk
     [x-el h ar]
     (- h x-el ar))
    
    :I-red-el
    (fnk
     [Ar Arc n' h x-el ar arc b]
     (+ (* 1/3 b (pow x-el 3))
        (* n'
           (+ (* Arc (pow (- x-el arc)  2))
              (* Ar  (pow (- h x-el ar) 2))))))
    
    :sigma-b
    (fnk
     [I-red-el x-el M-ser]
     (* (/ M-ser I-red-el) x-el))
    }))

(let [lazy-rect (lazy-compile rect-bending-flow)
      lazy-xi   (lazy-compile xi-flow)]
  (def rect-bending
    (fnk
     {:keys [Rc b h Rr Ar ar] :as args}
     (let [input (merge {:Rrc Rr :Arc ((pow si/m 2) 0) :arc (si/m 0)
                         :Ap ((pow si/m 2) 0) :Apc ((pow si/m 2) 0)
                         :Rp  (si/MPa 0) :Rpc (si/MPa 500)
                         :ap (si/m 0) :apc (si/m 0)
                         :sigma-pc1 (si/MPa 0)
                         :sigma-p (si/MPa 0)}
                        args)
           no-Arc  (lazy-rect (assoc input :Arc-ef ((pow si/m 2) 0)))
           all-Arc (lazy-rect (assoc input :Arc-ef (:Arc input)))]
       (lazy-xi (cond
                 (<  (:x no-Arc) (* 2 (:arc input))) (dissoc  no-Arc :M-max-sc)
                 (>= (:x all-Arc)(* 2 (:arc input))) (dissoc all-Arc :M-max-sc)
                 :else (assoc (dissoc all-Arc :M-max)
                              :x (* 2 (:arc input)))))))))
