
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
     [Rc b x h0 h01 Nrc arc sigma-pc Apc apc]
     (+ (* Rc b x
           (- h0 (* 0.5 x)))
        (* Nrc
           (- h01 arc))
        (* sigma-pc Apc
           (- h0 apc))))
    
    :x
    (fnk
     ""
     [Rc b Ap Apc Rp sigma-pc Nr Nrc]
     (/ (+ (* Rp Ap)
           Nr
           (* -1 Nrc)
           (* -1 sigma-pc Apc))
        (* Rc b)))
    
    :Ar
    (fnk
     [bottom-rebar]
     (apply + (pow (si/m 0) 2)
            (map (fnk
                  {:keys [n d] {:keys [A]} :rebar}
                  (* n (A d)))
                 bottom-rebar)))
    
    :Arc
    (fnk
     [top-rebar]
     (apply + (pow (si/m 0) 2)
            (map (fnk
                  {:keys [n d] {:keys [A]} :rebar}
                  (* n (A d)))
                 top-rebar)))
    
    :Nr
    (fnk
     [bottom-rebar]
     (apply + (si/N 0)
            (map (fnk
                  {:keys [n d] {:keys [A Rr]} :rebar}
                  (* n Rr (A d)))
                 bottom-rebar)))
    
    :Nrc
    (fnk
     [top-rebar]
     (apply + (si/N 0)
            (map (fnk
                  {:keys [n d] {:keys [A Rr]} :rebar}
                  (* n Rr (A d)))
                 top-rebar)))
    
    :ar
    (fnk
     [Nr bottom-rebar]
     (if (zero? Nr)
       (si/m 0)
       (/ (apply + (* si/N si/m)
                 (map (fnk
                       {:keys [ar n d] {:keys [A Rr]} :rebar}
                       (* ar n Rr (A d)))
                      bottom-rebar))
          Nr)))
    
    :arc
    (fnk
     [Nr top-rebar]
     (if (zero? Nr)
       (si/m 0)
       (/ (apply + (* si/N si/m)
                 (map (fnk
                       {:keys [ar n d] {:keys [A Rr]} :rebar}
                       (* ar n Rr (A d)))
                     top-rebar))
          Nr)))
    
    :h01
    (fnk
     [h ar]
     (- h ar))
    
    :h0
    (fnk
     [h h01 ap Nr Ap Rp]
     (if (zero? Ap)
       h01
       (/ (+ (* Nr h01)
             (* Ap Rp (- h ap)))
          (+ Nr
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
     [Rp Ap Nr h0 arc]
     (* (+ Nr (* Rp Ap))
        (- h0 arc)))
    
    :Ap  (fnk [] ((pow si/m 2) 0))
    :Apc (fnk [] ((pow si/m 2) 0))
    
    :Rp (fnk [] (si/MPa 0))
    :Rpc (fnk [] (si/MPa 500))
    
    :ap (fnk [] (si/m 0))
    :apc (fnk [] (si/m 0))
    
    :sigma-pc1 (fnk [] (si/MPa 0))
    :sigma-p (fnk [] (si/MPa 0))
    
    :top-rebar (fnk [] [])
    }))
(def T-bending-flow
  (merge rect-bending-flow
         (flow
          {
            :M-max
           (fnk
            ""
            [Rc b x bf hf h0 h01 Nrc arc sigma-pc Apc apc]
            (+ (* Rc b x
                  (- h0 (* 0.5 x)))
               (* Rc (- bf b) hf
                  (- h0 (* 0.5 hf)))
               (* Nrc
                  (- h01 arc))
               (* sigma-pc Apc
                  (- h0 apc))))
           
           :x
           (fnk
            ""
            [Rc b bf hf Ap Apc Rp sigma-pc Nr Nrc]
            (/ (+ (* Rp Ap)
                  Nr
                  (* -1 Nrc)
                  (* -1 sigma-pc Apc)
                  (* -1 Rc hf (- bf b)))
               (* Rc b)))
           })))
(def rect-crack-width-flow
  (merge
   (flow
    {
     :a-cr
     (fnk
      [sigma-r Er psi-cr]
      (* psi-cr (/ sigma-r Er)))
     
     :M-max-cr
     (fnk
      [delta-cr Er I-red-el n' Zr psi-cr]
      (/ (* delta-cr Er I-red-el)
         (* n' Zr psi-cr)))
     
     :psi-cr
     (fnk
      [R-cr]
      (si/cm (* 1.5 (sqrt (:magnitude (si/cm R-cr))))))
     
     :R-cr
     (fnk
      [A-cr bottom-rebar]
      (/ A-cr
         (apply + (si/m 0)
                (map (fnk
                      {:keys [d n beta-cr]}
                      (* beta-cr d n))
                     bottom-rebar))))
     
     :A-cr
     (fnk
      [b hr]
      (* b hr))
     
     :hr
     (fnk
      [ar-cr h x-el d-cr]
      (min (- h x-el)
           (+ ar-cr (* 6 d-cr))))
     
     [ar-cr d-cr]
     (fnk
      [bottom-rebar]
      (let [Ari (map (fnk
                      {:keys [n d] {:keys [A]} :rebar :as row}
                      (* n (A d)))
                     bottom-rebar)
            Ar-max (apply max Ari)
            bottom-rebar (map #(assoc %1 :Ar %2)
                              bottom-rebar Ari)]
        (map (reduce #(if (> (:ar %1) (:ar %2)) %1 %2)
                     (filter #(>= (:Ar %) (* 0.5 Ar-max))
                             bottom-rebar))
             [:ar :d])))
     
     :sigma-r
     (fnk
      [M-ser I-red-el Zr n']
      (* n' (/ M-ser I-red-el)
         Zr))
     
     :x-el
     (fnk
      [n' b h Ar-red Arc-red ar-red arc-red]
      (let [A (* 1/2 b)
            B (* n' (+ Ar-red Arc-red))
            C (* n' (- (* Ar-red ar-red)
                       (* Ar-red h)
                       (* Arc-red arc-red)))
            D (- (pow B 2) (* 4 A C))]
        (/ (+ (- B) (sqrt D))
           2 A)))
     
     :Zr
     (fnk
      [x-el h ar-red]
      (- h x-el ar-red))
     
     :I-red-el
     (fnk
      [Ar-red Arc-red n' b h x-el ar-red arc-red]
      (+ (* 1/3 b (pow x-el 3))
         (* n'
            (+ (* Arc-red (pow (- x-el arc-red)  2))
               (* Ar-red  (pow (- h x-el ar-red) 2))))))
     
     :Ar-red
     (fnk
      [bottom-rebar Er]
      (apply + (pow (si/m 0) 2)
             (map (fnk
                   {:keys [n d] {A :A Er' :Er} :rebar}
                   (* n (A d)
                       (/ Er' Er)))
                  bottom-rebar)))
     
     :Arc-red
     (fnk
      [top-rebar Er]
      (apply + (pow (si/m 0) 2)
             (map (fnk
                   {:keys [n d] {A :A Er' :Er} :rebar}
                   (* n (A d)
                       (/ Er' Er)))
                  top-rebar)))
     
     :ar-red
     (fnk
      [bottom-rebar Ar-red Er]
      (if (zero? Ar-red)
        (si/m 0)
        (/ (apply + (pow (si/m 0) 3)
                  (map (fnk
                        {:keys [n d ar] {A :A Er' :Er} :rebar}
                        (* ar n (A d)
                           (/ Er' Er)))
                       bottom-rebar))
           Ar-red)))
     
     :arc-red
     (fnk
      [top-rebar Arc-red Er]
      (if (zero? Arc-red)
        (si/m 0)
        (/ (apply + (pow (si/m 0) 3)
                  (map (fnk
                        {:keys [n d ar] {A :A Er' :Er} :rebar}
                        (* ar n (A d)
                           (/ Er' Er)))
                       top-rebar))
           Arc-red)))
     
     :sigma-c
     (fnk
      [I-red-el x-el M-ser]
      (* (/ M-ser I-red-el) x-el))
     
     :M-max-mc2
     (fnk
      [I-red-el x-el Rc-mc2]
      (/ (* Rc-mc2 I-red-el)
         x-el))
     })
   (select-keys rect-bending-flow [:Ar :Arc])))
(def T-crack-width-flow
  (merge
   rect-crack-width-flow
   (flow
    {
     :x-el
     (fnk
      [n' b h bf hf Ar-red Arc-red ar-red arc-red]
      (let [A (* 1/2 b)
            B (+ (* n' (+ Ar-red Arc-red))
                 (* hf (- bf b)))
            C (+ (* n' (- (* Ar-red ar-red)
                          (* Ar-red h)
                          (* Arc-red arc-red)))
                 (* -1/2 hf hf (- bf b)))
            D (- (pow B 2) (* 4 A C))]
        (/ (+ (- B) (sqrt D))
           2 A)))
     
     :Zr
     (fnk
      [x-el h ar-red]
      (- h x-el ar-red))
     
     :I-red-el
     (fnk
      [Ar-red Arc-red n' b h hf bf x-el ar-red arc-red]
      (+ (* 1/3 b (pow x-el 3))
         (* 1/12 (- bf b) (pow hf 3))
         (* hf (- bf b) (pow (- x-el (* 0.5 hf)) 2))
         (* n'
            (+ (* Arc-red (pow (- x-el arc-red)  2))
               (* Ar-red  (pow (- h x-el ar-red) 2))))))
     })))

(let [lazy-rect (lazy-compile (merge xi-flow rect-bending-flow))
      lazy-T    (lazy-compile (merge xi-flow T-bending-flow))]

  (defn cs-bending
    [flow input]
    (let [no-Arc  (flow (update-in input [:top-rebar]
                                   (fn [t] (mapv #(assoc % :n 0) t))))
          all-Arc (flow input)
          arc (:arc all-Arc)]
      (cond
       (<  (:x no-Arc) (* 2 arc)) (dissoc  no-Arc :M-max-sc)
       (>= (:x all-Arc)(* 2 arc)) (dissoc all-Arc :M-max-sc)
       :else (assoc (dissoc all-Arc :M-max)
               :x (* 2 arc)))))

  (def rect-bending
    (fnk
     {:keys [Rc b h bottom-rebar] :as args}
     (cs-bending lazy-rect args)))

  (def T-bending
    (fnk
     {:keys [Rc b h bf hf bottom-rebar] :as args}
     (cs-bending lazy-T args)))

  (def bending
    (fnk
     {:keys [Rc b h bottom-rebar] :as args}
     (if (and (:hf args)
              (:bf args))
       (let [as-rect (rect-bending (assoc args :b (:bf args)))
             as-T    (T-bending    args)]
         (if (<= (:x as-rect) (:hf args))
           (assoc as-rect :as "[]")
           (assoc as-T :as "T")))
       (assoc (rect-bending args) :as "[]")))))

(let [lazy-rect-cracking (lazy-compile rect-crack-width-flow)
      lazy-T-cracking    (lazy-compile T-crack-width-flow)]

  (def rect-cracking
    (fnk
     {:keys [h b bottom-rebar Er Rc-mc2] :as args}
     (lazy-rect-cracking args)))

  (def T-cracking
    (fnk
     {:keys [h b hf bf bottom-rebar Er Rc-mc2] :as args}
     (lazy-T-cracking args)))

  (def cracking
    (fnk
     {:keys [h b bottom-rebar Er Rc-mc2] :as args}
     (if (and (:hf args)
              (:bf args))
       (let [as-rect (rect-cracking (assoc args :b (:bf args)))
             as-T    (T-cracking    args)]
         (if (<= (:x-el as-T) (:hf args))
           (assoc as-rect :as "[]")
           (assoc as-T :as "T")))
       (assoc (rect-cracking args) :as "[]")))))
