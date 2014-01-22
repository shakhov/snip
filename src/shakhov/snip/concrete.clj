
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
    
    :Ari
    (fnk
     [reinf]
     (map (fn [{:keys [n d rebar]}]
            (* n ((:A rebar) d)))
          (:bottom reinf)))
    
    :Rri
    (fnk
     [reinf]
     (map #(get-in % [:rebar :Rr]) (:bottom reinf)))
    
    :ari
    (fnk
     [reinf]
     (map :ar (:bottom reinf)))
    
    :Ar
    (fnk
     [Ari]
     (apply + Ari))
    
    :Arci
    (fnk
     [reinf]
     (map (fn [{:keys [n d rebar]}]
            (* n ((:A rebar) d)))
          (:top reinf)))
    
    :Rrci
    (fnk
     [reinf]
     (map #(get-in % [:rebar :Rr]) (:top reinf)))
    
    :arci
    (fnk
     [reinf]
     (map :ar (:top reinf)))
    
    :Arc
    (fnk
     [Arci]
     (apply + Arci))
    
    :Nr
    (fnk
     [Ari Rri]
     (apply + (map * Ari Rri)))
    
    :Nrc
    (fnk
      [Arci Rrci]
      (apply + (map * Arci Rrci)))
    
    :ar
    (fnk
     [Nr Ari Rri ari]
     (/ (apply + (map * Ari Rri ari))
        Nr))
    
    :arc
    (fnk
     [Nrc Arci Rrci arci]
     (if (zero? Nrc)
       (si/cm 0)
       (/ (apply + (map * Arci Rrci arci))
          Nrc)))
    
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
      [A-cr reinf]
      (let [r (:bottom reinf)
            d (map :d r)
            n (map :n r)
            beta-cr (map :beta-cr r)]
        (/ A-cr
           (apply + (map * beta-cr d n)))))
     
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
      [reinf]
      (let [rb (:bottom reinf)
            rb (map (fn [{:keys [rebar n d] :as row}]
                      (assoc row :Ar (* n ((:A rebar) d))))
                    rb)
            Ar-max (apply max (map :Ar rb))]
        (map (apply max-key :ar (filter #(>= (:Ar %) (* 0.5 Ar-max))
                                        rb))
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
      [reinf Er]
      (apply + (map (fn [{:keys [n d rebar]}]
                      (* n ((:A rebar) d)
                         (/ (:Er rebar)
                            Er)))
                    (:bottom reinf))))
     
     :Arc-red
     (fnk
      [reinf Er]
      (apply + (map (fn [{:keys [n d rebar]}]
                      (* n ((:A rebar) d)
                         (/ (:Er rebar)
                            Er)))
                    (:top reinf))))
     
     :ar-red
     (fnk
      [reinf Ar-red Er]
      (/ (apply + (map (fn [{:keys [n d rebar ar]}]
                         (* ar
                            n ((:A rebar) d)
                            (/ (:Er rebar)
                               Er)))
                       (:bottom reinf)))
         Ar-red))
     
     :arc-red
     (fnk
      [reinf Arc-red Er]
      (/ (apply + (map (fn [{:keys [n d rebar ar]}]
                         (* ar
                            n ((:A rebar) d)
                            (/ (:Er rebar)
                               Er)))
                       (:top reinf)))
         Arc-red))
     
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

  (def rect-bending
    (fnk
     {:keys [Rc b h reinf] :as args}
     (let [input (merge {:Ap ((pow si/m 2) 0) :Apc ((pow si/m 2) 0)
                         :Rp  (si/MPa 0) :Rpc (si/MPa 500)
                         :ap (si/m 0) :apc (si/m 0)
                         :sigma-pc1 (si/MPa 0)
                         :sigma-p (si/MPa 0)}
                        args)
           no-Arc  (lazy-rect (update-in input [:reinf :top]
                                         (fn [t] (mapv #(assoc % :n 0) t))))
           all-Arc (lazy-rect input)
           arc (:arc all-Arc)]
       (cond
        (<  (:x no-Arc) (* 2 arc)) (dissoc  no-Arc :M-max-sc)
        (>= (:x all-Arc)(* 2 arc)) (dissoc all-Arc :M-max-sc)
        :else (assoc (dissoc all-Arc :M-max)
                :x (* 2 arc))))))

  (def T-bending
    (fnk
     {:keys [Rc b h bf hf reinf] :as args}
     (let [input (merge {:Ap ((pow si/m 2) 0) :Apc ((pow si/m 2) 0)
                         :Rp  (si/MPa 0) :Rpc (si/MPa 500)
                         :ap (si/m 0) :apc (si/m 0)
                         :sigma-pc1 (si/MPa 0)
                         :sigma-p (si/MPa 0)}
                        args)
           no-Arc  (lazy-T (update-in input [:reinf :top]
                                         (fn [t] (mapv #(assoc % :n 0) t))))
           all-Arc (lazy-T input)
           arc (:arc all-Arc)]
       (cond
        (<  (:x no-Arc) (* 2 arc)) (dissoc  no-Arc :M-max-sc)
        (>= (:x all-Arc)(* 2 arc)) (dissoc all-Arc :M-max-sc)
        :else (assoc (dissoc all-Arc :M-max)
                :x (* 2 arc))))))

  (def bending
    (fnk
     {:keys [Rc b h reinf] :as args}
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
     {:keys [h b reinf Er Rc-mc2] :as args}
     (lazy-rect-cracking args)))

  (def T-cracking
    (fnk
     {:keys [h b hf bf reinf Er Rc-mc2] :as args}
     (lazy-T-cracking args)))

  (def cracking
    (fnk
     {:keys [h b reinf Er Rc-mc2] :as args}
     (if (and (:hf args)
              (:bf args))
       (let [as-rect (rect-cracking (assoc args :b (:bf args)))
             as-T    (T-cracking    args)]
         (if (<= (:x-el as-T) (:hf args))
           (assoc as-rect :as "[]")
           (assoc as-T :as "T")))
       (assoc (rect-cracking args) :as "[]")))))
