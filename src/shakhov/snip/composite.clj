(ns shakhov.snip.composite
  (:refer-clojure :exclude [time force + - * / < > <= >= = zero? pos? neg? sgn abs
                            sin cos tan asin acos atan exp log min max])

  (:use [shakhov.flow.core]
        [shakhov.snip.utils])

  (:use [clojure.algo.generic.arithmetic :only [+ - * /]]
        [clojure.algo.generic.comparison :only [< > <= >= = zero? pos? neg? min max]]
        [clojure.algo.generic.math-functions :only [pow sqrt sgn abs sin cos tan
                                                    asin acos atan exp log]])
  (:require [shakhov.snip.dimensions :as dim]
            [shakhov.snip.units :as si]
            [shakhov.snip.steel :as steel]))

;;
;;  Utilities
;;

(def ^:private I-rect
  (fnk {:keys [b h dz] :or {dz nil}}
       (let [dz (or dz (h 0.0))]
         (+ (* 1/12 b h h h)
            (* b h dz dz)))))

;;
;; Steel concrete composite

;;

(def I-cs-geometry
  (flow
   {:nc
    (fnk
     [Est Ec]
     (/ Est Ec))

    :nr
    (fnk
     [Est Er]
     (/ Est Er))

    :Ec-creep
    (fnk
     [Ec phi-creep]
     (/ Ec (+ 1 phi-creep)))

    :Ec-shr
    (fnk
     [Ec]
     (* 0.5 Ec))

    :nc-creep
    (fnk
     [Est Ec-creep]
     (/ Est Ec-creep))

    :nc-shr
    (fnk
     [Est Ec-shr]
     (/ Est Ec-shr))

    :psi-crack
    (fnk [] 0.5)

    ;; Ai

    :Atf-i
    (fnk
     [top-flange]
     (map (fn [{:keys [w t]}] (* w t))
          top-flange))
    :Aw
    (fnk
     {{:keys [h t]} :web}
     (* h t))

    :Abf-i
    (fnk
     [bottom-flange]
     (map (fn [{:keys [w t]}] (* w t))
          bottom-flange))

    :Atf
    (fnk
     [Atf-i]
     (apply + Atf-i))

    :Abf
    (fnk
     [Abf-i]
     (apply + Abf-i))

    :Asl
    (fnk
     [slab Ar]
     (- (* (:b slab)
           (:h slab))
        Ar))

    :Asl-red
    (fnk
     [Asl nc]
     (/ Asl nc))

    :Asl-red-creep
    (fnk
     [Asl nc-creep]
     (/ Asl nc-creep))

    :Asl-red-shr
    (fnk
     [Asl nc-shr]
     (/ Asl nc-shr))

    :Ar
    (fnk {{:keys [n d]} :reinf}
         (* n 1/4 Math/PI d d))

    :Ar-red
    (fnk
     [Ar nr]
     (/ Ar nr))

    :Ar-crack
    (fnk
     [Ar-red psi-crack]
     (/ Ar-red psi-crack))

    ;; A

    :Ast
    (fnk
     [Atf Abf Aw]
     (+ Atf Aw Abf))

    :Astr
    (fnk
     [Ast Ar-red]
     (+ Ast Ar-red))

    :Astc
    (fnk
     [Asl-red Ast Ar-red]
     (+ Ast Asl-red Ar-red))

    :Astc-creep
    (fnk
     [Asl-red-creep Ast Ar-red]
     (+ Ast Asl-red-creep Ar-red))

    :Astc-shr
    (fnk
     [Asl-red-shr Ast Ar-red]
     (+ Ast Asl-red-shr Ar-red))

    :Astc-crack
    (fnk
     [Ast Ar-crack]
     (+ Ast Ar-crack))

    ;; Zw

    :Zw-tf
    (fnk
     [top-flange web]
     (let [z (reverse (reduce (fn [z {:keys [t]}]
                                (conj z (- (peek z) t)))
                              [(- (/ (:h web) 2))]
                              (reverse top-flange)))]
       (map #(/ (+ %1 %2) 2) z (next z))))

    :Zw-bf
    (fnk
     [bottom-flange web]
     (let [z (reduce (fn [z {:keys [t]}]
                       (conj z (+ (peek z) t)))
                     [(/ (:h web) 2)]
                     bottom-flange)]
       (map #(/ (+ %1 %2) 2) z (next z))))

    :Zw-sl
    (fnk
     [slab web top-flange]
     (- (+ (/ (:h web) 2)
           (apply + (map :t top-flange))
           (/ (:h slab) 2))))

    :Zw-r
    (fnk
     [Zw-sl]
     Zw-sl)

    ;; Sw

    :Sw-st
    (fnk
     [Atf-i Abf-i Zw-tf Zw-bf]
     (apply + (map (fn [ai zi] (* ai zi))
                   (concat Atf-i Abf-i)
                   (concat Zw-tf Zw-bf))))

    :Sw-str
    (fnk
     [Sw-st Ar-red Zw-r nr]
     (+ Sw-st
        (* Ar-red Zw-r)))

    :Sw-stc
    (fnk
     [Sw-st Asl-red Zw-sl Ar-red Zw-r]
     (+ Sw-st
        (* Ar-red Zw-r)
        (* Asl-red Zw-sl)))

    :Sw-stc-creep
    (fnk
     [Sw-st Asl-red-creep Zw-sl Ar-red Zw-r]
     (+ Sw-st
        (* Ar-red Zw-r)
        (* Asl-red-creep Zw-sl)))

    :Sw-stc-shr
    (fnk
     [Sw-st Asl-red-shr Zw-sl Ar-red Zw-r]
     (+ Sw-st
        (* Ar-red Zw-r)
        (* Asl-red-shr Zw-sl)))

    :Sw-stc-crack
    (fnk
     [Sw-st Ar-crack Zw-r]
     (+ Sw-st
        (* Ar-crack Zw-r)))

    ;; dZw

    :dZw-st
    (fnk
     [Sw-st Ast]
     (/ Sw-st Ast))

    :dZw-str
    (fnk
     [Sw-str Astr]
     (/ Sw-str Astr))

    :dZw-stc
    (fnk
     [Sw-stc Astc]
     (/ Sw-stc Astc))

    :dZw-stc-creep
    (fnk
     [Sw-stc-creep Astc-creep]
     (/ Sw-stc-creep Astc-creep))

    :dZw-stc-shr
    (fnk
     [Sw-stc-shr Astc-shr]
     (/ Sw-stc-shr Astc-shr))

    :dZw-stc-crack
    (fnk
     [Sw-stc-crack Astc-crack]
     (/ Sw-stc-crack Astc-crack))

    ;; Iw

    :Ic
    (fnk
     {{:keys [h b]} :slab}
     (I-rect {:h h :b b :dz (si/m 0)}))

    :Iw-st
    (fnk
     [top-flange web bottom-flange Zw-tf Zw-bf]
     (+ (I-rect {:b (:t web) :h (:h web)})
        (apply + (map (fn [{:keys [w t]} z]
                        (I-rect {:b w :h t :dz z}))
                      (concat top-flange bottom-flange)
                      (concat Zw-tf Zw-bf)))))

    :Iw-str
    (fnk
     [Iw-st Ar-red Zw-r]
     (+ Iw-st
        (* Ar-red Zw-r Zw-r)))

    :Iw-stc
    (fnk
     [Iw-st Zw-sl slab nr nc Ar Zw-r Ic]
     (let [{:keys [b h]} slab]
       (+ Iw-st
          (* Ar (- (/ nr) (/ nc)) Zw-r Zw-r)
          (/ (I-rect {:h h :b b :dz Zw-sl})
             nc))))

    :Iw-stc-creep
    (fnk
     [Iw-st Zw-sl slab nr nc-creep Ar Zw-r]
     (let [{:keys [b h]} slab]
       (+ Iw-st
          (* Ar (- (/ nr) (/ nc-creep)) Zw-r Zw-r)
          (/ (I-rect {:h h :b b :dz Zw-sl})
             nc-creep))))

    :Iw-stc-shr
    (fnk
     [Iw-st Zw-sl slab nr nc-shr Ar Zw-r]
     (let [{:keys [b h]} slab]
       (+ Iw-st
          (* Ar (- (/ nr) (/ nc-shr)) Zw-r Zw-r)
          (/ (I-rect {:h h :b b :dz Zw-sl})
             nc-shr))))

    :Iw-stc-crack
    (fnk
     [Iw-st Zw-r Ar-crack]
     (+ Iw-st
        (* Ar-crack Zw-r Zw-r)))

    ;; I

    :Ist
    (fnk
     [Iw-st dZw-st Ast]
     (- Iw-st (* Ast dZw-st dZw-st)))

    :Istr
    (fnk
     [Iw-str dZw-str Astr]
     (- Iw-str (* Astr dZw-str dZw-str)))

    :Istc
    (fnk
     [Iw-stc dZw-stc Astc]
     (- Iw-stc (* Astc dZw-stc dZw-stc)))

    :Istc-creep
    (fnk
     [Iw-stc-creep dZw-stc-creep Astc-creep]
     (- Iw-stc-creep (* Astc-creep dZw-stc-creep dZw-stc-creep)))

    :Istc-shr
    (fnk
     [Iw-stc-shr dZw-stc-shr Astc-shr]
     (- Iw-stc-shr (* Astc-shr dZw-stc-shr dZw-stc-shr)))

    :Istc-crack
    (fnk
     [Iw-stc-crack dZw-stc-crack Astc-crack]
     (- Iw-stc-crack (* Astc-crack dZw-stc-crack dZw-stc-crack)))

    ;; EI

    :EI-c
    (fnk
     [Ic Ec]
     (* Ec Ic))

    :EI-st
    (fnk
     [Ist Est]
     (* Ist Est))

    :EI-stc
    (fnk
     [Istc Est]
     (* Istc Est))

    :EI-stc-creep
    (fnk
     [Istc-creep Est]
     (* Istc-creep Est))

    :EI-stc-crack
    (fnk
     [Istc-crack Est]
     (* Istc-crack Est))

    ;; Z

    :Zs2-w
    (fnk
     [top-flange web]
     (- (apply + (/ (:h web) 2)
               (map :t top-flange))))

    :Zs1-w
    (fnk
     [bottom-flange web]
     (apply + (/ (:h web) 2)
            (map :t bottom-flange)))

    :Zs2-st
    (fnk
     [Zs2-w dZw-st]
     (- Zs2-w dZw-st))

    :Zs1-st
    (fnk
     [Zs1-w dZw-st]
     (- Zs1-w dZw-st))

    :Zs2-stc
    (fnk
     [Zs2-w dZw-stc]
     (- Zs2-w dZw-stc))

    :Zs1-stc
    (fnk
     [Zs1-w dZw-stc]
     (- Zs1-w dZw-stc))

    :Zs2-stc-crack
    (fnk
     [Zs2-w dZw-stc-crack]
     (- Zs2-w dZw-stc-crack))

    :Zs1-stc-crack
    (fnk
     [Zs1-w dZw-stc-crack]
     (- Zs1-w dZw-stc-crack))

    :Zs2-stc-creep
    (fnk
     [Zs2-w dZw-stc-creep]
     (- Zs2-w dZw-stc-creep))

    :Zs1-stc-creep
    (fnk
     [Zs1-w dZw-stc-creep]
     (- Zs1-w dZw-stc-creep))

    :Zc-st
    (fnk
     [Zw-sl dZw-st]
     (- Zw-sl dZw-st))

    :Zc-str
    (fnk
     [Zw-sl dZw-str]
     (- Zw-sl dZw-str))

    :Zc-stc
    (fnk
     [Zw-sl dZw-stc]
     (- Zw-sl dZw-stc))

    :Zr-stc
    (fnk
     [Zw-r dZw-stc]
     (- Zw-r dZw-stc))

    :Zc-stc-creep
    (fnk
     [Zw-sl dZw-stc-creep]
     (- Zw-sl dZw-stc-creep))

    :Zc-stc-shr
    (fnk
     [Zw-sl dZw-stc-shr]
     (- Zw-sl dZw-stc-shr))

    :Zr-stc-crack
    (fnk
     [Zw-r dZw-stc-crack]
     (- Zw-r dZw-stc-crack))

    :Zc-stc-crack
    (fnk
     [Zw-sl dZw-stc-crack]
     (- Zw-sl dZw-stc-crack))

    :Zr-st
    (fnk
     [Zw-r dZw-st]
     (- Zw-r dZw-st))

    ;; W

    ;; St
    :Wc-st
    (fnk
     [Zc-st Ist]
     (/ Ist Zc-st))

    :Ws1-st
    (fnk
     [Zs1-st Ist]
     (/ Ist Zs1-st))

    :Ws2-st
    (fnk
     [Zs2-st Ist]
     (/ Ist Zs2-st))

    ;; StC
    :Wc-stc
    (fnk
     [Zc-stc Istc]
     (/ Istc Zc-stc))

    :Wr-stc
    (fnk
     [Zr-stc Istc]
     (/ Istc Zr-stc))

    :Ws1-stc
    (fnk
     [Zs1-stc Istc]
     (/ Istc Zs1-stc))

    :Ws2-stc
    (fnk
     [Zs2-stc Istc]
     (/ Istc Zs2-stc))

    ;; Creep
    :Wc-stc-creep
    (fnk
     [Zc-stc-creep Istc-creep]
     (/ Istc-creep Zc-stc-creep))

    :Ws1-stc-creep
    (fnk
     [Zs1-stc-creep Istc-creep]
     (/ Istc-creep Zs1-stc-creep))

    :Ws2-stc-creep
    (fnk
     [Zs2-stc-creep Istc-creep]
     (/ Istc-creep Zs2-stc-creep))

    ;; Crack
    :Wr-stc-crack
    (fnk
     [Zr-stc-crack Istc-crack]
     (/ Istc-crack Zr-stc-crack))

    :Ws1-stc-crack
    (fnk
     [Zs1-stc-crack Istc-crack]
     (/ Istc-crack Zs1-stc-crack))

    :Ws2-stc-crack
    (fnk
     [Zs2-stc-crack Istc-crack]
     (/ Istc-crack Zs2-stc-crack))

    :Sst-w-top
    (fnk
     [Atf-i Zw-tf dZw-st]
     (apply + (map (fn [ai zi] (* ai (- zi dZw-st)))
                   Atf-i Zw-tf)))

    :Sst-w-bottom
    (fnk
     [Abf-i Zw-bf dZw-st]
     (apply + (map (fn [ai zi] (* ai (- zi dZw-st)))
                   Abf-i Zw-bf)))

    :Sst-w-max
    (fnk
     [web Sst-w-bottom dZw-st]
     (let [{:keys [h t]} web
           h (- (/ h 2) dZw-st)]
       (+ (* 1/2 t h h)
          Sst-w-bottom)))

    :Sst-w-min
    (fnk
     [Sst-w-top Sst-w-bottom]
     (min (abs Sst-w-top)
          (abs Sst-w-bottom)))

    :box?
    (fnk [] false)}))

  ;;
  ;;   ae, eta
  ;;


(def ^:private lazy-ae (lazy-compile steel/ae))

(def eta
  (fnk
   [Atf Abf Ast m Ry N]
   (let [kA (/ (min Atf Abf)
               (max Atf Abf))
         kN (/ (abs N)
               (* Ast m Ry))
         eta-above (table-2d
                    {:clip #{:l :r :t :b}
                     :xp [0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.7]
                     :yp [0 0.2 0.4 0.6 0.8 1.0]
                     :data [[1.00 1.00 1.00 1.00 1.00 1.00 0.99 0.98 0.96 0.95 0.92 0.88 0.83 0.87 0.63]
                            [1.00 1.00 1.00 1.02 1.03 1.04 1.05 1.06 1.07 1.06 1.05 1.02 0.99 0.90 0.75]
                            [1.00 1.04 1.08 1.12 1.14 1.16 1.19 1.20 1.21 1.20 1.18 1.16 1.13 1.09 1.04]
                            [1.00 1.10 1.19 1.28 1.35 1.40 1.44 1.46 1.47 1.46 1.45 1.42 1.39 1.35 1.30]
                            [1.00 1.20 1.39 1.55 1.70 1.83 1.93 1.98 2.00 2.02 2.01 1.99 1.97 1.91 1.84]
                            [1.00 1.29 1.63 2.04 2.47 2.86 3.20 3.38 3.49 3.56 3.57 3.53 3.43 3.29 3.05]]})
         eta-below (table-2d
                    {:clip #{:l :r :t :b}
                     :xp [0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.7]
                     :yp [0 0.2 0.4 0.6 0.8 1.0]
                     :data [[1.00 0.98 0.94 0.90 0.87 0.81 0.75 0.67 0.58 0.45 0.28 0.52 0.68 0.76 0.82]
                            [1.00 0.97 0.92 0.87 0.80 0.70 0.57 0.38 0.49 0.61 0.72 0.82 0.91 0.99 1.05]
                            [1.00 0.90 0.80 0.67 0.52 0.34 0.53 0.68 0.84 0.98 1.12 1.22 1.30 1.38 1.42]
                            [1.00 0.84 0.64 0.40 0.56 0.75 0.95 1.13 1.30 1.45 1.58 1.69 1.76 1.84 1.90]
                            [1.00 0.61 0.51 0.84 1.12 1.36 1.60 1.86 2.08 2.29 2.47 2.52 2.50 2.46 2.38]
                            [1.00 1.29 1.63 2.04 2.47 2.86 3.20 3.38 3.49 3.56 3.57 3.53 3.43 3.29 3.05]]})]
     (if (> Atf Abf)
       (eta-above kN kA)
       (eta-below kN kA)))))

;;
;;   Cases
;;

(def case-A
  (flow
    {:case (fnk [] :A)
     :valid-case
     (fnk
       [Rc mc sigma-c sigma-r mr Rr EI-st EI-c]
       (if (and (neg? sigma-c)
                (< (abs sigma-c) (* mc Rc))
                (neg? sigma-r)
                (< (abs sigma-r) (* mr Rr))
                (<= EI-c (* 0.2 EI-st)))
         true
         false))

     :sigma-c
     (fnk
       [M2 nc Wc-stc sigma-ci]
       {:post [(dim/stress? %)]}
       (+ (/ M2 nc Wc-stc)
          sigma-ci))

     :sigma-r
     (fnk
       [M2 Wr-stc sigma-ri nr]
       {:post [(dim/stress? %)]}
       (+ (/ M2 Wr-stc nr)
          sigma-ri))

     :M
     (fnk
       [M1 M2]
       (+ M1 M2))

     :sigma-s2
     (fnk [M Zc-st Ncr ae4 Ws2-st Ast]
          {:post [(dim/stress? %)]}
          (+ (/ (+ M (* Zc-st Ncr))
                (* ae4 Ws2-st))
             (/ Ncr Ast)))

     :sigma-s1
     (fnk [M Zc-st Ncr ae3 Ws1-st Ast]
          {:post [(dim/stress? %)]}
          (+ (/ (+ M (* Zc-st Ncr))
                (* ae3 Ws1-st))
             (/ Ncr Ast)))

     :Ncr
     (fnk [Asl sigma-c Ar sigma-r]
          (- (+ (* (- Asl Ar) sigma-c)
                (* Ar sigma-r))))

     :ae
     (fnk
       {:keys [gamma-m Abf Ast Sst-w-max Ryn Aw Ist m box? Q Sst-w-min Atf web] :as args}
       (:ae (lazy-ae {
                      :Abf Abf :Atf Atf :Ast Ast :Aw Aw
                      :Smax Sst-w-max :Smin Sst-w-min
                      :Ryn Ryn :m m :gamma-m gamma-m
                      :Ist Ist
                      :box? box?
                      :Q Q
                      :tw (:t web) :hw (:h web)})))

     :ae4
     (fnk
       [ae3 m1]
       (clip {:min 1.0} (/ ae3 m1)))

     :ae3
     (fnk
       [ae eta]
       (+ 1 (* eta (- ae 1))))

     :m1
     (fnk
       [m Ry mc Rc sigma-c Asl Atf]
       (clip {:max 1.2}
             (+ 1 (* (/ (- (* mc Rc) sigma-c)
                        (* m Ry))
                     (/ Asl Atf)))))

     :eta
     (fnk
       {:keys [Atf Abf Ast m Ry Ncr] :as input}
       (eta (merge input {:N Ncr})))

     :m1mRy
     (fnk
       [m1 m Ry]
       (* m1 m Ry))

     :mRy
     (fnk
       [m Ry]
       (* m Ry))}))

(def case-E
  (flow
    {:case (fnk [] :E)
     :valid-case
     (fnk
       [Rc mc sigma-c]
       (if (and (pos? sigma-c)
                (>= sigma-c (* 0.1 mc Rc)))
         true
         false))

     :sigma-c
     (fnk
       [M2 nc Wc-stc sigma-ci]
       {:post [(dim/stress? %)]}
       (+ (/ M2 nc Wc-stc)
          sigma-ci))

     :sigma-r
     (fnk
       [M2 sigma-ri sigma-ci
        Astc-crack Wr-stc-crack
        Zc-stc-crack Zr-stc-crack
        Asl psi-crack nr ]
       {:post [(dim/stress? %)]}
       (+ (/ (+ M2 (* Zc-stc-crack Asl sigma-ci))
             (* psi-crack nr Wr-stc-crack))
          (/ (* Asl sigma-ci)
             (* psi-crack nr Astc-crack))
          sigma-ri))

     :M
     (fnk
       [M1 M2]
       (+ M1 M2))

     :sigma-s2
     (fnk [M Zr-st NrR ae3-s2 Ws2-st Ast]
          {:post [(dim/stress? %)]}
          (- (/ (- M (* Zr-st NrR))
                (* ae3-s2 Ws2-st))
             (/ NrR Ast)))

     :sigma-s1
     (fnk [M Zr-st Nr ae3-s1 Ws1-st Ast]
          {:post [(dim/stress? %)]}
          (- (/ (- M (* Zr-st Nr))
                (* ae3-s1 Ws1-st))
             (/ Nr Ast)))

     :NrR
     (fnk [Ar Rr]
          (* Ar Rr))

     :Nr
     (fnk [Ar sigma-r Rr]
          (min (* Ar sigma-r)
               (* Ar Rr)))

     :ae
     (fnk
       {:keys [gamma-m Abf Ast Sst-w-max Ryn Aw Ist m box? Q Sst-w-min Atf web] :as args}
       (:ae (lazy-ae {
                      :Abf Abf :Atf Atf :Ast Ast :Aw Aw
                      :Smax Sst-w-max :Smin Sst-w-min
                      :Ryn Ryn :m m :gamma-m gamma-m
                      :Ist Ist
                      :box? box?
                      :Q Q
                      :tw (:t web) :hw (:h web)})))

     :ae3-s2
     (fnk
       [ae eta-s2]
       (+ 1 (* eta-s2 (- ae 1))))

     :ae3-s1
     (fnk
       [ae eta-s1]
       (+ 1 (* eta-s1 (- ae 1))))

     :eta-s2
     (fnk
       {:keys [Atf Abf Ast m Ry NrR] :as input}
       (eta (merge input {:N NrR})))

     :eta-s1
     (fnk
       {:keys [Atf Abf Ast m Ry Nr] :as input}
       (eta (merge input {:N Nr})))}))

;;
;; Creep, shrinkage, temperature
;;

(def creep
  (flow
    {:sigma-c-cr
     (fnk
       [alpha sigma-c]
       (- (* alpha sigma-c)))

     :alpha
     (fnk
       [phi-cr v]
       (/ phi-cr
          (+ (* 0.5 phi-cr) v 1)))

     :v
     (fnk
       [Asl-red Astr Zc-str Istr]
       (* Asl-red
          (+ (/ Astr)
             (/ (pow Zc-str 2)
                Istr))))

     :phi-cr
     (fnk
       [gamma-f cn-lim Ec]
       (* gamma-f Ec cn-lim))

     :gamma-f (fnk [] 1.1)

     :cn-lim
     (fnk
       [cn xi1 xi2 xi3 xi4]
       (* cn xi1 xi2 xi3 xi4))

     :xi1
     (fnk
       [Rcl:Rc]
       ((table-1d
          {:xp [0.5 0.6 0.7 0.8 0.9 1.0]
           :clip #{:l :r}
           :data [1.7 1.6 1.4 1.25 1.15 1.0]})
         Rcl:Rc))

     :xi2
     (fnk
       [load-age]
       ((table-1d
          {:xp [3 7 28 60 90 180 360]
           :clip #{:l :r}
           :data [1 1 1 0.8 0.7 0.6 0.5]})
         load-age))

     :xi3
     (fnk
       [Asl slab]
       ((table-1d
          {:xp [2.5 5.0 7.5 10.0 12.5 15.0 20.0]
           :clip #{:l :r}
           :data [1 0.85 0.76 0.72 0.69 0.67 0.64]})
         (:magnitude (si/cm (/ Asl (* 2 (+ (:h slab) (:b slab)))))) ))

     :xi4
     (fnk
       [H]
       ((table-1d
              {:xp [40 50 60 70 80 90 100]
               :clip #{:l :r}
               :data [1.33 1.25 1.15 1.0 0.85 0.7 0.51]})
         H))}))

(def shrinkage
  (flow
    {:eps-shr
     (fnk [] 2e-4)

     :sigma-c-shr
     (fnk
       [eps-shr Ec-shr Astr Astc-shr Sshr Istc-shr Zc-stc-shr]
       (* eps-shr Ec-shr
          (+ (/ Astr Astc-shr)
             (* (/ Sshr Istc-shr) Zc-stc-shr))))

     :Sshr
     (fnk
       [Astr Zstr-stc-shr]
       (* Astr Zstr-stc-shr))

     :Zstr-stc-shr
     (fnk
       [dZw-str dZw-stc-shr]
       (- dZw-str dZw-stc-shr))}))
