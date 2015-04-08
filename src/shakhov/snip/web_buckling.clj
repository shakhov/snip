(ns shakhov.snip.web-buckling
  (:refer-clojure :exclude [time force + - * / < > <= >= = zero? pos? neg? sgn abs
                            sin cos tan asin acos atan exp log min max])

  (:use [shakhov.flow.core]
        [shakhov.snip.utils])

  (:use [clojure.algo.generic.arithmetic :only [+ - * /]]
        [clojure.algo.generic.comparison :only [< > <= >= = zero? pos? neg? min max]]
        [clojure.algo.generic.math-functions :only [pow sqrt sgn abs sin cos tan
                                                    asin acos atan exp log]])
  (:require [shakhov.snip.dimensions :as dim]
            [shakhov.snip.units :as si])

  (:use [shakhov.snip.pprint]))

(def table-2
  (table-1d
    {:xp [0 0.5 1.0 1.5 2.0 3.0 4.0]
     :data [1.00 1.05 1.10 1.15 1.20 1.30 1.40]}))

(def table-4
  (table-1d
    {:xp [0.25 0.5 1.0 2.0 4.0 10.0 10.001]
     :data [1.21 1.33 1.46 1.55 1.60 1.63 1.65]
     :clip #{:r}}))

(def table-5
  (table-2d
    {:xp [0.4 0.5 0.6 0.67 0.75 0.8 0.9 1.0 1.2 2.0]
     :yp [0 0.67 0.80 1.0 1.33 2.0 3.0 4.0]
     :data [[8.41 6.25 5.14 4.75 4.36 4.2 4.04 4.0 4.34 4.0]
            [10.8 8.0 7.1 6.6 6.1 6.0 5.9 5.8 6.1 5.8]
            [13.3 9.6 8.3 7.7 7.1 6.9 6.7 6.6 7.1 6.6]
            [15.1 11.0 9.7 9.0 8.4 8.1 7.9 7.8 8.4 7.8]
            [18.7 14.2 12.9 12.0 11.0 11.2 11.1 11.0 11.5 11.0]
            [29.1 25.6 24.1 23.9 24.1 24.4 25.6 25.6 24.1 23.9]
            [54.3 54.5 58.0 53.8 53.8 53.8 53.8 53.8 53.8 53.8]
            [95.7 95.7 95.7 95.7 95.7 95.7 95.7 95.7 95.7 95.7]]}))

(def table-6
  (table-2d
    {:xp [0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.18 0.20 0.25 0.30 0.35]
     :yp [0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.5 2.0]
     :data [[1.70 1.67 1.65 1.63 1.61 1.60 1.60 1.60 1.60 1.60 1.60 1.60]
            [1.98 1.93 1.89 1.85 1.82 1.80 1.79 1.78 1.76 1.72 1.71 1.69]
            [2.23 2.17 2.11 2.06 2.02 1.98 1.96 1.93 1.89 1.82 1.79 1.76]
            [2.43 2.35 2.28 2.22 2.17 2.12 2.10 2.05 2.01 1.91 1.86 1.82]
            [2.61 2.51 2.43 2.36 2.30 2.24 2.21 2.16 2.11 1.98 1.92 1.87]
            [2.74 2.64 2.55 2.47 2.40 2.34 2.31 2.24 2.17 2.04 1.97 1.91]
            [2.79 2.68 2.59 2.51 2.43 2.37 2.33 2.26 2.19 2.05 1.98 1.91]
            [2.84 2.73 2.63 2.54 2.46 2.39 2.35 2.28 2.21 2.05 1.98 1.91]
            [2.86 2.75 2.65 2.56 2.48 2.41 2.37 2.30 2.25 2.07 1.99 1.91]
            [2.86 2.75 2.65 2.55 2.47 2.40 2.36 2.28 2.20 2.05 1.96 1.88]]}))

(def table-7
  (table-2d
    {:xp [0.4 0.6 0.8 1.0 1.5 2.0]
     :yp [0.25 0.5 1.0 4.0 10.0]
     :data [[1.19 1.19 1.20 1.20 1.19 1.18]
            [1.24 1.29 1.30 1.32 1.32 1.32]
            [1.28 1.36 1.41 1.47 1.52 1.56]
            [1.32 1.45 1.57 1.73 1.97 2.21]
            [1.34 1.49 1.65 1.88 2.51 2.95]]}))

(def table-8
  (table-1d
    {:xp [0.4 0.5 0.6 0.7 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.5]
     :data [4.88 5.12 5.37 5.59 5.80 6.26 6.87 7.69 8.69 9.86 11.21 15.28]
     :clip #{:r}}))

(def table-9
  (table-2d
    {:xp [0.5 0.67 1.0 2.0 2.5]
     :yp [0.25 0.5 1.0 2.0 5.0 10.0 10.001]
     :data [[1.014 1.063 1.166 1.170 1.192]
            [1.016 1.075 1.214 1.260 1.300]
            [1.017 1.081 1.252 1.358 1.416]
            [1.018 1.085 1.275 1.481 1.516]
            [1.018 1.088 1.292 1.496 1.602]
            [1.018 1.088 1.298 1.524 1.636]
            [1.018 1.089 1.303 1.552 1.680]]}))

(def table-10
  (table-1d
    {:xp [0.5 1.0 2.0 5.0 10.0]
     :data [1.16 1.22 1.27 1.31 1.35]
     :clip #{:r}}))

(def table-11
  (table-1d
    {:xp [0.5 0.8 1.0 1.5 2.0]
     :data [1.07 1.18 1.31 1.52 1.62]
     :clip #{:r}}))

(def table-12
  (table-2d
    {:xp [0.5 0.6 0.9 1.0 1.5 2.0 2.5 3.0]
     :yp [2.0 4.0]
     :data [[1.06 1.07 1.13 1.17 1.31 1.32 1.29 1.25]
            [1.06 1.07 1.14 1.19 1.38 1.44 1.43 1.39]]}))

(def table-13o
  (table-1d
    {:xp [0.4 0.5 0.6 0. 0.8 1.0 1.5 2.0]
     :data [1240 1380 1520 1650 1820 2240 3860 6300]}))

(def table-13i
  (table-1d
    {:xp [0.4 0.5 0.6 0. 0.8 1.0 1.5 2.0]
     :data [920 970 1020 1060 1100 1190 1530 2130]}))

;;
;;   Sigma-cr
;;

(def sigma-cr
  (fnk
    [sigma-cr-ef steel Est m]
    (let [[s1 s2 A B C D E] (case steel
                              :st-10HSND [229 591 -215.8 1.238 -0.1091e-3 0.03677 1.561e-3]
                              :st-15HSND [207 524 -201.2 1.024  0.0795e-3 0.03572 1.290e-3]
                              :else      [196 385 -170.7 0.6375 0.4048e-3 0.03114 0.9419e-3])]
      (cond
        (<= (si/MPa 0) sigma-cr-ef (si/MPa s1))
        (* 0.9 sigma-cr-ef m)

        (<= (si/MPa s1) sigma-cr-ef (si/MPa s2))
        (* (+ (* A (pow (/ sigma-cr-ef Est) 2))
              (* B (/ sigma-cr-ef Est))
              C)
           Est m)

        (<= (si/MPa s2) sigma-cr-ef)
        (* (+ (* D (/ sigma-cr-ef Est))
              E)
           Est m)))))

(def sigma-x-cr
  (fnk
    [sigma-x-cr-ef steel Est m]
    (sigma-cr {:sigma-cr-ef sigma-x-cr-ef :steel steel :Est Est :m m})))

(def sigma-y-cr
  (fnk
    [sigma-y-cr-ef steel Est m]
    (sigma-cr {:sigma-cr-ef sigma-y-cr-ef :steel steel :Est Est :m m})))

(def tau-xy-cr
  (fnk
    [tau-xy-cr-ef steel Est m]
    (* 0.6 (sigma-cr {:sigma-cr-ef ( / tau-xy-cr-ef 0.6) :steel steel :Est Est :m m}))))

;;
;;   Buckling
;;

(def omega-1
  (fnk
    [xi]
    (table-2 xi)))

(def xi
  (fnk
    [sigma-x-min sigma-x-max]
    (- 1 (/ sigma-x-min sigma-x-max))))

(def gamma
  (fnk
    [beta bf tf h-ef tw]
    (* beta (/ bf h-ef) (pow (/ tf tw) 3))))

(def mu
  (fnk
    [a h-ef]
    (/ a h-ef)))

(def sigma-x
  (fnk
    [sigma-x-max]
    sigma-x-max))

(def eps-x
  (fnk
    [mu xi]
    (table-5 mu xi)))

(def z-y
  (fnk

   [mu]
    (table-8 mu)))

(def chi-tau
  (fnk
    [mu gamma connection]
    (case connection
      :welded (table-9 mu gamma)
      :bolted 1.0)))

(def i-y
  (fnk
    [a h-ef mu]
    (cond
      (< 0.4 mu 0.7) 2.0
      (<= 0.7 mu)    1.0)))

(def sigma-x-cr-ef-11
  (fnk
    [chi-x eps-x Est tw h-ef]
    (* 9.05e-5 chi-x eps-x (pow (/ (* 100 tw) h-ef) 2) Est)))

(def sigma-y-cr-ef-12
  (fnk
    [dzeta-y z-y chi-y Est tw a]
    (* 9.05e-5 dzeta-y z-y chi-y (pow (/ (* 100 tw) a) 2) Est)))

(def tau-xy-cr-ef-13
  (fnk
    [a h-ef mu Est chi-tau tw]
    (let [d (min a h-ef)
          mu1 (if (< a h-ef) (/ mu) mu)]
      (* 0.476e-6 chi-tau (+ 1020 (/ 760 (pow mu1 2))) (pow (/ (* 100 tw) d) 2) Est))))

(def sigma-y-cr-ef-15
  (fnk
    [tw a mu chi-y i-y Est l-ef h-ef]
    (let [k  (if (not (zero? l-ef)) 1.55 1.00)
          mu (if (and (< 1.0 k)
                      (> a (* 2 (+ h-ef l-ef))))
               (/ (* 2 (+ h-ef l-ef)) h-ef)
               mu)]
      (* k 9.05e-5
         chi-y
         (pow (/ (+ 1 (pow (* mu i-y) 2))
                 (* mu i-y)) 2)
         (pow (/ (* 100 tw) a) 2)
         Est))))

(def sigma-y-cr-ef-17
  (fnk
    [delta tw a Est]
    (* 0.476e-6 delta (pow (/ (* 100 tw) a) 2) Est)))

(def buckling-10
  (fnk
    [sigma-x sigma-y tau-xy sigma-x-cr sigma-y-cr tau-xy-cr omega-1 omega-2]
    (sqrt (+ (pow (+ (/ sigma-x
                        omega-1 sigma-x-cr)
                     (/ sigma-y
                        sigma-y-cr))
                  2)
             (pow (/ (* 0.9 tau-xy)
                     omega-2 tau-xy-cr)
                  2)))))

(def buckling-14
  (fnk
    [sigma-x sigma-y tau-xy sigma-x-cr sigma-y-cr tau-xy-cr omega-1]
    (+ (/ sigma-x
          omega-1 sigma-x-cr)
       (/ sigma-y
          sigma-y-cr)
       (pow (/ (* 0.9 tau-xy)
               tau-xy-cr)
            2))))

(def buckling-16
  (fnk
    [sigma-y tau-xy sigma-y-cr tau-xy-cr]
    (sqrt (+ (/ sigma-y sigma-y-cr)
             (pow (/ (* 0.9 tau-xy)
                     tau-xy-cr)
                  2)))))

;;  Cases


(def single-cell
  (flow
    {:buckling buckling-10
     :omega-1 omega-1
     :omega-2
     (fnk
       [hw tw]
       (if (> (/ hw tw) 100)
         (+ 1.0 (* 0.5 (- (/ hw tw) 0.5)))
         1.0))

     :sigma-x sigma-x
     :xi xi
     :mu mu
     :gamma gamma
     :h-ef (fnk [hw] hw)

     :sigma-x-cr-ef sigma-x-cr-ef-11
     :eps-x eps-x
     :chi-x
     (fnk
       [gamma connection]
       (case connection
         :bolted 1.4
         (table-4 gamma)))

     :sigma-y-cr-ef sigma-y-cr-ef-12
     :z-y z-y
     :dzeta-y
     (fnk
       [l-ef h-ef mu]
       (if (not (zero? l-ef))
         (table-6 (* 1.04 (/ l-ef h-ef)) mu)
         1.0))

     :chi-y
     (fnk
       [mu gamma connection]
       (case connection
         :bolted 1.0
         (table-7 mu gamma)))

     :tau-xy-cr-ef tau-xy-cr-ef-13
     :chi-tau chi-tau

     :sigma-x-cr sigma-x-cr
     :sigma-y-cr sigma-y-cr
     :tau-xy-cr tau-xy-cr}))

;; Two cells

(def two-cells-compression
  (flow
    {:buckling buckling-14
     :omega-1 omega-1
     :sigma-x sigma-x
     :xi xi
     :mu mu
     :gamma gamma

     :sigma-x-cr-ef sigma-x-cr-ef-11
     :eps-x eps-x
     :chi-x
     (fnk
       [gamma connection]
       (case connection
         :bolted 1.3
         :concrete-slab 1.35
         (table-10 gamma)))

     :sigma-y-cr-ef sigma-y-cr-ef-15
     :i-y i-y
     :chi-y
     (fnk
       [gamma mu connection]
       (case connection
         :bolted (table-11 mu)
         :concrete-slab (table-11 mu)
         :welded (table-12 mu gamma)))

     :tau-xy-cr-ef tau-xy-cr-ef-13
    :chi-tau
     (fnk
       [mu gamma connection]
       (/ (+ 1 (chi-tau {:mu mu :gamma gamma :connection connection}))
          2))

     :sigma-x-cr sigma-x-cr
     :sigma-y-cr sigma-y-cr
     :tau-xy-cr tau-xy-cr}))

(def two-cells-tension
  (flow
    {:buckling buckling-10
     :omega-1 omega-1
     :omega-2 (fnk [] 1.0)
     :sigma-x sigma-x
     :xi xi
     :mu mu

     :sigma-x-cr-ef sigma-x-cr-ef-11
     :eps-x eps-x
     :chi-x (fnk [] 1.0)

     :sigma-y-cr-ef sigma-y-cr-ef-12
     :chi-y (fnk [] 1.0)
     :dzeta-y (fnk [mu] (table-6 0.35 mu))
     :z-y z-y

     :tau-xy-cr-ef tau-xy-cr-ef-13
     :chi-tau (fnk [] 1.0)

     :sigma-x-cr sigma-x-cr
     :sigma-y-cr sigma-y-cr
     :tau-xy-cr tau-xy-cr}))

;;; 3

(def three-cells-outer-compression
  (merge two-cells-compression
         {:chi-x
          (fnk
            [gamma connection]
            (case connection
              :bolted 1.4
              (table-4 gamma)))
          :chi-tau chi-tau}))

(def three-cells-inner-compression
  (merge three-cells-outer-compression
         {:chi-x (fnk [] 1.0)
          :chi-y (fnk [] 1.0)
          :chi-tau (fnk [] 1.0)
          :gamma (fnk [] 0.0)}))

(def three-cells-inner-tension-compression
  (merge two-cells-tension
         {:omega-1 (fnk [] 1.0)}))

;; Tension

(def three-cells-outer-tension
  (flow
    {:buckling buckling-16
     :mu mu

     :sigma-y-cr-ef sigma-y-cr-ef-17
     :delta (fnk [mu] (table-13o mu))

     :tau-xy-cr-ef
     (fnk
       [a h-ef mu Est tw]
       (let [d (min a h-ef)
             mu1 (if (< a h-ef) (/ mu) mu)]
         (* 0.476e-6 (+ 1250 (/ 950 (pow mu1 2))) (pow (/ (* 100 tw) d) 2) Est)))

     :sigma-y-cr sigma-y-cr
     :tau-xy-cr tau-xy-cr}))

(def three-cells-inner-tension
  (flow
    {:buckling buckling-16
     :mu mu

     :sigma-y-cr-ef sigma-y-cr-ef-17
     :delta (fnk [mu] (table-13i mu))

     :tau-xy-cr-ef
     (fnk
       [a h-ef mu Est tw]
       (let [d (min a h-ef)
             mu1 (if (< a h-ef) (/ mu) mu)]
         (* 0.476e-6 (+ 1020 (/ 760 (pow mu1 2))) (pow (/ (* 100 tw) d) 2) Est)))

     :sigma-y-cr sigma-y-cr
     :tau-xy-cr tau-xy-cr}))


(def cell-buckling
  {:one-cell    single-cell
   :two-cells   {:tension     two-cells-tension
                 :compression two-cells-compression}
   :three-cells {:outer {:tension     three-cells-outer-tension
                         :compression three-cells-outer-compression}
                 :inner {:tension     three-cells-inner-tension
                         :compression three-cells-inner-compression
                         :tension-compression three-cells-inner-tension-compression}}})

(def cell-stress
  (flow
    {:z-top
     (fnk
       [z1 dZw]
       (- z1 dZw))

     :z-bottom
     (fnk
       [z2 dZw]
       (- z2 dZw))

     :stress-state
     (fnk [sigma-x-1 sigma-x-2]
          (cond
            (and (pos? sigma-x-1) (pos? sigma-x-2)) :compression
            (and (neg? sigma-x-1) (neg? sigma-x-2)) :tension
            (neg? (* sigma-x-1 sigma-x-2)) :tension-compression))

     :sigma-x-1
     (fnk
       [z-top Ix M N A dZw]
       (- (+ (/ N A)
             (* (/ M Ix) z-top))))

     :sigma-x-2
     (fnk
       [z-bottom Ix M N A dZw]
       (- (+ (/ N A)
             (* (/ M Ix) z-bottom))))

     :sigma-x-max
     (fnk
       [sigma-x-1 sigma-x-2]
       (max sigma-x-1 sigma-x-2))

     :sigma-x-min
     (fnk
       [sigma-x-1 sigma-x-2]
       (min sigma-x-1 sigma-x-2))

     :tau-xy
     (fnk
       [tau-1 tau-2]
       (/ (+ tau-1 tau-2) 2))

     :tau-1
     (fnk
       [S1 tw Ix Q]
       (/ (* Q (abs S1))
          (* tw Ix)))

     :tau-2
     (fnk
       [S2 tw Ix Q]
       (/ (* Q (abs S2))
          (* tw Ix)))

     :tau-max
     (fnk
       [Sw-max tw Ix Q]
       (/ (* Q (abs Sw-max))
          (* tw Ix)))

     :S1
     (fnk
       [z-top Sw-bottom hw tw dZw]
       (let [Zw-bottom (- (/ hw 2) dZw)
             hw-z (- Zw-bottom z-top)]
         (+ Sw-bottom
            (* tw hw-z
               (/ (+ Zw-bottom z-top) 2)))))

     :S2
     (fnk
       [z-bottom Sw-bottom hw tw dZw]
       (let [Zw-bottom (- (/ hw 2) dZw)
             hw-z (- Zw-bottom z-bottom)]
         (+ Sw-bottom
            (* tw hw-z
               (/ (+ Zw-bottom z-bottom) 2)))))

     :sigma-y
     (fnk
       [sigma-y0 hw dZw l-ef z-top]
       (let [Zw-top (- (+ (/ hw 2) dZw))
             h0 (- z-top Zw-top)]
         (if (zero? h0)
           sigma-y0
           (if (zero? l-ef)
             (let [v (/ h0 hw)]
               (* sigma-y0
                  (+ 1
                     (* -3 (pow v 2))
                     (* +2 (pow v 3)))))
             (let [v (/ h0 hw)
                   alpha (* 0.5 (/ l-ef hw))]
               (* sigma-y0
                  (/ 2 Math/PI)
                  (- (atan (/ alpha v))
                     (* 3 (pow v 2)
                        (- 1 (* 2/3 v))
                        (atan alpha)))))))))}))
