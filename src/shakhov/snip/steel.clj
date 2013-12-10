(ns shakhov.snip.steel
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

;;;;;;; ae

(def ae
  (flow
    {:ae
     (fnk
       [tau-m m Rs ae1 alpha a b]
       (clip {:min 0.0 :max ae1}
             (if (<= tau-m (* 0.25 Rs))
               ae1
               (* ae1 (/ (+ (sqrt (- 1 (* alpha alpha)))
                            (* 2 a b))
                         (+ 1 (* 2 a)))))))

     :ae1
     (fnk
       [Atf Abf Aw Ast]
       (let [Af-min (min Atf Abf)
             x (/ (+ Af-min Aw) Ast)
             y (/ Af-min Aw)
             table-61 (table-2d
                        {:xp [0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]
                         :yp [0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0]
                         :clip #{:l :r :t :b}
                         :data [[1.243 1.248 1.253 1.258 1.264 1.269 1.274 1.279 1.283 1.267 1.243]
                                [1.187 1.191 1.195 1.199 1.202 1.206 1.209 1.212 1.214 1.160 1.160]
                                [1.152 1.155 1.158 1.162 1.165 1.168 1.170 1.172 1.150 1.150 1.150]
                                [1.128 1.131 1.133 1.136 1.139 1.142 1.144 1.145 1.097 1.097 1.097]
                                [1.110 1.113 1.115 1.118 1.120 1.123 1.125 1.126 1.069 1.069 1.069]
                                [1.097 1.099 1.102 1.104 1.106 1.109 1.110 1.106 1.061 1.061 1.061]
                                [1.087 1.089 1.091 1.093 1.095 1.097 1.099 1.079 1.079 1.079 1.079]
                                [1.078 1.080 1.082 1.084 1.086 1.088 1.090 1.055 1.055 1.055 1.055]
                                [1.071 1.073 1.075 1.077 1.079 1.081 1.082 1.044 1.044 1.044 1.044]
                                [1.065 1.067 1.069 1.071 1.073 1.074 1.076 1.036 1.036 1.036 1.036]
                                [1.060 1.062 1.064 1.066 1.067 1.069 1.071 1.031 1.031 1.031 1.031]
                                [1.035 1.036 1.037 1.038 1.039 1.040 1.019 1.019 1.019 1.019 1.019]
                                [1.024 1.025 1.026 1.027 1.028 1.029 1.017 1.017 1.017 1.017 1.017]
                                [1.019 1.019 1.020 1.021 1.021 1.022 1.015 1.015 1.015 1.015 1.015]
                                [1.015 1.015 1.016 1.017 1.018 1.018 1.018 1.018 1.018 1.018 1.018]]})]
         (table-61 x y)))

     :ae2
     (fnk
       [Smax Smin]
       (- 1.25 (* 0.25 (/ Smin Smax))))

     :Rs
     (fnk
       [Ryn gamma-m]
       (/ (* 0.58 Ryn)
          gamma-m))

     :tau-m
     (fnk
       [Q hw tw]
       (/ Q hw tw))

     :alpha
     (fnk
       [Q Qu]
       (/ Q Qu))

     :a
     (fnk
       [Atf Abf Aw]
       (/ (+ Atf Abf) Aw))

     :b
     (fnk
       [alpha box?]
       (Math/sqrt (- 1 (* (if box? 0.0625 0.25) alpha alpha))))

     :Qu
     (fnk
       [Rs m ae2 Ist tw Smax]
       (/ (* Rs m ae2 Ist tw)
          Smax))}))
