(ns shakhov.snip.materials.rebar
  (:refer-clojure :exclude [*])
  (:require [shakhov.snip.units :as si])
  (:use [clojure.algo.generic.arithmetic :only [*]]))

(defn- A
  [d]
  (* 1/4 Math/PI d d))

(def AIII
  {:Er (si/MPa 1.96e5)
   :Rr (si/MPa 350.0)
   :A A})
