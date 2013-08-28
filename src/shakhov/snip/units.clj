(ns shakhov.snip.units
  (:refer-clojure :exclude [time force * /])
  (:use shakhov.units.core)
  (:use shakhov.snip.dimensions)
  (:use [clojure.algo.generic.arithmetic :only [* /]]
            [clojure.algo.generic.math-functions :only [pow sqrt]]))

(def-unit-system 
  civil-units
  length m
  mass kg
  time s
  temperature degC)

(def-unit
  civil-units
  rad (as-unit civil-units angle)
  N (as-unit civil-units force)
  Pa (/ N m m)
  J (* N m))

(def-prefixed-units
  civil-units [m N Pa J]
  M 1e6
  k 1e3
  m 1e-3)

(def-unit civil-units 
  cm (* 0.01 m))

(let [g (* 9.80665 (/ m (pow s 2)))]
  (def-unit
    civil-units
    kgf  (* kg g)
    tonf (* 1000.0 kgf)))

(def-unit 
  civil-units
  kgf:cm3 (/ kgf cm cm cm)
  tonf:m3 (/ tonf m m m)
  kgf:cm2 (/ kgf cm cm)
  tonf:m2 (/ tonf m m)
  kgf:cm (/ kgf cm)
  tonf:m (/ tonf m)
  kN*m (* kN m)
  MN*m (* MN m)
  tonf*m (* tonf m))

(def-unit 
  civil-units
  deg (* (/ Math/PI 180.0) rad))
