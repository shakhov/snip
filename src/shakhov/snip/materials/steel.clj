(ns shakhov.snip.materials.steel
  (:require [shakhov.snip.units :as si]))

(def st-15HSND
  {:steel :st-15HSND
   :Est (si/MPa 2.06e5)
   :Ry  (si/MPa 285.0)
   :Ru  (si/MPa 400.0)
   :Ryn (si/MPa 330.0)
   :Run (si/MPa 470.0)
   :gamma-m 1.165})

(def st-10HSND
  {:steel :st-10HSND
   :Est (si/MPa 2.06e5)
   :Ry  (si/MPa 350.0)
   :Ru  (si/MPa 450.0)
   :Ryn (si/MPa 390.0)
   :Run (si/MPa 510.0)
   :gamma-m 1.125})