(ns shakhov.snip.materials.concrete
  (:refer-clojure :exclude [time force + - * /])
  (:require [shakhov.snip.units :as si])
  (:use [clojure.algo.generic.arithmetic :only [+ - * /]]))
  
  

(def B20 {:Rc      (si/MPa 10.5)  
          :Rct     (si/MPa 0.85)
          :Rc-ser  (si/MPa 15.0)
          :Rct-ser (si/MPa 1.40)
          :Rc-sh   (si/MPa 1.95)
          :Rc-mc1  nil
          :Rc-mc2  (si/MPa 8.8)
          :Ec      (si/MPa 27.0e3)
          :n'      22.5
          :cn      ((/ si/MPa) 115e-6)})
         
(def B22.5 {:Rc      (si/MPa 11.75)  
            :Rct     (si/MPa 0.90)
            :Rc-ser  (si/MPa 16.8)
            :Rct-ser (si/MPa 1.50)
            :Rc-sh   (si/MPa 2.30)
            :Rc-mc1  nil
            :Rc-mc2  (si/MPa 10.3)
            :Ec      (si/MPa 28.5e3)
            :n'      20.0
            :cn      ((/ si/MPa) 107e-6)})
         
(def B25 {:Rc      (si/MPa 13.0)  
          :Rct     (si/MPa 0.95)
          :Rc-ser  (si/MPa 18.5)
          :Rct-ser (si/MPa 1.60)
          :Rc-sh   (si/MPa 2.50)
          :Rc-mc1  (si/MPa 13.7)
          :Rc-mc2  (si/MPa 11.8)
          :Ec      (si/MPa 30.0e3)
          :n'      20.0
          :cn      ((/ si/MPa) 100e-6)})
         
(def B27.5 {:Rc      (si/MPa 14.3)  
            :Rct     (si/MPa 1.05)
            :Rc-ser  (si/MPa 20.5)
            :Rct-ser (si/MPa 1.70)
            :Rc-sh   (si/MPa 2.75)
            :Rc-mc1  (si/MPa 15.2)
            :Rc-mc2  (si/MPa 13.2)
            :Ec      (si/MPa 31.5e3)
            :n'      17.0
            :cn      ((/ si/MPa) 92e-6)})
         
(def B30 {:Rc      (si/MPa 15.5)  
          :Rct     (si/MPa 1.10)
          :Rc-ser  (si/MPa 22.0)
          :Rct-ser (si/MPa 1.80)
          :Rc-sh   (si/MPa 2.90)
          :Rc-mc1  (si/MPa 16.7)
          :Rc-mc2  (si/MPa 14.6)
          :Ec      (si/MPa 32.5e3)
          :n'      15.0
          :cn      ((/ si/MPa) 84e-6)})
         
(def B35 {:Rc      (si/MPa 17.5)  
          :Rct     (si/MPa 1.15)
          :Rc-ser  (si/MPa 25.5)
          :Rct-ser (si/MPa 1.95)
          :Rc-sh   (si/MPa 3.25)
          :Rc-mc1  (si/MPa 19.6)
          :Rc-mc2  (si/MPa 16.7)
          :Ec      (si/MPa 34.5e3)
          :n'      15.0
          :cn      ((/ si/MPa) 75e-6)})
         
(def B40 {:Rc      (si/MPa 20.0)  
          :Rct     (si/MPa 1.25)
          :Rc-ser  (si/MPa 29.0)
          :Rct-ser (si/MPa 2.10)
          :Rc-sh   (si/MPa 3.60)
          :Rc-mc1  (si/MPa 23.0)
          :Rc-mc2  (si/MPa 19.6)
          :Ec      (si/MPa 36.0e3)
          :n'      10.0
          :cn      ((/ si/MPa) 67e-6)})
         
(def B45 {:Rc      (si/MPa 22.0)  
          :Rct     (si/MPa 1.30)
          :Rc-ser  (si/MPa 32.0)
          :Rct-ser (si/MPa 2.20)
          :Rc-sh   (si/MPa 3.80)
          :Rc-mc1  (si/MPa 26.0)
          :Rc-mc2  (si/MPa 22.0)
          :Ec      (si/MPa 37.5e3)
          :n'      10.0
          :cn      ((/ si/MPa) 55e-6)})
         
(def B50 {:Rc      (si/MPa 25.0)  
          :Rct     (si/MPa 1.40)
          :Rc-ser  (si/MPa 36.0)
          :Rct-ser (si/MPa 2.30)
          :Rc-sh   (si/MPa 4.15)
          :Rc-mc1  (si/MPa 29.9)
          :Rc-mc2  (si/MPa 25.0)
          :Ec      (si/MPa 39.0e3)
          :n'      10.0
          :cn      ((/ si/MPa) 50e-6)})
         
(def B55 {:Rc      (si/MPa 27.5)  
          :Rct     (si/MPa 1.45)
          :Rc-ser  (si/MPa 39.5)
          :Rct-ser (si/MPa 2.40)
          :Rc-sh   (si/MPa 4.45)
          :Rc-mc1  (si/MPa 32.8)
          :Rc-mc2  (si/MPa 27.5)
          :Ec      (si/MPa 39.5e3)
          :n'      10.0
          :cn      ((/ si/MPa) 41e-6)})
         
(def B60 {:Rc      (si/MPa 30.0)  
          :Rct     (si/MPa 1.50)
          :Rc-ser  (si/MPa 43.0)
          :Rct-ser (si/MPa 2.50)
          :Rc-sh   (si/MPa 4.75)
          :Rc-mc1  (si/MPa 36.2)
          :Rc-mc2  (si/MPa 30.0)
          :Ec      (si/MPa 40.0e3)
          :n'      10.0
          :cn      ((/ si/MPa) 39e-6)})
