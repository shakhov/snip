(ns shakhov.snip.concrete-test
  (:refer-clojure :exclude [time force + - * / < > <= >= = zero? pos? neg? sgn abs
                            sin cos tan asin acos atan exp log min max])
  (:use clojure.test
        shakhov.snip.concrete
        shakhov.flow.core
        shakhov.snip.pprint)

  (:use [clojure.algo.generic.arithmetic :only [+ - * /]]
        [clojure.algo.generic.comparison :only [< > <= >= = zero? pos? neg? min max]]
        [clojure.algo.generic.math-functions :only [pow sqrt sgn abs sin cos tan
                                                    asin acos atan exp log]])
  (:require [shakhov.snip.dimensions :as dim]
            [shakhov.snip.units :as si])

  (:require [shakhov.snip.materials.concrete :as concrete]
            [shakhov.snip.materials.steel :as steel]
            [shakhov.snip.materials.rebar :as rebar]))

(defn flip-cs
  [cs]
  (update-in cs [:reinf] #(merge %
                                 {:top (:bottom %)}
                                 {:bottom (:top %)})))

(defn test-I
  [& cs]
  (let [res (mapv bending cs)]
    (println "Расчет на прочность:")
    (println (vtable {:table {:width 64}
                      :rows [{:title "Тип сечения (расч.)"
                              :key :as}
                             {:title "Высота - h, cм"
                              :key :h
                              :format (format-units si/cm "%.1f")}
                             {:title "Ширина ребра - b, cм"
                              :key :b
                              :format (format-units si/cm "%.1f")}
                             {:title "Толщина плиты - hf, cм"
                              :key :hf
                              :format (format-units si/cm "%.1f")}
                             {:title "Ширина плиты - bf, cм"
                              :key :bf
                              :format (format-units si/cm "%.1f")}
                             {:title "Рабочая высота - h0, cм"
                              :key :h0
                              :format (format-units si/cm "%.1f")}
                             {:title "Высота сжатой зоны - x, см"
                              :key :x
                              :format (format-units si/cm "%.2f")}
                             {:title "Растянутая арматура - As, см2"
                              :key :Ar
                              :format (format-units (pow si/cm 2) "%.2f")}
                             {:title "Центр тяжести - as, см"
                              :key :ar
                              :format (format-units si/cm "%.2f")}
                             {:title "Учтено сжатой арматуры - Asc, см2"
                              :key :Arc
                              :format (format-units (pow si/cm 2) "%.2f")}
                             {:title "Относительная высота сж. зоны - xi"
                              :key :xi
                              :format #(format "%.3f" %)}
                             {:title "Предельная высота сж. зоны - xi-y"
                              :key :xi-y
                              :format #(format "%.3f" %)}
                             {:title "Предельный момент - M, тм"
                              :key :M-max
                              :format (format-units si/tonf*m "%.1f")}
                             {:title "Предельный момент - Msc, тм"
                              :key :M-max-sc
                              :format (format-units si/tonf*m "%.1f")}]}
                     res))))

(defn test-II
  [& cs]
  (let [res (mapv cracking cs)]
    (println "Расчет на раскрытие трещин:")
    (println (vtable {:table {:width 64}
                      :rows [{:title "Тип сечения (расч.)"
                              :key :as}
                             {:title "Высота сжатой зоны - x, см"
                              :key :x-el
                              :format (format-units si/cm "%.2f")}
                             {:title "Растянутая арматура - As-red, см2"
                              :key :Ar-red
                              :format (format-units (pow si/cm 2) "%.2f")}
                             {:title "Сжатая арматура - Asc-red, см2"
                              :key :Arc-red
                              :format (format-units (pow si/cm 2) "%.2f")}
                             {:title "Коэффициент n'"
                              :key :n'}
                             {:title "Момент инерции - Iel, м4"
                              :key :I-red-el
                              :format (format-units (pow si/cm 4) "%.3G")}
                             ;; {:title "Напряжения в бетоне - sb, кг/cм2"
                             ;;  :key :sigma-c
                             ;;  :format (format-units si/kgf:cm2 "%.2f")}
                             {:title "Предельные напряж. - Rb,mc2, кг/cм2"
                              :key :Rc-mc2
                              :format (format-units si/kgf:cm2 "%.2f")}
                             ;; {:title "Напряжения в арматуре - sr, кг/cм2"
                             ;;  :key :sigma-r
                             ;;  :format (format-units si/kgf:cm2 "%.2f")}
                             {:title "К-т раскрытия трещин - psi, см"
                              :key :psi-cr
                              :format (format-units si/cm "%.3f")}
                             ;; {:title "Ширина раскрытия трещин a-cr, см"
                             ;;  :key :a-cr
                             ;;  :format (format-units si/cm "%.3f")}
                             {:title "Предельный момент Mmax,cr, тм"
                              :key :M-max-cr
                              :format (format-units si/tonf*m "%.1f")}
                             {:title "Предельный момент Mmax,mc2, тм"
                              :key :M-max-mc2
                              :format (format-units si/tonf*m "%.1f")}]}
                     res))))

(defn test-all
  [& cs]
  (apply test-I cs)
  (println)
  (apply test-II cs))

(comment
  (test-all
   (merge {:h (si/cm 76)
           :b (si/cm 156)
           :bottom-rebar [{:rebar rebar/AIII
                           :ar (si/mm 70)
                           :d (si/mm 18)
                           :n 10
                           :beta-cr 1.0}]
           :top-rebar [{:rebar rebar/AIII
                        :ar (si/mm 70)
                        :d (si/mm 20)
                        :n 10
                        :beta-cr 1.0}]}
          rebar/AIII concrete/B30
          {:delta-cr (si/cm 0.02)})))
