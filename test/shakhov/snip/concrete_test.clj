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

(defn circle-area
  [d]
  (* 1/4 Math/PI d d))

(comment (binding [*logger* {:pre (fn [f k i] (prn k))
                             :post (fn [f k i o] (prn k "=" o))}]))

(let [pile-cap-long-bottom {:h (si/cm 100.0)
                            :b (si/cm 100.0)
                            :d (si/mm 18)
                            :nd 6.0
                            :reinf {:bottom [{:rebar rebar/AIII
                                              :z (si/mm 120)
                                              :d (si/mm 18)
                                              :n 6}
                                             {:rebar rebar/AIII
                                              :z (si/mm 70)
                                              :d (si/mm 18)
                                              :n 10}]}
                            :Arc (* 6.0 (circle-area (si/mm 20)))
                            :arc (si/mm 70.0)
                            :M-ser (si/tonf*m 24.92)}
      csI (mapv #(rect-bending (merge % concrete/B30 rebar/AIII))
                [pile-cap-long-bottom
                 ])
      ;; csII (mapv #((lazy-compile rect-crack-width-flow)
      ;;              (merge % {:beta-cr 1.0} concrete/B30 rebar/AIII))
      ;;            [pile-cap-long-bottom])
      ]
  (defn pile-cap-I
    []
    (println "Расчет на прочность:")
    (println (vtable {:table {:width 64}
                      :cols [{:title "Вдоль, низ"
                              :key #(nth % 1)}]
                      :rows [{:title "Рабочая высота - h0, cм"
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
                             {:title "Сжатая арматура Asc, см2"
                              :key :Arc
                              :format (format-units (pow si/cm 2) "%.2f")}
                             {:title "Учтено сжатой арматуры - Asc, см2"
                              :key :Arc-ef
                              :format (format-units (pow si/cm 2) "%.2f")}
                             {:title "Относительная высота сж. зоны - xi"
                              :key :xi
                              :format #(format "%.3f" %)}
                             {:title "Предельная высота сж. зоны - xi-y"
                              :key :xi-y
                              :format #(format "%.3f" %)}
                             {:title "Предельный момент - M, тм"
                              :key :M-max
                              :format (format-units si/tonf*m "%.1f")}]}
                     csI)))
  ;; (defn pile-cap-II
  ;;   []
  ;;   (println "Расчет на раскрытие трещин:")
  ;;   (println (vtable {:table {:width 64}
  ;;                     :cols [{:title "Вдоль, низ"
  ;;                            :key #(nth % 1)}
  ;;                           {:title "Вдоль, верх"
  ;;                            :key #(nth % 2)}
  ;;                           {:title "Поперек, низ"
  ;;                            :key #(nth % 3)}
  ;;                           {:title "Поперек, верх"
  ;;                            :key #(nth % 4)}]
  ;;                    :rows [{:title "Высота сжатой зоны - x, см"
  ;;                            :key :x-el
  ;;                            :format (format-units si/cm "%.2f")}
  ;;                           {:title "Растянутая арматура - As, см2"
  ;;                            :key :Ar
  ;;                            :format (format-units (pow si/cm 2) "%.2f")}
  ;;                           {:title "Сжатая арматура Asc, см2"
  ;;                            :key :Arc
  ;;                            :format (format-units (pow si/cm 2) "%.2f")}
  ;;                           {:title "Момент инерции - Iel, м4"
  ;;                            :key :I-red-el
  ;;                            :format (format-units (pow si/cm 4) "%.3G")}
  ;;                           {:title "Напряжения в бетоне - sb, кг/cм2"
  ;;                            :key :sigma-b
  ;;                            :format (format-units si/kgf:cm2 "%.2f")}
  ;;                           {:title "Предельные напряж. - Rb,mc2, кг/cм2"
  ;;                            :key :Rc-mc2
  ;;                            :format (format-units si/kgf:cm2 "%.2f")}
  ;;                           {:title "Напряжения в арматуре - sr, кг/cм2"
  ;;                            :key :sigma-r
  ;;                            :format (format-units si/kgf:cm2 "%.2f")}
  ;;                           {:title "Ширина раскрытия трещин a-cr, см"
  ;;                            :key :a-cr
  ;;                            :format (format-units si/cm "%.3f")}]}
  ;;                    csII)))
  )
