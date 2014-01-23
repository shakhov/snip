(ns shakhov.snip.pprint
  (:use [clojure.string :only [join]]
        [shakhov.units.core :only [magnitude]]))

(defn cell [c & args]
  (let [{:keys [width align fill]} args
        s (str c)
        s (if(<= (count s) width)
                 s
                 (str (subs s 0 (- width 3)) "..."))
        fill #(apply str (repeat % fill))
        f (- width (count s))
        s (case align
            :left   (str s (fill f))
            :right  (str (fill f) s)
            :center (str (fill (Math/ceil (/ f 2)))
                         s
                         (fill (Math/floor (/ f 2)))))]
    s))

(defn col-width
  [cols data table-width]
  (let [minwidth 3
        width* (fn [{:keys [key title format] :or {format identity}}]
                 (apply max minwidth
                        (count (str title))
                        (map (comp count str format)
                             (map key data))))
        cols (map #(merge % {:minwidth (or (:width %) minwidth)
                             :maxwidth (or (:width %) (width* %))})
                  cols)
        Wmin (apply + (map :minwidth cols))
        Wmax (apply + (map :maxwidth cols))]
    (cond
     (nil? table-width)    (map #(assoc % :width (or (:width %) (:maxwidth %))) cols)
     (<= table-width Wmin) (map #(assoc % :width (:minwidth %)) cols)
     :else (let [D (/ (- table-width Wmin)
                      (- Wmax Wmin))
                 spare (atom 0)]
             (map (fn [c]
                    (assoc c :width
                           (if-let [w (:width c)]
                             (max w minwidth)
                             (let [dd (* D (- (:maxwidth c) (:minwidth c)))
                                   di (int dd)
                                   sp (if (<= 1 (swap! spare + (- dd di)))
                                        (do (swap! spare dec) 1)
                                        0)]
                               (+ (:minwidth c) di sp)))))
                  cols)))))

(defn table*
  ([spec data]
    (let [{:keys [table cols]} spec
          cols (map #(assoc % :title (or (:title %) (str (name (:key %))))) cols)
          cols (col-width cols data (:width table))
          header (for [{:keys [title width head-align head-fill]} cols]
                   (cell title
                         :width  width
                         :align (or head-align :center)
                         :fill  (or head-fill " ")))
          body   (for [row data]
                   (map (fn [{:keys [key width format align body-fill]
                              :or {format identity}}]
                          (cell (format (key row))
                                :width width
                                :align  (or align :center)
                                :fill   (or body-fill " ")))
                        cols))
          spacer (str "+-" (join "-+-" (map #(apply str (repeat (:width %) "-")) cols)) "-+")]
      (concat
        [spacer
         (str "| " (apply str (join " | " header)) " |")
         spacer]
        (for [tr body]
          (str "| " (apply str (join " | " tr)) " |"))
        [spacer]))))

(defn table
  ([data]
    (table {:cols (vec (map #(assoc {} :key %)
                             (into (sorted-set-by #(compare (str %1) (str %2)))
                                   (mapcat keys data))))}
           data))
  ([spec data]
    (apply str (join "\n" (table* spec data)))))

(defn vtable
  ([data]
    (vtable {:rows (vec (map #(assoc {} :key %)
                             (into (sorted-set-by #(compare (str %1) (str %2)))
                                   (mapcat keys data))))}
             data))
  ([spec data]
     (let [vspec {:cols (concat [{:key #(nth % 0) :title "" :align :left}]
                                (vec (or (:cols spec)
                                         (map-indexed (fn [i r] {:key #(nth % (inc i))
                                                                :title (str (inc i))})
                                                      data))))}
           vdata (for [{:keys [key format title] :or {format identity}} (:rows spec)]
                   (conj (map #(format (get % key)) data)
                         (str (name (or title key)))))]
       (table (merge vspec (dissoc spec :rows :cols)) vdata))))

(defn format-units
  [u f]
  (fn [x]
    (if x
      (format f (double (magnitude (u x))))
      "N/A")))
