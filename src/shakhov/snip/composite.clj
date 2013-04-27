(ns shakhov.snip.composite
  (:use [shakhov.flow.core]))

(def ^:private I-rect
  (fnk {:keys [b h dz] :or {dz 0.0}}
       (+ (* 1/12 b h h h)
          (* b h dz dz))))

(def I-cs-geometry
  (flow 
    {:nc (fnk 
           [Est Ec]
           (/ Est Ec))
     
     :Ec-creep (fnk
                 [Ec phi-creep]
                 (/ Ec (+ 1 phi-creep)))
     
     :nc-creep (fnk
                 [Est Ec-creep]
                 (/ Est Ec-creep))
     
     :psi-crack (fnk [] 0.5)
     
     ;; Ai
     
     :Atf-i (fnk
              [top-flange]
              (map (fn [{:keys [w t]}] (* w t))
                   top-flange))
     :Aw (fnk
           {{:keys [h t]} :web}
           (* h t))
     
     :Abf-i (fnk
              [bottom-flange]
              (map (fn [{:keys [w t]}] (* w t))
                   bottom-flange))
     
     :Atf (fnk
            [Atf-i]
            (apply + Atf-i))
     
     :Abf (fnk
            [Abf-i]
            (apply + Abf-i))
     
     :Asl (fnk
            [slab Ar]
            (* (:b slab)
                  (:h slab)))
     
     :Asl-red (fnk
                [Asl nc]
                (/ Asl nc))
     
     :Asl-red-creep (fnk
                      [Asl nc-creep]
                      (/ Asl nc-creep))
     
     :Ar (fnk {{:keys [n d]} :reinf}
              (* n 1/4 Math/PI d d))
     
     :Ar-stc (fnk
               [Ar nc]
               (* Ar (- 1 (/ nc))))
     
     :Ar-stc-creep (fnk
                     [Ar nc-creep]
                     (* Ar (- 1 (/ nc-creep))))
     
     :Ar-crack (fnk 
                 [Ar psi-crack]
                 (/ Ar psi-crack))
     
     ;; A
     
     :Ast (fnk
            [Atf Abf Aw]
            (+ Atf Aw Abf))
     
     
     :Astc (fnk
             [Asl-red Ast Ar-stc]
             (+ Ast Asl-red Ar-stc))
     
     :Astc-creep (fnk
                   [Asl-red-creep Ast Ar-stc-creep]
                   (+ Ast Asl-red-creep Ar-stc-creep))
     
     :Astc-crack (fnk
                   [Ast Ar-crack]
                   (+ Ast Ar-crack))
     
     ;; Zw
     
     :Zw-tf (fnk
              [top-flange web]
              (let [z (reverse (reduce (fn [z {:keys [t]}]
                                         (conj z (- (peek z) t)))
                                       [(- (/ (:h web) 2))] 
                                       (reverse top-flange)))]
                (map #(/ (+ %1 %2) 2) z (next z))))
     
     :Zw-bf (fnk
              [bottom-flange web]
              (let [z (reduce (fn [z {:keys [t]}]
                                (conj z (+ (peek z) t)))
                              [(/ (:h web) 2)] 
                              bottom-flange)]
                (map #(/ (+ %1 %2) 2) z (next z))))
     
     :Zw-sl (fnk
              [slab web top-flange]
              (- (+ (/ (:h web) 2) 
                    (apply + (map :t top-flange))
                    (/ (:h slab) 2))))
     
     :Zw-r (fnk
              [Zw-sl]
             Zw-sl)
     ;; Sw
     
     :Sw-st (fnk
              [Atf-i Abf-i Zw-tf Zw-bf]
              (apply + (map (fn [ai zi] (* ai zi)) 
                            (concat Atf-i Abf-i)
                            (concat Zw-tf Zw-bf))))
     
     :Sw-stc (fnk
               [Sw-st Asl-red Zw-sl Ar-stc Zw-r]
               (+ Sw-st
                  (* Ar-stc Zw-r)
                  (* Asl-red Zw-sl)))
     
     
     :Sw-stc-creep (fnk
                     [Sw-st Asl-red-creep Zw-sl Ar-stc-creep Zw-r]
                     (+ Sw-st
                        (* Ar-stc-creep Zw-r)
                        (* Asl-red-creep Zw-sl)))
     
     :Sw-stc-crack (fnk
                     [Sw-st Ar-crack Zw-r]
                     (+ Sw-st
                        (* Ar-crack Zw-r)))
     
     ;; dZw
     
     :dZw-st (fnk
               [Sw-st Ast]
               (/ Sw-st Ast))
     
     
     :dZw-stc (fnk
                [Sw-stc Astc]
                (/ Sw-stc Astc))
     
     :dZw-stc-creep (fnk
                      [Sw-stc-creep Astc-creep]
                      (/ Sw-stc-creep Astc-creep))
     
     :dZw-stc-crack (fnk
                      [Sw-stc-crack Astc-crack]
                      (/ Sw-stc-crack Astc-crack))
     
     ;; Iw
     
     :Iw-st (fnk
              [top-flange web bottom-flange Zw-tf Zw-bf]
              (+ (I-rect {:b (:t web) :h (:h web)})
                 (apply + (map (fn [{:keys [w t]} z]
                                 (I-rect {:b w :h t :dz z}))
                               (concat top-flange bottom-flange)
                               (concat Zw-tf Zw-bf))))) 
     
     :Iw-stc (fnk
               [Iw-st Zw-sl slab nc Ar-stc Zw-r Asl]
               (let [{:keys [b h]} slab]
                 (+ Iw-st
                    (* Ar-stc Zw-r Zw-r)
                    (/ (I-rect {:h h :b b :dz Zw-sl})
                       nc))))
     
     :Iw-stc-creep (fnk
                     [Iw-st Zw-sl slab nc-creep Ar-stc-creep Zw-r Asl]
                     (let [{:keys [b h]} slab]
                       (+ Iw-st
                          (* Ar-stc-creep Zw-r Zw-r)
                          (/ (I-rect {:h h :b b :dz Zw-sl})
                             nc-creep))))
     
     :Iw-stc-crack (fnk
                     [Iw-st Zw-r Ar-crack]
                       (+ Iw-st
                          (* Ar-crack Zw-r Zw-r)))
     ;; Ist
     
     :Ist (fnk
            [Iw-st dZw-st Ast]
            (- Iw-st (* Ast dZw-st dZw-st)))
     
     :Istc (fnk
             [Iw-stc dZw-stc Astc]
             (- Iw-stc (* Astc dZw-stc dZw-stc)))
     
     :Istc-creep (fnk
                   [Iw-stc-creep dZw-stc-creep Astc-creep]
                   (- Iw-stc-creep (* Astc-creep dZw-stc-creep dZw-stc-creep)))
     
     :Istc-crack (fnk
                   [Iw-stc-crack dZw-stc-crack Astc-crack]
                   (- Iw-stc-crack (* Astc-crack dZw-stc-crack dZw-stc-crack)))
     
     :EI-st (fnk [Ist Est] (* Ist Est))
     :EI-stc (fnk [Istc Est] (* Istc Est))
     :EI-stc-creep (fnk [Istc-creep Est] (* Istc-creep Est))
     :EI-stc-crack (fnk [Istc-crack Est] (* Istc-crack Est))
     
     ;; Z
     
     :Zs2-st (fnk 
               [top-flange web dZw-st]
               (- (- (apply + (/ (:h web) 2)
                            (map :t top-flange)))
                  dZw-st))
     
     :Zs1-st (fnk 
               [bottom-flange web dZw-st]
               (- (apply + (/ (:h web) 2)
                         (map :t bottom-flange))
                  dZw-st))
     
     :Zs2-stc (fnk 
                [top-flange web dZw-stc]
                (- (- (apply + (/ (:h web) 2)
                             (map :t top-flange)))
                   dZw-stc))
     
     :Zs1-stc(fnk 
               [bottom-flange web dZw-stc]
               (- (apply + (/ (:h web) 2)
                         (map :t bottom-flange))
                  dZw-stc))
     
     :Zs2-stc-crack (fnk 
                      [top-flange web dZw-stc-crack]
                      (- (- (apply + (/ (:h web) 2)
                                   (map :t top-flange)))
                         dZw-stc-crack))
     
     :Zs1-stc-crack (fnk 
                      [bottom-flange web dZw-stc-crack]
                      (- (apply + (/ (:h web) 2)
                                (map :t bottom-flange))
                         dZw-stc-crack))
     
     :Zc-st  (fnk 
               [Zw-sl dZw-st]
               (- Zw-sl dZw-st))
     
     :Zc-stc (fnk
               [Zw-sl dZw-stc]
               (- Zw-sl dZw-stc))
     
     :Zc-stc-creep  (fnk 
                       [Zw-sl dZw-stc-creep]
                       (- Zw-sl dZw-stc-creep))
     
     :Zr-stc-crack (fnk
                      [Zw-r dZw-stc-crack]
                      (- Zw-r dZw-stc-crack))
     
     ;; W
     
     :Wc-st (fnk
               [Zc-st Ist]
               (/ Ist Zc-st))
     
     :Ws1-st (fnk
               [Zs1-st Ist]
               (/ Ist Zs1-st))
     
     :Ws2-st (fnk
               [Zs2-st Ist]
               (/ Ist Zs2-st))
     
     :Wc-stc (fnk
               [Zc-stc Istc]
                (/ Istc Zc-stc))
     
     :Ws1-stc (fnk
               [Zs1-stc Istc]
                (/ Istc Zs1-stc))
     
     :Ws2-stc (fnk
               [Zs2-stc Istc]
               (/ Istc Zs2-stc))
     
     :Wr-stc-crack (fnk
                      [Zr-stc-crack Istc-crack]
                      (/ Istc-crack Zr-stc-crack))
     
     :Ws1-stc-crack (fnk
                      [Zs1-stc-crack Istc-crack]
                      (/ Istc-crack Zs1-stc-crack))
     
     :Ws2-stc-crack (fnk
               [Zs2-stc-crack Istc-crack]
               (/ Istc-crack Zs2-stc-crack))}))






