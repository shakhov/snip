(ns shakhov.snip.dimensions
  (:refer-clojure :exclude [time force])
  (:use shakhov.units.core))

(def-dimension-system civil
   length mass time temperature)

(def-dimension 
  civil
  angle [length 1 length -1]
  area [length 2]
  volume [length 3]
  density [mass 1 volume -1]
  surface-density [mass 1 area -1]
  speed [length 1 time -1]
  acceleration [speed 1 time -1]
  force [mass 1 acceleration 1]
;  pressure [force 1 area -1]
  stress [force 1 area -1]
  frequency [time -1]
  energy [force 1 length 1]
  moment [force 1 length 1]
  area-static-moment [length 3]
  area-moment-of-inertia [length 4]
  section-modulus [area-moment-of-inertia 1 length -1]
  elastic-modulus [stress 1]
  linear-load-intensity [force 1 length -1]
  area-load-intensity [force 1 area -1]
  volume-load-intensity [force 1 volume -1])