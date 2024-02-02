(reset-meep)
(define-param eps-averaging? false)
(define-param sx 40) 
(define-param sy 40) 
(set! geometry-lattice (make lattice (size  sx sy no-size)))
(define-param pad 5) ; padding distance between waveguide and cell edge
(define-param w 6) ; width of waveguide 
(define-param lam 5) ; wavelength
(define-param thk 1) ; thickness of the hollow metal waveguide
(define-param ctime (* 100 lam) ) ; computational time 10 periods
(define wvg-ycen (* 0.5 (- sy w (* 2 pad)))) ; y center of horiz. wvg
(define wvg-xcen (* 0.5 (- sx w (* 2 pad)))) ; x center of vert. wvg
(set! geometry (list 
  (make block (center (* -0.5 pad) wvg-ycen) (size (- sx pad) w infinity) (material metal))
  (make block (center  wvg-xcen (* -0.5 pad)) (size w (- sy pad) infinity) (material metal))
  (make block (center (- (* -0.5 pad) (* 0.5 thk)) wvg-ycen) (size (- sx (+ pad thk)) (- w (* 2 thk)) infinity) (material air))
  (make block (center  wvg-xcen (- (* -0.5 pad) (* 0.5 thk))) (size (- w (* 2 thk)) (- sx (+ pad thk)) infinity) (material air))
))
(set! pml-layers ( list  (make pml (thickness 1.0))))
(set! resolution 10)
(set! sources ( list (make source
  (src (make continuous-src (wavelength lam))) (component Ez) (center (+ 1 (* -0.5 sx)) wvg-ycen) (size 0 w)
)))
(run-until ctime
    (at-beginning output-epsilon)
;    (to-appended "ez" (at-every 0.6 output-efield-z))
    (at-every 1 (output-png Ez "-Zc /usr/share/h5utils/colormaps/dkbluered"))
)