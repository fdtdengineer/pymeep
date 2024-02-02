(set! geometry-lattice (make lattice (size 1 1 no-size)
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))

(define-param kz 0) ; use non-zero kz to consider vertical propagation

(set! k-points (list (vector3 0 0 kz)          ; Gamma
                     (vector3 (/ -3) (/ 3) kz) ; K
                     (vector3 0 0.5 kz)        ; M
                     (vector3 0 0 kz)))        ; Gamma
(define-param k-interp 4)
(set! k-points (interpolate k-interp k-points))

; Now, define the geometry, etcetera:

(define-param eps 11.7) ; the dielectric constant of the background
(define-param r 0.3) ; the hole radius

(set! default-material (make dielectric (epsilon eps)))
(set! geometry (list (make cylinder (center 0) (material air)
			   (radius r) (height infinity))))

(set-param! resolution 32)
(set-param! num-bands 8)

(if (= kz 0)
    (begin
      (run-te)
      (run-tm))
    (run)) ; if kz != 0 there are no purely te and tm bands