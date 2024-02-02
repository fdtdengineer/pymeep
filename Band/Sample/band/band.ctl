(set-param! num-bands 10) ;固有周波数の数

(define-param k-interp 4) ; number of k points to interpolate

(define Gamma (vector3 0 0 0))
(define X (vector3 0.5 0 0))
(define M (vector3 0.5 0.5 0))
(set! k-points (interpolate k-interp (list Gamma X M Gamma)))
;Γ点,X点,K点の指定


(set! geometry (list (make cylinder (center 0 0 0) (radius 0.2) (height infinity)
            (material(make dielectric (epsilon 12))))))

(set! geometry-lattice (make lattice (size 1 1 no-size))) ;単位格子サイズ

(set-param! resolution 32) #メッシュの細かさ的な

(begin-time
 "total time for both TE and TM bands: "
 (run-te)
 (run-tm))

(display-eigensolver-stats)