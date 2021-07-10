;;    The Computer Language Benchmarks Game
;;    https://salsa.debian.org/benchmarksgame-team/benchmarksgame/
;;
;;    Adapted from the C (gcc) code by Sebastien Loisel
;;
;;    Contributed by Christopher Neufeld
;;    Modified by Juho Snellman 2005-10-26
;;      * Use SIMPLE-ARRAY instead of ARRAY in declarations
;;      * Use TRUNCATE instead of / for fixnum division
;;      * Rearrange EVAL-A to make it more readable and a bit faster
;;    Modified by Andy Hefner 2008-09-18
;;      * Eliminate array consing
;;      * Clean up type declarations in eval-A
;;      * Distribute work across multiple cores on SBCL
;;    Modified by Witali Kusnezow 2008-12-02
;;      * use right shift instead of truncate for division in eval-A
;;      * redefine eval-A as a macro
;;    Modified by Bela Pecsek
;;      * Using AVX calculations
;;      * Improvement in type declarations
;;      * Changed code to be compatible with sb-simd
;;      * Eliminated mixing VEX and non-VEX instructions as far as possible
(declaim (optimize (speed 3) (safety 0) (space 0) (debug 0)))

(asdf:load-system :sb-simd)

(defpackage #:spectralnorm2
  (:use #:cl :sb-simd-avx) 
  (:nicknames #:sn2)
  (:import-from #:cl-user #:define-alien-routine
                          #:long
                          #:int)
  (:export #:main))

(in-package #:spectralnorm2)

(declaim (ftype (function (f64.4 f64.4) f64.4) eval-A))
(define-inline eval-A (%i %j)
  (let* ((%i+1   (f64.4+ %i (f64.4 1d0)))
         (%i+j   (f64.4+ %i %j))
         (%i+j+1 (f64.4+ %i+1 %j)))
    (f64.4+ (f64.4* %i+j %i+j+1 (f64.4 0.5d0)) %i+1)))

(declaim (ftype (function (f64vec f64vec u32 u32 u32) null) eval-A-times-u))
(defun eval-A-times-u (src dst begin end length)
  (loop for i from begin below end by 4
        with %src-0 of-type f64.4 = (f64.4 (aref src 0))
        with %i+    of-type f64.4 = (make-f64.4 0d0 1d0 2d0 3d0)
        do (let* ((%ti  (f64.4+ (f64.4 i) %i+))
                  (%eA  (eval-A %ti (f64.4 0d0)))
		  (%sum (f64.4/ %src-0 %eA)))
	     (loop for j from 1 below length
                   for src-j of-type f64 = (aref src j)
		   do (let ((%idx (f64.4+ %eA %ti (f64.4 j))))
			(setf %eA %idx)
			(f64.4-incf %sum (f64.4/ (f64.4 src-j) %idx))))
	     (setf (f64.4-aref dst i) %sum))))

(declaim (ftype (function (f64vec f64vec u32 u32 u32) null) eval-At-times-u))
(defun eval-At-times-u (src dst begin end length)
  (loop for i from begin below end by 4
        with %src-0 of-type f64.4 = (f64.4 (aref src 0))
        with %i++   of-type f64.4 = (make-f64.4 1d0 2d0 3d0 4d0)
        do (let* ((%ti  (f64.4+ (f64.4 i) %i++))
                  (%eAt (eval-A (f64.4 0d0) (f64.4- %ti (f64.4 1))))
		  (%sum (f64.4/ %src-0 %eAt)))
	     (loop for j of-type u32 from 1 below length
                   for src-j of-type f64 = (aref src j)
                   do (let ((%idx (f64.4+ %eAt %ti (f64.4 j))))
			(setf %eAt %idx)
			(f64.4-incf %sum (f64.4/ (f64.4 src-j) %idx))))
	     (setf (f64.4-aref dst i) %sum))))

(declaim (ftype (function () (integer 1 256)) get-thread-count))
#+sb-thread
(defun get-thread-count ()
  (progn (define-alien-routine sysconf long (name int))
         (sysconf 84)))

(declaim (ftype (function (u32 u32 function) null) execute-parallel))
#+sb-thread
(defun execute-parallel (start end function)
  (declare (optimize (speed 0)))
  (let* ((n    (truncate (- end start) (get-thread-count)))
         (step (- n (mod n 2))))
    (loop for i from start below end by step
          collecting (let ((start i)
                           (end (min end (+ i step))))
                       (sb-thread:make-thread
			(lambda () (funcall function start end))))
            into threads
          finally (mapcar #'sb-thread:join-thread threads))))

#-sb-thread
(defun execute-parallel (start end function)
  (funcall function start end))

(declaim (ftype (function (f64vec f64vec f64vec u32 u32 u32) null)
                EvalAtATimesU))
(defun eval-AtA-times-u (src dst tmp start end N)
      (progn
	(execute-parallel start end (lambda (start end)
				      (eval-A-times-u src tmp start end N)))
	(execute-parallel start end (lambda (start end)
				      (eval-At-times-u tmp dst start end N)))))

(declaim (ftype (function (u32) f64) spectralnorm))
(defun spectralnorm (n)
  (let ((u   (make-array (+ n 3) :element-type 'f64 :initial-element 1.0d0))
        (v   (make-array (+ n 3) :element-type 'f64))
        (tmp (make-array (+ n 3) :element-type 'f64)))
    (declare (type f64vec u v tmp))
    (loop repeat 10 do
      (eval-AtA-times-u u v tmp 0 n n)
      (eval-AtA-times-u v u tmp 0 n n))
    (sqrt (/ (f64.4-vdot u v) (f64.4-vdot v v)))))

(declaim (ftype (function (&optional u32) null) main))
(defun main (&optional (n-supplied 5500))
  (let ((n (or n-supplied (parse-integer (second sb-ext::*posix-argv*)))))
    (if (< n 8)
        (error "The supplied value of 'n' must be at least 8"))
    (format t "~11,9F~%" (spectralnorm n))))
