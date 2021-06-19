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
;;      * redefine eval-A as a macroc
;;    Modified by Bela Pecsek/Tomas Wain
;;      * Substantial speedup compared to sbcl-9 of Shubhamkar Ayare
;;      * Using SSE calculation in two lanes
;;      * Improvement in type declarations
(declaim (optimize (speed 3) (safety 0) (space 0) (debug 0)))

(ql:quickload :sb-simd :silent t)

(defpackage #:spectralnorm3
  (:use #:cl #:sb-simd-avx) 
  (:nicknames #:sn3)
  (:import-from #:cl-user #:define-alien-routine
                          #:long
                          #:int)
  (:export #:main
           #:spectralnorm))

(in-package #:spectralnorm3)

(deftype uint31 (&optional (bits 31)) `(unsigned-byte ,bits))

(declaim (ftype (function (f64.2 f64.2) f64.2) eval-A)
         (inline eval-A))
(defun eval-A (%i %j)
  (let* ((%i+1   (f64.2+ %i (f64.2-broadcast 1d0)))
         (%i+j   (f64.2+ %i %j))
         (%i+j+1 (f64.2+ %i+1 %j)))
    (f64.2+ (f64.2* (f64.2* %i+j %i+j+1) (f64.2-broadcast 0.5d0)) %i+1)))

(declaim (ftype (function (f64vec f64vec uint31 uint31 uint31) null)
                eval-A-times-u eval-At-times-u)
         (inline eval-A-times-u eval-At-times-u))
(defun eval-A-times-u (src dst begin end length)
  (loop with %src0 of-type f64.2 = (f64.2-broadcast (aref src 0))
        with %0.0  of-type f64.2 = (f64.2-broadcast 0d0)
	for i of-type uint31 from begin below end by 4
	do (let* ((%eA0   (eval-A (make-f64.2 (+ i 0) (+ i 1)) %0.0))
		  (%eA1   (eval-A (make-f64.2 (+ i 2) (+ i 3)) %0.0))
		  (%sum0  (f64.2/ %src0 %eA0))
		  (%sum1  (f64.2/ %src0 %eA1))
		  (%ti0   (make-f64.2 (+ i 0) (+ i 1)))
		  (%ti1   (make-f64.2 (+ i 2) (+ i 3)))
		  (%last0 %eA0)
		  (%last1 %eA1))
	     (loop for j of-type uint31 from 1 below length
		   do (let* ((%j     (f64.2-broadcast j))
			     (src-j  (aref src j))
			     (%src-j (f64.2-broadcast src-j))
			     (%idx0  (f64.2+ %last0 %ti0 %j))
			     (%idx1  (f64.2+ %last1 %ti1 %j)))
			(setf %last0 %idx0
                              %last1 %idx1)
			(f64.2-incf %sum0 (f64.2/ %src-j %idx0))
			(f64.2-incf %sum1 (f64.2/ %src-j %idx1))))
	     (setf (f64.2-aref dst (+ i 0)) %sum0
                   (f64.2-aref dst (+ i 2)) %sum1))))

(defun eval-At-times-u (src dst begin end length)
  (loop with %src0 of-type f64.2 = (f64.2-broadcast (aref src 0))
        with %0.0  of-type f64.2 = (f64.2-broadcast 0d0)
	for i of-type uint31 from begin below end by 4
        do (let* ((%eAt0  (eval-A %0.0 (make-f64.2 (+ i 0) (+ i 1))))
		  (%eAt1  (eval-A %0.0 (make-f64.2 (+ i 2) (+ i 3))))
                  (%sum0  (f64.2/ %src0 %eAt0))
		  (%sum1  (f64.2/ %src0 %eAt1))
                  (%ti0   (make-f64.2 (+ i 1) (+ i 2)))
                  (%ti1   (make-f64.2 (+ i 3) (+ i 4)))
                  (%last0 %eAt0)
                  (%last1 %eAt1))
	     (loop for j of-type uint31 from 1 below length
                   do (let* ((%j     (f64.2-broadcast j))
			     (src-j  (aref src j))
			     (%src-j (f64.2-broadcast src-j))
			     (%idx0  (f64.2+ %last0 %ti0 %j))
			     (%idx1  (f64.2+ %last1 %ti1 %j)))
			(setf %last0 %idx0
                              %last1 %idx1)
			(f64.2-incf %sum0 (f64.2/ %src-j %idx0))
			(f64.2-incf %sum1 (f64.2/ %src-j %idx1))))
	     (setf (f64.2-aref dst (+ i 0)) %sum0
                   (f64.2-aref dst (+ i 2)) %sum1))))

(declaim (ftype (function () (integer 1 256)) get-thread-count))
#+sb-thread
(defun get-thread-count ()
  (progn (define-alien-routine sysconf long (name int))
         (sysconf 84)))

(declaim (ftype (function (uint31 uint31 function) null) execute-parallel))
#+sb-thread
(defun execute-parallel (start end function)
  (declare (optimize (speed 0)))
  (let* ((n (truncate (- end start) (get-thread-count)))
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

(declaim (ftype (function (f64vec f64vec f64vec uint31 uint31 uint31) null)
                EvalAtATimesU))
(defun eval-AtA-times-u (src dst tmp start end N)
      (progn
	(execute-parallel start end (lambda (start end)
				      (eval-A-times-u src tmp start end N)))
	(execute-parallel start end (lambda (start end)
				      (eval-At-times-u tmp dst start end N)))))

(declaim (ftype (function (uint31) f64) spectralnorm))
(defun spectralnorm (n)
  (declare (optimize (speed 3) (safety 0) (space 0) (debug 0)))
  (let ((u (make-array (+ n 3) :element-type 'f64 :initial-element 1.0d0))
        (v (make-array (+ n 3) :element-type 'f64))
        (tmp (make-array (+ n 3) :element-type 'f64)))
    (declare (type f64vec u v tmp))
    (loop repeat 10 do
      (eval-AtA-times-u u v tmp 0 n n)
      (eval-AtA-times-u v u tmp 0 n n))
    (sqrt (the f64 (/ (f64.4-vdot u v)
                      (f64.4-vdot v v))))))

(declaim (ftype (function (&optional uint31) null) main))
(defun main (&optional (n-supplied 5500))
  (let ((n (or n-supplied (parse-integer (second sb-ext::*posix-argv*)))))
    (declare (type uint31 n)) 
    (if (< n 8) (error "The supplied value of 'n' bust be at least 8")
        (format t "~11,9F~%" (spectralnorm N)))))
