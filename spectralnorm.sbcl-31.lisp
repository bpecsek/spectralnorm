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
;;    Modified by Bela Pecsek
;;      * Using SSE calculation on two lanes
;;      * Improvement in type declarations
;;      * Changed code to be compatible with sb-simd
;;      * Eliminated mixing VEX and non-VEX instructions as far as possible
;;        in the hot loops
(declaim (optimize (speed 3) (safety 0) (space 0) (debug 0)))

(asdf:load-system :sb-simd)

(defpackage #:spectralnorm31
  (:use #:cl #:sb-simd-avx) 
  (:nicknames #:sn31)
  (:import-from #:cl-user #:define-alien-routine
                          #:long
                          #:int)
  (:export #:main))

(in-package #:spectralnorm31)

(defmacro eval-A (%i %j)
  `(let* ((%i+1   (f64.4+ ,%i (f64.4 1)))
          (%i+j   (f64.4+ ,%i ,%j))
          (%i+j+1 (f64.4+ %i+1 ,%j))
          (evala  (f64.4+ (f64.4* %i+j %i+j+1 (f64.4 0.5)) %i+1)))
     (values (f64.4-extract128 evala 0)
             (f64.4-extract128 evala 1))))

(declaim (ftype (function (f64vec f64vec u32 u32 u32) null)
                eval-A-times-u eval-At-times-u))
(defun eval-A-times-u (src dst begin end length)
  (loop for i of-type u32 from begin below end by 4
	with %src-0 of-type f64.2 = (f64.2 (aref src 0))
        with %inc of-type f64.4 = (make-f64.4 0 1 2 3)
	do (multiple-value-bind (%eA0 %eA1)
               (eval-A (f64.4+ (f64.4 i) %inc) (f64.4 0))
             (let* ((%ti0   (make-f64.2 (+ i 0) (+ i 1)))
		    (%ti1   (make-f64.2 (+ i 2) (+ i 3)))
		    (%sum0  (f64.2/ %src-0 %eA0))
		    (%sum1  (f64.2/ %src-0 %eA1)))
	       (loop for j of-type u32 from 1 below length
                     for src-j of-type f64 = (aref src j)
		     do (let* ((%j     (f64.2 j))
                               (%src-j (f64.2 src-j))
			       (%idx0  (f64.2+ %eA0 %ti0 %j))
			       (%idx1  (f64.2+ %eA1 %ti1 %j)))
			  (setf %eA0 %idx0
                                %eA1 %idx1)
			  (f64.2-incf %sum0 (f64.2/ %src-j %idx0))
			  (f64.2-incf %sum1 (f64.2/ %src-j %idx1))))
	       (setf (f64.2-aref dst (+ i 0)) %sum0
                     (f64.2-aref dst (+ i 2)) %sum1)))))

(defun eval-At-times-u (src dst begin end length)
  (loop for i of-type u32 from begin below end by 4
        with %src-0 of-type f64.2 = (f64.2 (aref src 0))
        with %inc of-type f64.4 = (make-f64.4 0 1 2 3)
	do (multiple-value-bind (%eAt0 %eAt1)
               (eval-A (f64.4 0) (f64.4+ (f64.4 i) %inc))
             (let* ((%ti0   (make-f64.2 (+ i 1) (+ i 2)))
                    (%ti1   (make-f64.2 (+ i 3) (+ i 4)))
                    (%sum0  (f64.2/ %src-0 %eAt0))
		    (%sum1  (f64.2/ %src-0 %eAt1)))
	       (loop for j of-type u32 from 1 below length
                     for src-j of-type f64 = (aref src j)
                     do (let* ((%j     (f64.2 j))
                               (%src-j (f64.2 src-j))
			       (%idx0  (f64.2+ %eAt0 %ti0 %j))
			       (%idx1  (f64.2+ %eAt1 %ti1 %j)))
			  (setf %eAt0 %idx0
                                %eAt1 %idx1)
			  (f64.2-incf %sum0 (f64.2/ %src-j %idx0))
			  (f64.2-incf %sum1 (f64.2/ %src-j %idx1))))
	       (setf (f64.2-aref dst (+ i 0)) %sum0
                     (f64.2-aref dst (+ i 2)) %sum1)))))

(declaim (ftype (function () (integer 1 256)) get-thread-count))
#+sb-thread
(defun get-thread-count ()
  (progn (define-alien-routine sysconf long (name int))
         (sysconf 84)))

(declaim (ftype (function (u32 u32 function) null) execute-parallel))
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

(declaim (ftype (function (f64vec f64vec f64vec u32 u32 u32) null)
                eval-AtA-times-u))
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
    (sqrt (/ (f64.2-vdot u v) (f64.2-vdot v v)))))

(declaim (ftype (function (&optional u32) null) main))
(defun main (&optional (n-supplied 5500))
  (let ((n (or n-supplied (parse-integer (second sb-ext::*posix-argv*)))))
    (declare (type u32 n)) 
    (if (< n 8) (error "The supplied value of 'n' must be at least 8")
        (format t "~11,9F~%" (spectralnorm N)))))