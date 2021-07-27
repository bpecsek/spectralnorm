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
;;      * Using SSE calculation
;;      * Improvement in type declarations
;;      * Changed code to be compatible with sb-simd
;;      * Eliminated mixing VEX and non-VEX instructions as far as possible
;;        in the hot loops
(declaim (optimize (speed 3) (safety 0) (debug 0)))

(asdf:load-system :sb-simd)

(defpackage #:spectralnorm1
  (:use #:cl #:sb-simd-avx) 
  (:nicknames #:sn1)
  (:import-from #:cl-user #:define-alien-routine
                          #:long
                          #:int)
  (:export #:main))

(in-package #:spectralnorm1)

(declaim (ftype (function (f64.2 f64.2) f64.2) eval-A))
(define-inline eval-A (%i %j)
  (let* ((%i+1   (f64.2+ %i (f64.2 1)))
         (%i+j   (f64.2+ %i %j))
         (%i+j+1 (f64.2+ %i+1 %j)))
    (f64.2+ (f64.2* %i+j %i+j+1 (f64.2 0.5)) %i+1)))

(declaim (ftype (function (f64vec f64vec u32 u32 u32) null)
                eval-A-times-u eval-At-times-u))
(defun eval-A-times-u (src dst begin end length)
  (loop for i of-type u32 from begin below end by 2
        do (let* ((%src-0 (f64.2 (aref src 0)))
                  (%ti  (f64.2+ (f64.2 i) (make-f64.2 0 1)))
                  (%eA  (eval-A %ti (f64.2 0)))
		  (%sum (f64.2/ %src-0 %eA)))
	     (loop for j of-type u32 from 1 below length
		   do (let* ((src-j (aref src j))
                             (%idx (f64.2+ %eA %ti (f64.2 j))))
			(setf %eA %idx)
			(f64.2-incf %sum (f64.2/ (f64.2 src-j) %idx))))
	     (setf (f64.2-aref dst i) %sum))))

(defun eval-At-times-u (src dst begin end length)
  (loop for i of-type u32 from begin below end by 2
        do (let* ((%src-0 (f64.2 (aref src 0)))
                  (%ti  (f64.2+ (f64.2 i) (make-f64.2 1 2)))
                  (%eAt (eval-A (f64.2 0) (f64.2- %ti)))
                  (%sum (f64.2/ %src-0 %eAt)))
	     (loop for j of-type u32 from 1 below length
                   do (let* ((src-j (aref src j))
                             (%idx (f64.2+ %eAt %ti (f64.2 j))))
			(setf %eAt %idx)
			(f64.2-incf %sum (f64.2/ (f64.2 src-j) %idx))))
	     (setf (f64.2-aref dst i) %sum))))

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
                eval-AtA-times-u))
(defun eval-AtA-times-u (src dst tmp start end N)
      (progn
	(execute-parallel start end (lambda (start end)
				      (eval-A-times-u src tmp start end N)))
	(execute-parallel start end (lambda (start end)
				      (eval-At-times-u tmp dst start end N)))))

(declaim (ftype (function (u32) null) spectralnorm))
(defun spectralnorm (n)
  (let ((u   (make-array (1+ n) :element-type 'f64 :initial-element 1d0))
        (v   (make-array (1+ n) :element-type 'f64))
        (tmp (make-array (1+ n) :element-type 'f64)))
    (declare (type f64vec u v tmp))
    (loop repeat 10 do
      (eval-AtA-times-u u v tmp 0 N N)
      (eval-AtA-times-u v u tmp 0 N N))
    (let ((sumvb 0d0)
          (sumvv 0d0))
      (loop for i below n
            for aref-v-i of-type f64 = (aref v i)
            do (incf sumvb (the f64 (* (aref u i) aref-v-i)))
               (incf sumvv (the f64 (* aref-v-i aref-v-i))))
      (format t "~11,9F~%" (sqrt (the f64 (/ sumvb sumvv)))))))

(declaim (ftype (function (&optional u32) null) main))
(defun main (&optional (n-supplied 5500))
  (let ((n (or n-supplied (parse-integer (second sb-ext::*posix-argv*)))))
    (if (< n 8)
        (error "The supplied value of 'n' must be at least 8")
        (spectralnorm n))))
