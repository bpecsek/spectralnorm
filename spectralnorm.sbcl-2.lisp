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
;;    Modified by Bela Pecsek/Tomas Wain
;;      * Substantial speedup compared to sbcl-9 of Shubhamkar Ayare
;;      * Using SSE calculation in two lanes
;;      * Improvement in type declarations
(declaim (optimize (speed 3) (safety 0) (space 0) (debug 0)))

(ql:quickload :sb-simd :silent t)
(use-package '(:sb-simd))

(deftype int31 (&optional (bits 31)) `(signed-byte ,bits))
(deftype d+array () '(simple-array (double-float 0d0) (*)))
(deftype d+ () '(double-float 0d0))

(defmacro eval-A (%i %j)
  `(let* ((%i+1   (f64.4+ ,%i (make-f64.4 1 1 1 1)))
          (%i+j   (f64.4+ ,%i ,%j))
          (%i+j+1 (f64.4+ %i+1 ,%j)))
     (f64.4+ (f64.4/ (f64.4* %i+j %i+j+1) (make-f64.4 2 2 2 2)) %i+1)))

(declaim (ftype (function (d+array d+array int31 int31 int31) null) eval-A-times-u eval-At-times-u))
(defun eval-A-times-u (src dst begin end length)
  (loop with %src0 of-type f64.4 = (make-f64.4 (aref src 0) (aref src 0) (aref src 0) (aref src 0))
	for i from begin below end by 4
	do (let* ((%eAt  (eval-A (make-f64.4 (+ i 0) (+ i 1) (+ i 2) (+ i 3)) (f64.4-zeros)))
		  (%sum  (f64.4/ %src0 %eAt))
		  (%ti   (make-f64.4 (+ i 0) (+ i 1) (+ i 2) (+ i 3)))
		  (%last %eAt))
	     (loop for j from 1 below length
		   do (let* ((%j     (make-f64.4 j j j j))
			     (src-j  (aref src j))
			     (%src-j (make-f64.4 src-j src-j src-j src-j))
			     (%idx   (f64.4+ %last %ti %j)))
			(setf %last %idx)
			(f64.4-incf %sum (f64.4/ %src-j %idx))))
	     (setf (f64.4-ref dst i) %sum))))

(defun eval-At-times-u (src dst begin end length)
  (loop with %src0 of-type f64.4 = (make-f64.4 (aref src 0) (aref src 0) (aref src 0) (aref src 0))
	for i from begin below end by 4
        do (let* ((%eA   (eval-A (f64.4-zeros) (make-f64.4 (+ i 0) (+ i 1) (+ i 2) (+ i 3))))
		  (%sum  (f64.4/ %src0 %eA))
		  (%ti   (make-f64.4 (+ i 1) (+ i 2) (+ i 3) (+ i 4)))
		  (%last %eA))
	     (loop for j from 1 below length
                   do (let* ((%j     (make-f64.4 j j j j))
			     (src-j  (aref src j))
			     (%src-j (make-f64.4 src-j src-j src-j src-j))
			     (%idx   (f64.4+ %last %ti %j)))
			(setf %last %idx)
			(f64.4-incf %sum (f64.4/ %src-j %idx))))
	     (setf (f64.4-ref dst i) %sum))))

(declaim (ftype (function () (integer 1 256)) GetThreadCount))
#+sb-thread
(defun get-thread-count ()
  (progn (define-alien-routine sysconf long (name int))
         (sysconf 84)))

(declaim (ftype (function (int31 int31 function) null) execute-parallel))
#+sb-thread
(defun execute-parallel (start end function)
  (declare (optimize (speed 0)))
  (let ((step (multiple-value-bind (n _)(truncate (- end start) (get-thread-count))
		(declare (ignore _)) (- n (mod n 2)))))
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

(declaim (ftype (function (d+array d+array d+array int31 int31 int31) null)
                EvalAtATimesU))
(defun eval-AtA-times-u (src dst tmp start end N)
      (progn
	(execute-parallel start end (lambda (start end)
				      (eval-A-times-u src tmp start end N)))
	(execute-parallel start end (lambda (start end)
				      (eval-At-times-u tmp dst start end N)))))

(declaim (ftype (function (int31) d+) spectralnorm))
(defun spectralnorm (n)
  (let ((u (make-array (+ n 3) :element-type 'd+ :initial-element 1.0d0))
        (v (make-array (+ n 3) :element-type 'd+))
        (tmp (make-array (+ n 3) :element-type 'd+)))
    (declare (type d+array u v tmp))
    (loop repeat 10 do
      (eval-AtA-times-u u v tmp 0 n n)
      (eval-AtA-times-u v u tmp 0 n n))
    (let ((sumvb 0d0)
          (sumvv 0d0))
      (loop for i below n
            for aref-v-i of-type d+ = (aref v i)
            do (incf sumvb (the d+ (* (aref u i) aref-v-i)))
               (incf sumvv (the d+ (* aref-v-i aref-v-i))))
      (sqrt (the d+ (/ sumvb sumvv))))))

(declaim (ftype (function (&optional int31) null) main))
(defun main (&optional n-supplied)
  (let ((n (or n-supplied (parse-integer (or (second sb-ext::*posix-argv*)
                                             "5500")))))
    (declare (type int31 n)) 
    (if (< n 8)
        (error "The supplied value of 'n' bust be at least 8"))
    (format t "~11,9F~%" (spectralnorm n))))
