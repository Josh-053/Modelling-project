;; 1. Based on: run01
;; 2. Description: PMX001 1CMT LINEAR M7+
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 One Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX=DROP AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA      ../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step
$MODEL      NCOMP=1 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
$PK
;;; VALB-DEFINITION START
VALB = ( 1 + THETA(8)*(ALB - 31.85))
;;; VALB-DEFINITION END


;;; VCREAT-DEFINITION START
VCREAT = ( 1 + THETA(7)*(CREAT - 1.03))
;;; VCREAT-DEFINITION END


;;; VBW-DEFINITION START
IF(BW.LE.65.93) VBW = ( 1 + THETA(5)*(BW - 65.93))
IF(BW.GT.65.93) VBW = ( 1 + THETA(6)*(BW - 65.93))
;;; VBW-DEFINITION END


;;; VADA-DEFINITION START
IF(ADA.EQ.0) VADA = 1  ; Most common
IF(ADA.EQ.1) VADA = ( 1 + THETA(4))
;;; VADA-DEFINITION END

;;; V-RELATION START
VCOV=VADA*VBW*VCREAT*VALB
;;; V-RELATION END


;;; CLADA-DEFINITION START
IF(ADA.EQ.0) CLADA = 1  ; Most common
IF(ADA.EQ.1) CLADA = ( 1 + THETA(3))
;;; CLADA-DEFINITION END

;;; CL-RELATION START
CLCOV=CLADA
;;; CL-RELATION END


TVCL = THETA(1)

TVCL = CLCOV*TVCL
  CL = TVCL * EXP(ETA(1))
 TVV = THETA(2)

TVV = VCOV*TVV
   V = TVV  * EXP(ETA(2))
 K10 = CL/V
  S1 = V ; scale prediction based on DOSE (mmol) and DV (mmol/L)

$THETA  (0.001,1.38084,8) ; [1]   CL (L/d)
 (0.001,7.67529,10) ; [2,3] V  (L)
 ; values are determined in 3 iterations
$THETA  (-1,4.15573,5) ; CLADA1
$THETA  (-1,-0.279594,5) ; VADA1
$THETA  (-100000,0.0144423,0.042) ; VBW1
 (-0.055,-0.0497278,100000) ; VBW2
$THETA  (-1.250,0.264683,1.190) ; VCREAT1
$THETA  (-0.134,0.019778,0.084) ; VALB1
$DES DADT(1) = -K10*A(1)  ; ODE for central compartment

$ERROR
IPRED = F
LLOQ = 1  ; (mg/L)

; --- Residual error SDs ---
PROP_SD = SIGMA(1)
ADD_SD  = SIGMA(2)

; --- BLQ inflation (gentle, stable) ---
IF (BLQ.EQ.1) ADD_SD = ADD_SD + LLOQ

; --- Final residual error ---
W = SQRT(ADD_SD**2+PROP_SD**2)

IF (W.LE.0.000001) W=0.000001 ; Protective code

IRES=DV-IPRED
IWRES=IRES/W
Y = IPRED * (1 + EPS(1)) + EPS(2)

$OMEGA  0.598566  ;     IIV CL
 0.342736  ;      IIV V
$SIGMA  0.282155  ; EPS(1), proportional
 5E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
 ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            SIG=3 SIGL=3 PRINT=5 MAXEVALS=220 CTYPE=4

