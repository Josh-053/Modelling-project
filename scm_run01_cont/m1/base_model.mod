;; 1. Based on: run01
;; 2. Description: PMX001 1CMT LINEAR M7+
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 One Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA      ../../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step
$MODEL      NCOMP=1 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
$PK
TVCL = THETA(1)
  CL = TVCL * EXP(ETA(1))
 TVV = THETA(2)
   V = TVV  * EXP(ETA(2))
 K10 = CL/V
  S1 = V ; scale prediction based on DOSE (mmol) and DV (mmol/L)

$THETA  (0.001,1.46,8) ; [1]   CL (L/d)
 (0.001,4.88,10) ; [2,3] V  (L)
 ; values are determined in 3 iterations
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

$OMEGA  0.584  ;     IIV CL
 0.193  ;      IIV V
$SIGMA  0.295  ; EPS(1), proportional
 5E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
 ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            SIG=3 SIGL=3 PRINT=5 MAXEVALS=220 CTYPE=4

