;; 1. Based on: final_backward_3
;; 2. Description: PMX001 1CMT LINEAR M7+
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 One Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA      dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step
$MODEL      NCOMP=1 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
$PK
;; time after dose
IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE
;; covariates on CL
IF(ADA.EQ.0) CLADA = 1  ; Most common
IF(ADA.EQ.1) CLADA = ( 1 + THETA(3))

CLCOV=CLADA

;; covariates on V
IF(ADA.EQ.0) VADA = 1  ; Most common
IF(ADA.EQ.1) VADA = ( 1 + THETA(4))
VALB   = ( 1 + THETA(5)*(ALB - 31.85))
VCREAT = ( 1 + THETA(6)*(CREAT - 1.03))

VCOV=VADA*VALB*VCREAT

TVCL = THETA(1)
  CL = TVCL * EXP(ETA(1)) * CLCOV
 TVV = THETA(2)
   V = TVV  * EXP(ETA(2)) * VCOV
 K10 = CL/V
  S1 = V ; scale prediction based on DOSE (mmol) and DV (mmol/L)

$THETA   ; values are determined in 3 iterations
(0.001,   1.4,        8) ; CL (L/d) [1]
(0.001,   4.96,      15) ; V (L) [2,3]
(0.001,   4.02,       9) ; CLADA
(-1,     -0.263,      2) ; VADA
(-0.138, -0.125,  0.084) ; VALB
(-0.4,    0.184,    0.4) ; VCREAT

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

$OMEGA
0.604 ; IIV CL (value obtained from SCM+)
0.37 ; IIV V

$SIGMA    ; residual variability
0.283  ; EPS(1), proportional (value obtained from SCM+)
5E-13 FIX ; EPS(2)  ; additive, required by M7+ censoring method

$ESTIMATION
METHOD=1 INTERACTION ; FOCE-I
SIG=3
PRINT=5
MAXEVALS=9999

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

