;; 1. Based on: run02
;; 2. Description: PMX001 2CMT LINEAR M7+ (BASE MODEL)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA
            BW=DROP ALB=DROP SEX=DROP AGE=DROP CREAT COL EVID BLQ
            CONC_M7=DV
$DATA      ../../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=4 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=7 -> DES will aim for accuracy of 1E-7 for each integration step
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
;;; V1CREAT-DEFINITION START
V1CREAT = ( 1 + THETA(6)*(CREAT - 1.05))
;;; V1CREAT-DEFINITION END


;;; V1ADA-DEFINITION START
IF(ADA.EQ.0) V1ADA = 1  ; Most common
IF(ADA.EQ.1) V1ADA = ( 1 + THETA(5))
;;; V1ADA-DEFINITION END

;;; V1-RELATION START
V1COV=V1ADA*V1CREAT
;;; V1-RELATION END


IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

TVCL = THETA(1)
TVV1 = THETA(2)

TVV1 = V1COV*TVV1
TVV2 = THETA(3)
 TVQ = THETA(4)
   
  CL = TVCL
  V1 = TVV1 * EXP(ETA(1))
  V2 = TVV2
   Q = TVQ
 
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
 
  S1 = V1/1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2/1

$THETA  (0.001,0.979091,10) ; [1]   CL (L/d)
 (0.001,2.75345,10) ; [2,3] V1  (L)
 (0.01,0.0400275,10) ; V2
 (0.001,1.34461,10) ; Q
 ; values are determined in 3 iterations
$THETA  (-1,-0.879661,5) ; V1ADA1
$THETA  (-1.282,0.001163,1.163) ; V1CREAT1
$DES DADT(1) = -K10*A(1) -K12*A(1) +K21*A(2) ; ODE for central    compartment
     DADT(2) =            K12*A(1) -K21*A(2) ; ODE for peripheral compartment
     
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

$OMEGA  0.549216  ;     IIV V1
;0.573 ; IIV CL
$SIGMA  1.44474  ; EPS(1), proportional
 5E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            SIG=3 SIGL=3 PRINT=5 MAXEVALS=359 CTYPE=4

