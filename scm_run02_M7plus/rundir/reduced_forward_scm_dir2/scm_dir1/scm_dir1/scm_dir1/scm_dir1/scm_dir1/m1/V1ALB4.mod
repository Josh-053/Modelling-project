;; 1. Based on: run01
;; 2. Description: PMX001 2CMT LINEAR M7+ (BASE MODEL)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX=DROP AGE=DROP CREAT=DROP COL EVID BLQ CONC_M7=DV
$DATA      ../../../../../../../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=7 -> DES will aim for accuracy of 1E-7 for each integration step
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
;;; V1ALB-DEFINITION START
   V1ALB = EXP(THETA(9)*(ALB - 31.85))
;;; V1ALB-DEFINITION END


;;; V1BW-DEFINITION START
IF(BW.LE.65.93) V1BW = ( 1 + THETA(7)*(BW - 65.93))
IF(BW.GT.65.93) V1BW = ( 1 + THETA(8)*(BW - 65.93))
;;; V1BW-DEFINITION END


;;; V1ADA-DEFINITION START
IF(ADA.EQ.0) V1ADA = 1  ; Most common
IF(ADA.EQ.1) V1ADA = ( 1 + THETA(6))
;;; V1ADA-DEFINITION END

;;; V1-RELATION START
V1COV=V1ADA*V1BW*V1ALB
;;; V1-RELATION END


;;; CLADA-DEFINITION START
IF(ADA.EQ.0) CLADA = 1  ; Most common
IF(ADA.EQ.1) CLADA = ( 1 + THETA(5))
;;; CLADA-DEFINITION END

;;; CL-RELATION START
CLCOV=CLADA
;;; CL-RELATION END


IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

TVCL = THETA(1)

TVCL = CLCOV*TVCL
TVV1 = THETA(2)

TVV1 = V1COV*TVV1
TVV2 = THETA(3)
 TVQ = THETA(4)
   
  CL = TVCL * EXP(ETA(1))
  V1 = TVV1 * EXP(ETA(2))
  V2 = TVV2
   Q = TVQ
 
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
 
  S1 = V1/1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2/1

$THETA  (0.001,1.41059,8) ; [1]   CL (L/d)
 (0.001,7.39099,10) ; [2,3] V1  (L)
 (0.01,0.182501,6) ; V2
 (0.001,1.328,6) ; Q
 ; values are determined in 3 iterations
$THETA  (-1,4.19345,5) ; CLADA1
$THETA  (-1,-0.270135,5) ; V1ADA1
$THETA  (-100000,0.0128256,0.042) ; V1BW1
 (-0.055,-0.0509691,100000) ; V1BW2
$THETA  (-0.388950184627373,0.0127266,0.388950184627373) ; V1ALB1
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

$OMEGA  0.605041  ;     IIV CL
 0.378531  ;      IIV V
$SIGMA  0.279705  ; EPS(1), proportional
 5E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            SIG=3 SIGL=4 PRINT=5 MAXEVALS=9999 CTYPE=4

