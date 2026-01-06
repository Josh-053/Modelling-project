;; 1. Based on: final_backward
;; 2. Description: PMX001 2CMT LINEAR M7+ (COV MODEL)
;; Joshua I., Ali K.
;; 2026-01-05
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA      dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=7 -> DES will aim for accuracy of 1E-7 for each integration step
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
;;; V1BW-DEFINITION START
IF(BW.LE.65.93) V1BW = ( 1 + THETA(9)*(BW - 65.93))
IF(BW.GT.65.93) V1BW = ( 1 + THETA(10)*(BW - 65.93))
;;; V1BW-DEFINITION END


;;; V1ALB-DEFINITION START
V1ALB = ( 1 + THETA(8)*(ALB - 31.85))
;;; V1ALB-DEFINITION END

;;; V1-RELATION START
V1COV=V1ALB*V1BW
;;; V1-RELATION END


;;; CLALB-DEFINITION START
IF(ALB.LE.31.88) CLALB = ( 1 + THETA(6)*(ALB - 31.88))
IF(ALB.GT.31.88) CLALB = ( 1 + THETA(7)*(ALB - 31.88))
;;; CLALB-DEFINITION END


;;; CLADA-DEFINITION START
IF(ADA.EQ.0) CLADA = 1  ; Most common
IF(ADA.EQ.1) CLADA = ( 1 + THETA(5))
;;; CLADA-DEFINITION END

;;; CL-RELATION START
CLCOV=CLADA*CLALB
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

$THETA ; values are determined in 3 iterations
(0.001,   0.997,    4) ; [1]   CL (L/d)
(0.001,   6.66,    10) ; [2,3] V1  (L)
(0.01,    0.251,    2) ; V2
(0.001,   1.93,     3) ; Q
(-1,      4.61,    10) ; CLADA
(-2,     -0.127,  0.4) ; CLALB1
(-0.135,  0.313,   10) ; CLALB2
(-0.134, -0.123,  0.4) ; V1ALB
(-2,      0.0206, 0.4) ; V1BW1
(-0.055, -0.0544, 0.4) ; V1BW2
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

$OMEGA
0.694  ; IIV CL
0.522  ; IIV V
$SIGMA    ; residual variability
0.285     ; EPS(1); proportional
5E-13 FIX ; EPS(2); additive, required by M7+ censoring method

    
$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME TAD DV EVID PRED IPRED WRES IWRES RES IRES CWRES NOPRINT ONEHEADER FILE=backward2_sdtab

$TABLE ; output table for PK parameters
ID CL V1 V2 Q K10 K12 K21 NOPRINT NOAPPEND ONEHEADER FILE=backward2_patab

$TABLE ; output table for categorical covariates
ID ADA SEX COL NOPRINT NOAPPEND ONEHEADER FILE=backward2_catab

$TABLE ; output table for continuous covariates
ID BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=backward2_cotab

$TABLE ; output table for parameters and covariates

