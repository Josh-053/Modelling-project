;; 1. Based on: run01
;; 2. Description: PMX001 2CMT LINEAR M7+ (BASE MODEL)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA
            BW=DROP ALB SEX=DROP AGE=DROP CREAT=DROP COL EVID BLQ
            CONC_M7=DV
$DATA      ../../../../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=7 -> DES will aim for accuracy of 1E-7 for each integration step
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
;;; CLALB-DEFINITION START
CLALB = ( 1 + THETA(8)*(ALB - 31.88))
;;; CLALB-DEFINITION END


;;; V1ALB-DEFINITION START
IF(ALB.LE.31.85) V1ALB = ( 1 + THETA(6)*(ALB - 31.85))
IF(ALB.GT.31.85) V1ALB = ( 1 + THETA(7)*(ALB - 31.85))
;;; V1ALB-DEFINITION END

;;; V1-RELATION START
V1COV=V1ALB
;;; V1-RELATION END


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

$THETA  (0.001,1.42126,8) ; [1]   CL (L/d)
 (0.001,5.53619,10) ; [2,3] V1  (L)
 (0.01,0.205842,6) ; V2
 (0.001,1.47189,6) ; Q
 ; values are determined in 3 iterations
$THETA  (-1,4.67297,5) ; CLADA1
$THETA  (-100000,0.00668625,0.084) ; V1ALB1
 (-0.134,-0.126597,100000) ; V1ALB2
$THETA  (-0.135,8.4E-05,0.084) ; CLALB1
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

$OMEGA  0.592834  ;     IIV CL
 0.28635  ;      IIV V
$SIGMA  0.280484  ; EPS(1), proportional
 5E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            SIG=3 SIGL=4 PRINT=5 MAXEVALS=9999 CTYPE=4

