;; 1. Based on: run02
;; 2. Description: PMX001 2CMT LINEAR M7+ (BASE MODEL)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX=DROP AGE=DROP CREAT=DROP COL EVID BLQ CONC_M7=DV
$DATA      ../../../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=4 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=7 -> DES will aim for accuracy of 1E-7 for each integration step
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
;;; V1BW-DEFINITION START
V1BW = ( 1 + THETA(7)*(BW - 65.91))
;;; V1BW-DEFINITION END


;;; V1ALB-DEFINITION START
V1ALB = ( 1 + THETA(6)*(ALB - 31.84))
;;; V1ALB-DEFINITION END


;;; V1ADA-DEFINITION START
IF(ADA.EQ.0) V1ADA = 1  ; Most common
IF(ADA.EQ.1) V1ADA = ( 1 + THETA(5))
;;; V1ADA-DEFINITION END

;;; V1-RELATION START
V1COV=V1ADA*V1ALB*V1BW
;;; V1-RELATION END


IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

TVCL = THETA(1)
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

$THETA  (0.001,1.53862,4.5) ; [1]   CL (L/d)
 (0.001,4.74949,8) ; [2,3] V1  (L)
 (0.01,0.107995,3) ; V2
 (0.001,1.24253,6) ; Q
 ; values are determined in 3 iterations
$THETA  (-1,-0.268132,5) ; V1ADA1
$THETA  (-0.134,-0.127597,0.085) ; V1ALB1
$THETA  (-0.054,-0.0535351,0.042) ; V1BW1
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

$OMEGA  0.608261  ;     IIV CL
 0.64993  ;     IIV V1
$SIGMA  0.28955  ; EPS(1), proportional
 1E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            SIG=3 SIGL=3 PRINT=5 MAXEVALS=9999 CTYPE=4

