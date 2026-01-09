;; 1. Based on: run13
;; 2. Description: PMX001 2CMT TMDD M7+ (BASE)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment TMDD Model
;; 2 assumptions:

;; 1. Free ligand concentration is much larger than total target concentration.

;; 2. Vmax = cst => Rtot = cst => target degradation rate (kout) = ligand-target complex internalization rate (kint)
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA=DROP
            BW=DROP ALB SEX=DROP AGE=DROP CREAT=DROP COL EVID BLQ
            CONC_M7=DV
$DATA      ../../../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=5 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=7 -> DES will aim for accuracy of 1E-7 for each integration step
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
;;; CLALB-DEFINITION START
   CLALB = ((ALB/31.85)**THETA(7))
;;; CLALB-DEFINITION END

;;; CL-RELATION START
CLCOV=CLALB
;;; CL-RELATION END


IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

  TVCL = THETA(1)

TVCL = CLCOV*TVCL
TVVMAX = THETA(2)
  TVKM = THETA(3)
  TVV1 = THETA(4)
  TVV2 = THETA(5)
   TVQ = THETA(6)
   
  CL = TVCL   * EXP(ETA(1))
VMAX = TVVMAX
  KM = TVKM
  V1 = TVV1
  V2 = TVV2
   Q = TVQ
 
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
 
  S1 = V1/1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2/1

$THETA  (0.01,1.0365,2.3) ; [1] CL (L/d)
 (0.5,1.38119,3) ; VMAX
 (0.001,0.00787185,0.02) ; KM
 (0.5,4.1052,9) ; [2,3] V1 (L)
 (0.00001,0.00029715,0.0006) ; V2
 (0.1,13.4136,30) ; Q
 ; values are determined in 3 iterations
$THETA  (-100,-1.41014,100000) ; CLALB1
$DES DADT(1) = -K10*A(1) -K12*A(1) +K21*A(2) -VMAX*A(1)/(KM*V1 + A(1)) ; ODE central    compartment
     DADT(2) =            K12*A(1) -K21*A(2)                           ; ODE peripheral compartment
     
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

$OMEGA  0.364473  ;     IIV CL
$SIGMA  0.330411  ; EPS(1), proportional
 1E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            SIG=3 SIGL=3 PRINT=5 MAXEVALS=303 CTYPE=4

