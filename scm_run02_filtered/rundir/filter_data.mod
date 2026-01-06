;; 1. Based on: run01
;; 2. Description: PMX001 2CMT LINEAR M7+ (BASE MODEL)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA      ../../dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO COMRES=6
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=4 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=7 -> DES will aim for accuracy of 1E-7 for each integration step
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

TVCL = THETA(1)
TVV1 = THETA(2)
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

$THETA  (0.001,1.48291,8) ; [1]   CL (L/d)
 (0.001,4.61209,10) ; [2,3] V1  (L)
 (0.01,0.349516,6) ; V2
 (0.001,2.71608,6) ; Q
 ; values are determined in 3 iterations
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

"LAST
"  COM(1)=ABS(G(1,1))
"  COM(2)=ALB*ABS(G(1,1))
"  COM(3)=ABS(G(2,1))
"  COM(4)=BW*ABS(G(2,1))
"  COM(5)=ALB*ABS(G(2,1))
"  COM(6)=CREAT*ABS(G(2,1))
$OMEGA  0.573074  ;     IIV CL
 0.2141  ;      IIV V
$SIGMA  0.290424  ; EPS(1), proportional
 5E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION SIGL=3 PRINT=5 CTYPE=4 MAXEVALS=0
$TABLE      MDV ID COM(1)=CLRATIO COM(2)=CLALBNUM DV CIPREDI CIWRESI
            NOAPPEND NOPRINT ONEHEADER FORMAT=s1PE23.16
            FILE=CL_time_varying.dta
$TABLE      MDV ID COM(3)=V1RATIO COM(4)=V1BWNUM COM(5)=V1ALBNUM
            COM(6)=V1CREATNUM DV CIPREDI CIWRESI NOAPPEND NOPRINT
            ONEHEADER FORMAT=s1PE23.16 FILE=V1_time_varying.dta

