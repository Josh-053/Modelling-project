;; 1. Based on: run02
;; 2. Description: PMX001 2CMT LINEAR M7+ (BASE MODEL)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA      dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO COMRES=4
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
   
  CL = TVCL
  V1 = TVV1 * EXP(ETA(1))
  V2 = TVV2
   Q = TVQ
 
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
 
  S1 = V1/1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2/1

$THETA  (0.001,0.934671,10) ; [1]   CL (L/d)
 (0.001,2.1793,10) ; [2,3] V1  (L)
 (0.01,0.0395178,10) ; V2
 (0.001,1.97824,10) ; Q
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
"  COM(2)=BW*ABS(G(1,1))
"  COM(3)=ALB*ABS(G(1,1))
"  COM(4)=CREAT*ABS(G(1,1))
$OMEGA  0.438106  ;     IIV V1
;0.573 ; IIV CL
$SIGMA  1.32012  ; EPS(1), proportional
 5E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION SIGL=3 PRINT=5 CTYPE=4 MAXEVALS=0
$TABLE      MDV ID COM(1)=V1RATIO COM(2)=V1BWNUM COM(3)=V1ALBNUM
            COM(4)=V1CREATNUM DV CIPREDI CIWRESI NOAPPEND NOPRINT
            ONEHEADER FORMAT=s1PE23.16 FILE=V1_time_varying.dta

