$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB
            SEX AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA      bs_pr1_607.dta IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=5 ; tolerance for $DES,higher TOL: more accurate,but slower computation
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

IF(ADA.EQ.0) CLADA = 1  ; Most common
IF(ADA.EQ.1) CLADA = ( 1 + THETA(5))
V1ALB = ( 1 + THETA(6)*(ALB - 31.84))
   
  CL = TVCL * EXP(ETA(1)) * CLADA
  V1 = TVV1 * EXP(ETA(2)) * V1ALB
  V2 = TVV2
   Q = TVQ
 
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
 
  S1 = V1/1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2/1

$THETA  (0.001,1.39726,3) ; [1]   CL (L/d)
 (0.001,4.70664,7) ; [2,3] V1  (L)
 (0.01,0.228756,0.8) ; V2
 (0.001,2.1274,5) ; Q
 (0.001,5.78029,8) ; CLADA
 (-0.140,-0.127136,-0.001) ; V1ALB1
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

$OMEGA  0.607047  ;     IIV CL
 0.375377  ;     IIV V1
$SIGMA  0.289024  ; EPS(1), proportional
 1E-13  FIX  ;     EPS(2)  ; additive, required by M7+ censoring method
    ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            MAXEVAL=9999 SIG=3 PRINT=5

