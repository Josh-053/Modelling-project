;; 1. Based on: run02
;; 2. Description: PMX001 2CMT LINEAR M7+ (COV MODEL)
;; Joshua I., Ali K.
;; 2025-12-29

$PROBLEM PMX001 Two Compartment

$INPUT DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB SEX AGE CREAT COL EVID BLQ CONC_M7=DV

$DATA dataset_clean.csv IGNORE=@

$ABBR DERIV2=NO

$SUBROUTINES
ADVAN13 ; general nonlinear model
TOL=9   ; tolerance for $DES, higher TOL: more accurate, but slower computation
        ; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step

$MODEL
NCOMP=2                      ; number of compartments
COMP=(DOSE, DEFDOSE, DEFOBS) ; first (central) compartment
COMP=(PERIPH)                ; second (peripheral) compartment

$PK
;; time after last dose
IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

;; typical values
 TVCL = THETA(1)
 TVV1 = THETA(2)
 TVV2 = THETA(3)
  TVQ = THETA(4)

;; covariate determined from SCM+ based on base model
IF(ADA.EQ.0) CLADA = 1               ; Most common
IF(ADA.EQ.1) CLADA = ( 1 + THETA(5))
CLCOV = CLADA  

;; individual values
   CL = TVCL * EXP(ETA(1)) * CLCOV
   V1 = TVV1 * EXP(ETA(2))
   V2 = TVV2
    Q = TVQ

;; derived parameters
  K10 = CL/V1
  K12 = Q/V1
  K21 = Q/V2

;; scaling factors
   S1 = V1/1 ; scale prediction based on DOSE (mmol) and CONC (mmol/L)
   S2 = V2/1

$THETA ; values are determined in 3 iterations
(0.01,  1.32,   3)   ; [1]   TVCL (L/h)
(0.01,  4.68,   8)   ; [2,3] TVV1 (L)
(0.01,  0.226,  1)   ; V2
(0.001, 0.0546, 1.5) ; Q
(-1,    0.466,  5)   ; CLADA

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

$OMEGA ; interindividual variability
0.539  ; IIV CL
0.235  ; IIV V

$SIGMA          ; residual variability
0.131           ; EPS(1), proportional
0.000000005 FIX ; EPS(2); additive, required by M7+ censoring method

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table
ID TIME TAD DV EVID PRED IPRED WRES IWRES RES IRES CWRES CL V1 V2 Q K10 K12 K21 ADA SEX COL BW ALB AGE CREAT NOPRINT ONEHEADER FILE=run06

;; REFERENCES
;; 1) https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2) https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3) https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

