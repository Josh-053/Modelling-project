;; 1. Based on: run06
;; 2. Description: PMX001 2CMT LINEAR M3 (COV MODEL)
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
(0.01,  0.617,  10)   ; [1]   TVCL (L/d)
(0.01,  4.63,   15)   ; [2,3] TVV1 (L)
(0.01,  29.8,  60)   ; V2
(0.001, 0.819, 3) ; Q
(-5,    -0.283, 5)   ; CLADA

$DES DADT(1) = -K10*A(1) -K12*A(1) +K21*A(2) ; ODE for central    compartment
     DADT(2) =            K12*A(1) -K21*A(2) ; ODE for peripheral compartment
     
$ERROR
IPRED = F
LLOQ = 1  ; (mg/L)

;; residual error standard deviations
PROP_SD = SIGMA(1)
ADD_SD  = SIGMA(2)

;; final residual error
W = SQRT(ADD_SD**2+PROP_SD**2)
IF (W.LE.0.000001) W=0.000001 ; Protective code

IRES=CONC-IPRED
IWRES=IRES/W
Y = IPRED * (1 + EPS(1)) + EPS(2)

;; M3 censoring
IF (BLQ.EQ.1) THEN
  Y    = LLOQ
  CENS = 1
ELSE
  CENS = 0
  Y    = IPRED + W*EPS(3)    ; EPS(3) ~ N(0,1) (see $SIGMA third entry FIXED at 1)
ENDIF

$OMEGA ; interindividual variability
1.14  ; IIV CL
0.32  ; IIV V

$SIGMA      ; residual variability
0.0851       ; EPS(1), proportional
0.000000005 FIX ; EPS(2); additive
1 FIX       ; EPS(3); M3: normal distribution

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table
ID TIME TAD DV EVID PRED IPRED WRES IWRES RES IRES CWRES CL V1 V2 Q K10 K12 K21 ADA SEX COL BW ALB AGE CREAT NOPRINT ONEHEADER FILE=run07

;; REFERENCES
;; 1) https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2) https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3) https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

