;; 1. Based on: run09
;; 2. Description: PMX001 2CMT LINEAR M7+ (COV MODEL)
;; Joshua I., Ali K.
;; 2025-12-09

$PROBLEM PMX001 Two Compartment

$INPUT DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB SEX AGE CREAT COL EVID BLQ CONC_M7=DV

$DATA dataset_clean.csv IGNORE=@

$ABBR DERIV2=NO

$SUBROUTINES
ADVAN13 ; general nonlinear model
TOL=5   ; tolerance for $DES, higher TOL: more accurate, but slower computation

$MODEL
NCOMP=2                      ; number of compartments
COMP=(DOSE, DEFDOSE, DEFOBS) ; first (central) compartment
COMP=(PERIPH)                ; second (peripheral) compartment

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

$THETA ; values are determined in 3 iterations
(0.001, 1.41,    3) ; [1]   CL (L/d)
(0.001, 4.73,    7) ; [2,3] V1  (L)
(0.01,  0.231, 0.8) ; V2
(0.001, 2.15,    5) ; Q
(0.001, 5.78,    8) ; CLADA
(-0.140,-0.127,-0.001) ; V1ALB1

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
0.608 ; IIV CL
0.375 ; IIV V1

$SIGMA    ; residual variability
0.289     ; EPS(1), proportional
1E-13 FIX ; EPS(2); additive, required by M7+ censoring method

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME TAD DV EVID PRED IPRED WRES IWRES RES IRES CWRES NOPRINT ONEHEADER FILE=run12_sdtab

$TABLE ; output table for PK parameters
ID CL V1 V2 Q K10 K12 K21 NOPRINT NOAPPEND ONEHEADER FILE=run12_patab

$TABLE ; output table for categorical covariates
ID ADA SEX COL NOPRINT NOAPPEND ONEHEADER FILE=run12_catab

$TABLE ; output table for continuous covariates
ID BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=run12_cotab

$TABLE ; output table for parameters and covariates
ID CL V1 V2 Q K10 K12 K21 ADA SEX COL BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=run12_pa_cov

;; REFERENCES
;; 1) https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2) https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3) https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

