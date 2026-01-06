;; 1. Based on: run02
;; 2. Description: PMX001 2CMT LINEAR M7+ (BASE MODEL)
;; Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 Two Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW
            ALB SEX AGE CREAT COL EVID BLQ CONC_M7=DV
$DATA       dataset_clean.csv IGNORE=@
$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate, but slower computation
$MODEL      NCOMP=2 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
            COMP=(PERIPH) ; second (peripheral) compartment
$PK
;; time after dose
IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

;; typical values
TVCL = THETA(1)
TVV1 = THETA(2)
TVV2 = THETA(3)
 TVQ = THETA(4)

;; covariates on V1
IF(ADA.EQ.0)    V1ADA =  1                          ; ADA=0: most common
IF(ADA.EQ.1)    V1ADA = (1 + THETA(5))              ; ADA=1
IF(BW.LE.65.97) V1BW  = (1 + THETA(6)*(BW - 65.97)) ; BW below median
IF(BW.GT.65.97) V1BW  = (1 + THETA(7)*(BW - 65.97)) ; BW above median

;; individual values
  CL = TVCL
  V1 = TVV1 * EXP(ETA(1)) * V1ADA * V1BW
  V2 = TVV2
   Q = TVQ
 
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
 
  S1 = V1/1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2/1

$THETA ; values are determined in 3 iterations
(0.001,0.919,5) ; [1] TVCL (L/d)
(0.001,3.1,  8) ; [2,3] TVV1 (L)
(0.01, 0.154,3) ; TVV2
(0.001,20.4, 30) ; TVQ
(-1,-0.888,5) ; V1ADA
(-100,0.0285,0.042) ; V1BW1
(-0.055,-0.0544,100) ; V1BW2
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

$OMEGA  0.768  ;     IIV V1
$SIGMA    ; residual variability
1.27   ; EPS(1) ; proportional
5E-13 FIX ; EPS(2) ; additive, required by M7+ censoring method

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
SIGL=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME TAD DV EVID PRED IPRED WRES IWRES RES IRES CWRES NOPRINT ONEHEADER FILE=finalbw1_sdtab

$TABLE ; output table for PK parameters
ID CL V1 V2 Q K10 K12 K21 NOPRINT NOAPPEND ONEHEADER FILE=finalbw1_patab

$TABLE ; output table for categorical covariates
ID ADA SEX COL NOPRINT NOAPPEND ONEHEADER FILE=finalbw1_catab

$TABLE ; output table for continuous covariates
ID BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=finalbw1_cotab

$TABLE ; output table for parameters and covariates
ID CL V1 V2 Q K10 K12 K21 ADA SEX COL BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=finalbw1_pa_cov

;; REFERENCES
;; 1) https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2) https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3) https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

