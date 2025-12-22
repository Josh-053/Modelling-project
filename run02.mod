;; 1. Based on: run01
;; 2. Description: PMX001 2CMT LINEAR
;; x1. Joshua I., Ali K.
;; 2025-12-09

$PROBLEM PMX001 Two Compartment

$INPUT DUMMY=DROP ID TIME WEEK DOSE=DROP RATE CONC=DV ADA BW ALB SEX AGE CREAT COL EVID AMT

$DATA
dataset_clean.csv
IGNORE=@

$SUBROUTINES
ADVAN13 ; general nonlinear model
TOL=9   ; tolerance for $DES, higher TOL: more accurate, but slower computation
        ; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step

$MODEL
NCOMP=2                      ; number of compartments
COMP=(DOSE, DEFDOSE, DEFOBS) ; first (central) compartment
COMP=(PERIPH)                ; second (peripheral) compartment

$PK
TVCL = THETA(1)
  CL = TVCL * EXP(ETA(1))
TVV1 = THETA(2)
  V1 = TVV1 * EXP(ETA(2))
TVV2 = THETA(3)
  V2 = TVV2
 TVQ = THETA(4)
   Q = TVQ
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
  S1 = V1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2

$THETA
(0, 0.012) ; CL [1]
(0, 4.45)  ; V1 [2,3]
(0, 10)    ; V2
(0,1)      ; Q

$OMEGA
0.05 ; IIV CL
0.05 ; IIV V1

$SIGMA ; residual variability
0.1    ; proportional error

$DES DADT(1) = -K10*A(1) -K12*A(1) +K21*A(2) ; ODE for central    compartment
     DADT(2) =            K12*A(1) -K21*A(2) ; ODE for peripheral compartment
     
$ERROR
IPRED = F
Y = IPRED + IPRED*ERR(1)
W = SQRT(IPRED**2*SIGMA(1,1)) ;; new
IRES  = CONC - IPRED
IWRES = IRES / W

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME DV EVID PRED IPRED WRES IWRES RES IRES CWRES NOPRINT ONEHEADER FILE=run02_sdtab

$TABLE ; output table for PK parameters
ID CL V1 V2 Q K10 K12 K21 NOPRINT NOAPPEND ONEHEADER FILE=run02_patab

$TABLE ; output table for categorical covariates
ID ADA SEX COL NOPRINT NOAPPEND ONEHEADER FILE=run02_catab

$TABLE ; output table for continuous covariates
ID BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=run02_cotab

;; REFERENCES
;; 1) https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2) https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3) https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

