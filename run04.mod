;; 1. Based on: run03
;; 2. Description: PMX001 2CMT MM M7+
;; Joshua I., Ali K.
;; 2025-12-22

$PROBLEM PMX001 Two compartment MM M7+

$INPUT DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB SEX AGE CREAT COL EVID BLQ CONC_M7=DV

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
TVVMAX = THETA(1)
  VMAX = TVVMAX * EXP(ETA(1))
  TVKM = THETA(2)
    KM = TVKM
  TVV1 = THETA(3)
    V1 = TVV1 * EXP(ETA(2))
  TVV2 = THETA(4)
    V2 = TVV2 * EXP(ETA(3))
   TVQ = THETA(5)
     Q = TVQ
   K12 = Q/V1
   K21 = Q/V2
    S1 = V1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
    S2 = V2
    
$THETA
(0, 0.75)     ; VMAX [2]
(1E-9, 0.04) ; KM [2]
(0, 3.23)     ; V1 [5]
(0, 3.68)     ; V2 [5]
(0, 0.02)    ; Q  [5]

$DES DADT(1) = -VMAX*A(1)/(KM*V1 +A(1)) -K12*A(1) +K21*A(2) ; ODE for central    compartment
     DADT(2) =                           K12*A(1) -K21*A(2) ; ODE for peripheral compartment

$ERROR
IPRED = F
LLOQ = 1      ; mg/L

Y = IPRED + IPRED*ERR(1) + ERR(2)
W = SQRT(IPRED**2*SIGMA(1,1) + SIGMA(2,2))
IRES  = DV - IPRED
IWRES = IRES / W

$OMEGA
0.05     ; IIV VMAX
;0.05     ; IIV KM, high shrinkage
0.05     ; IIV V1 [5]
0.05     ; IIV V2 [5]
;0.0005   ; IIV Q 

$SIGMA ; residual variability
0.1    ; ERR(1); EPS(1), proportional
0.05   ; ERR(2); additive, required by M7+ censoring method

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME DV EVID PRED IPRED WRES IWRES RES IRES CWRES NOPRINT ONEHEADER FILE=run04_sdtab

$TABLE ; output table for PK parameters
ID VMAX KM V1 V2 Q K12 K21 NOPRINT NOAPPEND ONEHEADER FILE=run04_patab

$TABLE ; output table for categorical covariates
ID ADA SEX COL NOPRINT NOAPPEND ONEHEADER FILE=run04_catab

$TABLE ; output table for continuous covariates
ID BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=run04_cotab

;; REFERENCES
;; 1] https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2] https://www.tandfonline.com/doi/pdf/10.1080/19420862.2025.2512217?needAccess=true
;; 3] https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 4] https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212
;; 5] https://accp1.onlinelibrary.wiley.com/doi/10.1177/0091270006298188

