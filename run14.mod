;; 1. Based on: run11
;; 2. Description: PMX001 2CMT TMDD M7+
;; Joshua I., Ali K.
;; 2025-12-09

$PROBLEM PMX001 Two Compartment TMDD Model
;; Two assumptions:
;; 1. Free ligand concentration is much larger than total target concentration.
;; 2. Vmax = cst => Rtot = cst => target degradation rate (kout) = ligand-target complex internalization rate (kint)


$INPUT DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB SEX AGE CREAT COL EVID BLQ CONC_M7=DV

$DATA dataset_clean.csv IGNORE=@

$ABBR DERIV2=NO

$SUBROUTINES
ADVAN13 ; general nonlinear model
TOL=5   ; tolerance for $DES, higher TOL: more accurate, but slower computation
        ; TOL=5 -> accurate to 1E-5, cannot increase accuracy without risking run failures & too high RSEs

$MODEL
NCOMP=2                      ; number of compartments
COMP=(DOSE, DEFDOSE, DEFOBS) ; first (central) compartment
COMP=(PERIPH)                ; second (peripheral) compartment

$PK
IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

  TVCL = THETA(1)
TVVMAX = THETA(2)
  TVKM = THETA(3)
  TVV1 = THETA(4)
  TVV2 = THETA(5)
   TVQ = THETA(6)

;IF(ADA.EQ.0) CLADA = 1  ; Most common
;IF(ADA.EQ.1) CLADA = ( 1 + THETA(7))
   
  CL = TVCL   * EXP(ETA(1))
VMAX = TVVMAX
  KM = TVKM
  V1 = TVV1   * EXP(ETA(2))
  V2 = TVV2
   Q = TVQ
 
 K10 = CL/V1
 K12 = Q/V1
 K21 = Q/V2
 
  S1 = V1/1 ; scale prediction based on DOSE (mmol) and DV (mmol/L)
  S2 = V2/1

$THETA ; values are determined in 3 iterations
(0.01, 0.524, 1) ; [1] CL (L/d)
(0.5, 2.56, 6)  ; VMAX
(0.0001, 0.0011, 0.1) ; KM
(0.5, 6.2, 12)   ; [2,3] V1 (L)
(0.01,  0.0463, 0.2)  ; V2
(0.1, 1.1, 2)   ; Q

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

$OMEGA
0.626 ; IIV CL
0.431 ; IIV V1

$SIGMA    ; residual variability
0.12     ; EPS(1), proportional
1E-13 FIX ; EPS(2); additive, required by M7+ censoring method

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
SIGL=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S


;; REFERENCES
;; 1) https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2) https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3) https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

