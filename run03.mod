;; 1. Based on: run01
;; 2. Description: PMX001 1CMT MM M7+
;; Joshua I., Ali K.
;; 2025-12-09

$PROBLEM PMX001 One Compartment, Michaelis-Menten

$INPUT DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB SEX AGE CREAT COL EVID BLQ CONC_M7=DV

$ABBR DERIV2=NO

$DATA dataset_clean.csv IGNORE=@

$SUBROUTINES
ADVAN13 ; general nonlinear model
TOL=9   ; tolerance for $DES, higher TOL: more accurate, but slower computation
        ; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step

$MODEL
NCOMP=1                      ; number of compartments
COMP=(DOSE, DEFDOSE, DEFOBS) ; first (central) compartment

$PK
TVVMAX = THETA(1)
  VMAX = TVVMAX
  TVKM = THETA(2)
    KM = TVKM
   TVV = THETA(3)
     V = TVV * EXP(ETA(1))
    S1 = V ; scale prediction based on DOSE (mmol) and DV (mmol/L)

$THETA
(0,    0.75)    ; VMAX [1]
(1E-9, 0.04)    ; KM [1]
(0, 5.16, 10)   ; [2,3] TVV  (L)

$DES
DADT(1) = -VMAX*A(1)/(KM*V + A(1))  ; ODE for central compartment, Michaelis-Menten elimination

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

$OMEGA ; IIV/BSV, fix to 0 to exclude
0.237 ; IIV/BSV V

$SIGMA    ; residual variability
0.167    ; ERR(1); EPS(1), proportional
0.000000005   ; ERR(2); additive, required by M7+ censoring method

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME DV EVID PRED IPRED WRES IWRES RES IRES CWRES VMAX KM V ADA SEX COL BW ALB AGE CREAT NOPRINT ONEHEADER FILE=run02

;; REFERENCES
;; 1] https://www.tandfonline.com/doi/pdf/10.1080/19420862.2025.2512217?needAccess=true
;; 2] https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3] https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

