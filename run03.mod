;; 1. Based on: run01
;; 2. Description: PMX001 1CMT MM M7+
;; x1. Joshua I., Ali K.
;; 2025-12-09

$PROBLEM PMX001 One Compartment, Michaelis-Menten

$INPUT DUMMY=DROP ID TIME WEEK DOSE=DROP RATE CONC=DV ADA BW ALB SEX AGE CREAT COL EVID BLQ AMT

$DATA
dataset_clean.csv
IGNORE=@

$SUBROUTINES
ADVAN13 ; general nonlinear model
TOL=9   ; tolerance for $DES, higher TOL: more accurate, but slower computation
        ; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step

$MODEL
NCOMP=1                      ; number of compartments
COMP=(DOSE, DEFDOSE, DEFOBS) ; first (central) compartment

$PK
TVVMAX = THETA(1)
  VMAX = TVVMAX   * EXP(ETA(1))
  TVKM = THETA(2)
    KM = TVKM     * EXP(ETA(2))
   TVV = THETA(3)
     V = TVV      * EXP(ETA(3))
    S1 = V ; scale prediction based on DOSE (mmol) and DV (mmol/L)

$DES
DADT(1) = -VMAX*A(1)/(KM*V + A(1))  ; ODE for central compartment, Michaelis-Menten elimination

$ERROR
IPRED = F
LLOQ = 1

PROP = IPRED*THETA(4)
ADD = THETA(5)
IF (BLQ.EQ.1) ADD = ADD + LLOQ

W = SQRT(ADD**2+PROP**2)

IF (W.LE.0.000001) W=0.000001 ; Protective code

IRES=CONC-IPRED
IWRES=IRES/W

Y = IPRED + W*ERR(1)

$THETA
(0, 0.012)    ; Vm
(0, 1)        ; Km
(0, 4.45)     ; V
(0, 0.2, 0.5) ; PROP error ()
(0, 0.2, 10)  ; ADD (mg/L)


$OMEGA; IIV/BSV, fix to 0 to exclude
0.05 ; VMAX
0.05 ; KM
0.05 ; V

$SIGMA ; residual variability
1 FIX ; errors already adjusted under $ERROR as THETAs

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME DV EVID PRED IPRED WRES IWRES RES IRES CWRES VMAX KM V ADA SEX COL BW ALB AGE CREAT NOPRINT ONEHEADER FILE=run02

;; REFERENCES
;; 1] https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2] https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3] https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

