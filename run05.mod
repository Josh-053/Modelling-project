;; 1. Based on: run01
;; 2. Description: PMX001 1CMT LIN+MM M7+
;; Joshua I., Ali K.
;; 2025-12-29

$PROBLEM PMX001 One Compartment IV infusion with linear + Michaelis–Menten eliminations

$INPUT DUMMY=DROP ID TIME WEEK DOSE=AMT RATE CONC=DROP ADA BW ALB SEX AGE CREAT COL EVID BLQ CONC_M7=DV

$ABBR DERIV2=NO

$DATA dataset_clean.csv IGNORE=@

$SUBROUTINES
ADVAN13 ; general nonlinear model
TOL=9   ; tolerance for $DES, higher TOL: more accurate, but slower computation
        ; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step

$MODEL
NCOMP = 1                       ; number of compartments
COMP  = (DOSE, DEFDOSE, DEFOBS) ; first (central) compartment

$PK
IF (EVID.EQ.1) LASTDOSE = TIME
TAD = TIME - LASTDOSE

; Typical (population) parameters
TVCL   = THETA(1)   ; Linear clearance [L/d]
TVV    = THETA(2)   ; Central volume [L]
TVVMAX = THETA(3)   ; Vmax for MM elimination [amount/d]
TVKM   = THETA(4)   ; Km for MM elimination [conc units, e.g., mg/L]

; Individual parameters
CL   = TVCL  * EXP(ETA(1))   
V    = TVV   * EXP(ETA(2)) 
VMAX = TVVMAX
KM   = TVKM

K  = CL/V

; Scaling
S1 = V

$THETA; values were iteratively determined
(0, 1.19, 20)   ; [1] TVCL (L/d)
(0, 4.45, 10)  ; [2] TVV  (L)
(0.00001, 0.977,20) ; [3] VMAX (mg/d) 
(0.00001, 3.65, 10)  ; [3] KM   (mg/L)

$DES ; ODE for 1CMT with combined elimination
DADT(1) = - (VMAX*A(1)) / (KM*V + A(1)) - K*A(1)
 
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

$OMEGA; values were iteratively determined
0.615 ; IIV CL
0.217 ; IIV V

$SIGMA    ; residual variability
0.123     ; proportional error
5E-13 FIX ; additive error

$EST
METHOD=1 INTERACTION; FOCE-I
MAXEVAL=9999
SIG=3
SIGL=3
PRINT=5

$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S

$TABLE ; output table for standard outcomes
ID TIME DV EVID PRED IPRED WRES RES CWRES NOPRINT ONEHEADER FILE=run05_sdtab

$TABLE ; output table for PK parameters
ID CL V K NOPRINT NOAPPEND ONEHEADER FILE=run05_patab

$TABLE ; output table for categorical covariates
ID ADA SEX COL NOPRINT NOAPPEND ONEHEADER FILE=run05_catab

$TABLE ; output table for continuous covariates
ID BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=run05_cotab

$TABLE ; output table for parameters and covariates
ID CL V K ADA SEX COL BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER FILE=run05_pa_cov

;; REFERENCES
;; 1) https://www.ncbi.nlm.nih.gov/books/NBK557889/
;; 2) https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/
;; 3) https://www.tandfonline.com/doi/pdf/10.1080/19420862.2025.2512217?needAccess=true

