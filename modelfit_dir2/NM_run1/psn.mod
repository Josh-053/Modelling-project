;; 1. Based on: run01
;; 2. Description: PMX001 One Compartment
;; x1. Joshua I., Ali K.
;; 2025-12-09
$PROBLEM    PMX001 One Compartment
$INPUT      DUMMY=DROP ID TIME WEEK DOSE=DROP RATE CONC=DV ADA BW ALB
            SEX AGE CREAT COL EVID AMT
$DATA      dataset_clean.csv IGNORE=@
$SUBROUTINE ADVAN13 ; general nonlinear model
            TOL=9 ; tolerance for $DES,higher TOL: more accurate,but slower computation

; TOL=9 -> DES will aim for accuracy of 1E-9 for each integration step
$MODEL      NCOMP=1 ; number of compartments
            COMP=(DOSE,DEFDOSE,DEFOBS) ; first (central) compartment
$PK
TVCL = THETA(1)
  CL = TVCL * EXP(ETA(1))
 TVV = THETA(2)
   V = TVV  * EXP(ETA(2))
 K10 = CL/V
  S1 = V ; scale prediction based on DOSE (mmol) and DV (mmol/L)

$DES DADT(1) = -K10*A(1)  ; ODE for central compartment

$ERROR
IPRED = F
LLOQ = 1

PROP = IPRED*THETA(3)
ADD = THETA(4)
IF (BLQ.EQ.1) ADD = ADD + LLOQ

W = SQRT(ADD**2+PROP**2)

IF (W.LE.0.000001) W=0.000001 ; Protective code

IRES=CONC-IPRED
IWRES=IRES/W

Y = IPRED + W*ERR(1)

$THETA  (0,0.012) ; CL [1]
 (0,4.45) ; V [2,3]
 (0,0.2,0.5) ; PROP error ()
 (0,0.2,10) ; ADD (mg/L)
$OMEGA  0.05  ; IIV/BSV CL, fix to 0 to exclude
 0.05  ;  IIV/BSV V
$SIGMA  1  FIX  ; errors already adjusted under $ERROR as THETAs
 ; residual variability
$ESTIMATION METHOD=1 INTERACTION ; FOCE-I
            MAXEVAL=9999 SIG=3 PRINT=5
$COVARIANCE PRINT=E UNCONDITIONAL MATRIX=S
 ; output table for standard outcomes
$TABLE      ID TIME DV EVID PRED IPRED WRES IWRES RES IRES CWRES
            NOPRINT ONEHEADER FILE=run01_sdtab
 ; output table for PK parameters
$TABLE      ID CL V K10 NOPRINT NOAPPEND ONEHEADER FILE=run01_patab
 ; output table for categorical covariates
$TABLE      ID ADA SEX COL NOPRINT NOAPPEND ONEHEADER FILE=run01_catab
 ; output table for continuous covariates
$TABLE      ID BW ALB AGE CREAT NOPRINT NOAPPEND ONEHEADER
            FILE=run01_cotab
;; REFERENCES

;; 1] https://www.ncbi.nlm.nih.gov/books/NBK557889/

;; 2] https://pmc.ncbi.nlm.nih.gov/articles/PMC8301575/

;; 3] https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1016/j.clpt.2004.12.212

