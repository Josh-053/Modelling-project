model=run02.mod
directory=scm_run02
search_direction=both
p_forward=0.05
p_backward=0.01
continuous_covariates=BW,ALB,AGE,CREAT
categorical_covariates=ADA
time_varying=BW,ALB,CREAT

missing_data_token=-99
linearize=0
nm_version=default
retries=1
threads=6
tweak_inits=1
picky=0

[test_relations]
CL=ALB,AGE,ADA
V1=BW,ALB,CREAT,ADA

[valid_states]
continuous=1,2,3,5
categorical=1,2

[included_relations]
CL=ADA-2,AGE-2,ALB-3
V1=ADA-2,ALB-3,BW-3,CREAT-2
