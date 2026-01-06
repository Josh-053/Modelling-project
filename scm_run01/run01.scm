model=run01.mod
directory=scm_run01_cont
search_direction=both
p_forward=0.05
p_backward=0.01
continuous_covariates=BW,ALB,CREAT
categorical_covariates=ADA,SEX
time_varying=BW,ALB,CREAT

missing_data_token=-99
linearize=0
nm_version=default
retries=1
threads=6
tweak_inits=1
picky=0

[test_relations]
CL=BW,ALB,CREAT,ADA,SEX
V=BW,ALB,CREAT,ADA

[valid_states]
continuous=1,2,3
categorical=1,2

[included_relations]
CL=ADA-2
V=ALB-3,BW-2
