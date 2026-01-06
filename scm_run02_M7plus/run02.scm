model=run02.mod
directory=scm_run02
search_direction=both
p_forward=0.05
p_backward=0.01
continuous_covariates=BW,ALB,AGE,CREAT
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
CL=BW,ALB,AGE,CREAT,ADA,SEX
V1=BW,ALB,AGE,CREAT,ADA,SEX

[valid_states]
continuous=1,2,3,4,5
categorical=1,2
