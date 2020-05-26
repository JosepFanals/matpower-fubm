function [results, iBeqz, iBeqv] = test_fubm
%% TEST_FUBM Tests FUBM 
clear all;
clc; 

t_jacobian_fubm;
t_hessian_fubm;
t_pf_acdc_fubm;
