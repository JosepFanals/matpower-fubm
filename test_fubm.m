function [results, iBeqz, iBeqv] = test_fubm
%% TEST_FUBM Tests the FUBM Derivatives against Finite Differences Method
clear all;
clc; 

t_jacobian_fubm;
t_hessian_fubm;

