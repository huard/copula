% TEST SUITE FOR COPULA RELATED FUNCTIONS
% Some tests check that the functions return matrices of the correct shape.
% Other functions check that the actual values are vaild, by comparing
% results with those obtained using Maple. See file analytic.mws to see how
% these computations were performed.

% Modified by:
%% D. Huard,  Nov. 8, 2006. Added TEST CHECK_TAU, TEST CHECK_ALPHA
%% D. Huard, Nov, 9, 2006. Added TEST COPULAPDF, TEST COPULACDF, COPULA
%%  TEST COPULAPARAM, TEST COPULASTAT
%%  Found bug in copulacdf('clayton')
%% D. Huard, Nov, 12, 2006. Added TAUJACOBIAN
%%  Found bug in copulacdf('fgm')

warning('off');%, COPULA:BadParameter)

% TEST COPULA_LIKE
run test_posterior
run test_check_tau
run test_check_alpha
run test_copulastat
run test_copulaparam
run test_copulapdf
run test_copulacdf
run test_taujacobian
run test_copularnd
%run test_copularnd2

fprintf('Everything looks fine.\n')

% TEST BCS
fprintf('\nRunning the example...\n')
run example

