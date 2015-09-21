%% ChirpQualityCriteriaScatter
% Compares quality criteria used to evaluate units for their responisveness to the Chirp-stimulus
% Does so by (a) correlating the criteria-values with each other (b) comparing criteria-value with
% PSTH-profiles (visual inspection)

%% General admin
clear all;

%% Inter-criteria-correlation
% Note to self: to access all elements of a field within a structure (and not just the first one),
% and get a vector, put the argument in square brackets (and transpose if necessary), as shown below.
% Likewise, to produce a cell array, use curly braces.

load('units_for_chirp_sorted2.mat');
% Exclude some weird data points form unsorted units_for_chirp

scatter([units_for_chirp_sorted.ranksum]', [units_for_chirp_sorted.qi]');
xlabel({'Miro ranksum criterion' ' -the smaller the better-'});
ylabel({'Philipp qi criterion' '-the larger the better?- '});
title('ChirpQualityCriterionCorrelation: ranksum vs qi');

