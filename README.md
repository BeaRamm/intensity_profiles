# intensity_profiles
%%%%% function analyze_bathtubs_profiles(folder) %%%%%
% function analyzes time-averaged intensity profiles in
% microcompartments. 
% INPUT: folder is a string for the target folder in which the 
% source profiles are located. Source files are saved as csv files. These profiles are obtained 
% using FIJI
%
% Extracted data sets are stored in a data structure
%
% Note the convention for filenames: 
% *C2* corresponds to MinD-GFP, *C1* correspopnds to the  red % channel; for streptavidin it is the inverse
% For detailed information: Ramm et al, Nat. Commun. 2018
%
% Jonas Mücksch 2018
