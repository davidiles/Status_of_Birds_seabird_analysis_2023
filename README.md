# Status_of_Seabirds

 Analysis of seabird colony trends for ECCC's status and trend assessments, used for 2022 State of Canada's Birds report.
 
 Analysis was conducted separately for each region (Arctic, Atlantic, Pacific).
 
 Analysis fits Generalized Additive Mixed Models (GAMMS) with Bayesian methods in JAGS.  Shape of colony-level trajectories are potentially shared among colonies, but intercepts are fit separately to each colony.  Regional trajectories (and trends) are calculated by summing colony-level trajectories.  This requires careful selection of which colonies and time periods to sum - surveys are very infrequent, and uncertainty can be substantial outside the range of survey data.
 - The choice of colonies to include in the regional summation, and the time periods to include, is currently performed somewhat subjectively based on colony size (attempting to include the largest colonies) and survey frequency (colonies are preferentially included if they are frequently surveyed).  In some cases, large colonies must be omitted because they were not surveyed frequently enough; their estimates would be dangerous extrapolations if included in a regional summary.  This component of the analysis is a priority for future improvement.
 
 For Arctic and Pacific data, some colonies are surveyed repeatedly through time with "fixed" plots on each island that provide information about trend, but are not on the same scale as "census" data where an estimate of island-level abundance is provided.  In those cases, the trajectory is subsequently rescaled based on the most recent estimate of island-level abundance, so that trajectory can be weighted correctly at the regional level.
