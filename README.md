# case_control_logistic
This repository houses the basic R functions and examples to run the case control logistic regression model introduced in: Smith Jeffrey A., McPherson, Miller, and Lynn Smith-Lovin.  2014.“Social Distance in the United States: Sex, Race, Religion, Age and Education Homophily among Confidants, 1985-2004.” American Sociological Review 79:432-456. The model takes basic ego network data, including number of partners, ego attributes and alter attributes and estimates the strength of homophily (relative to chance) along the demographic dimensions of interest.

The R tutorial can be found in: example_estimate_homophily_ego_net_data_case_control_model.txt. The R code depends on the following packages: 
biglm, ergm and doParallel. 
