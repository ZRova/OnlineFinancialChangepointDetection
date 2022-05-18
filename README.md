# OnlineFinancialChangepointDetection

This project looked at an online changepoint detection method proposed in Fearnhead Liu 2007 (https://doi.org/10.1111/j.1467-9868.2007.00601.x). 

Strengths and weaknesses were discussed, and several suggestions were made to improve the reliability of the method, including better estimation of mean, and updating of distribution of changepoints using a Beta posterior distribution. 
In addition, using Median Absolute Deviation (MAD) to estimate variance the method could be made entirely non parametric. 

These suggestions were then implemented, and can be found in Code\FLFunction, where parameters can be passed to determine which estimations should be used. Calling 
"FL(data)" will run the method with no adaptations. 

The project concluded by running a simulation study comparing the original method with with the adapted method, where there was a significant reduction in false positives was found. 
