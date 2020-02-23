# Disclaimer: 
This Repository is still very much work in progress. Therefore it is not yet thoroughly tested, contains only minimal documentation and is not optimized.
<br>
# Multiple Treatments for Uplift Modeling
This repository looks at uplift modeling specifically in cases where we want to campare multiple treatments. To do this several uplift methods are implemented and/or adapted for the multiple treatment case.
<br>
# Implemented Models<br>
## Tree after Rzepakowski & Jorszewicz
The decision tree is based on the one proposed by Piotr Rzepakowski & Szymon Jaroszewicz ["Decision trees for uplift modeling with single
and multiple treatments"](https://core.ac.uk/download/pdf/81899141.pdf)
<br>
## Causal Tree & Causal Forest
These models, which can be found at (https://github.com/susanathey/causalTree) are adapted for the multiple treatment case. To do that we build one model for each treatment and then compare the predicted uplifts.
<br>
## Contextual Treatment Selection (CTS)
This implementation follows the approach outlined by Yan Zhao, Xiao Fang and David Simchi-Levi in their paper "Uplift Modeling with Multiple Treatments and General Response Types" (https://arxiv.org/pdf/1705.08492.pdf). 
<br>

# Example
An example of the code applied to the Hillstrom data set (https://blog.minethatdata.com/2008/03/minethatdata-e-mail-analytics-and-data.html) can be found in "Hillstrom Evaluation Spend.R"
<br>
<br>
<br>

This repository was first created by Matthias Becher and Jan Krol as part of the Applied Predictive Analytics Seminar at Humboldt-University of Berlin. <br>
Currently it is being maintained and improved by Matthias Becher as part of his Master Thesis.