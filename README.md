# Multiple Treatments for Uplift Modeling
# Disclaimer: 
This Repository is still work in progress and not yet thoroughly tested. If you find any mistakes or have any questions, feel free to write me at becherma@hu-berlin.de
<br>
#Introduction
Uplift modeling is a predictive modeling technique which is concerned with directly estimating the effect of a treatment on a persons behaviour. Specifiacally it tries to estimate the conditional average treatment effect (CATE) of a treatment based on the characteristics of a given person (age, gender, income, ...). While the type of predictive modeling has applications in a wide variety of field, we will focus on a marketing framework. Therefore when we talk about "treatments" we usually refer to marketing campaigns, ads, etc. <br>
So far most research in this field has been done on a single treatment case and on cases with a binary outcome. We aim to fill that gap by looking exclusively at cases where multiple possible treatments are compared and focus on continuous outcome variables.
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