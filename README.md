# Multiple Treatments for Uplift Modeling
# Disclaimer: 
This repository was first created by Matthias Becher and Jan Krol as part of the Applied Predictive Analytics Seminar at Humboldt-University of Berlin. <br>
Afterwards I continued to work on it for my master's thesis.

I currently do not actively work on this repository. However, if you find any mistakes or have any questions, feel free to write me at matthias.becher@charite.de
<br>
# Abstract
The evaluation and selection of so-called treatments (e.g. marketing actions, advertisements) is one of the key challenges for marketing practitioners today. Uplift modeling is a machine learning tool for causal modeling and commonly used for this task. This thesis looks at uplift modeling in cases where multiple treatments are evaluated at the same time.<br>
The thesis first gives an overview over the current literature. Then, a new method for uplift modeling is proposed called the Difference in Outcome Model (DOM), which is designed specifically with the multiple treatment case in mind. This method is compared to some of the most popular current methods using data from real world marketing campaigns.<br>
The results show that the newly proposed method performs competitively or outperforms the other models in the cases tested in this thesis. However, they also show that the performance of uplift models can be somewhat unstable and is highly dependent on the data.
Therefore, further testing on more and different data sets needs to be done, in order to confirm the good performance of the new model. Based on this thesis it is recommended that practitioners test several models to find the one which performs best on their data. The new method could be one of those models due to its overall good performance. <br>

# Introduction
Uplift modeling is a predictive modeling technique which is concerned with directly estimating the effect of a treatment on a person's behavior. Specifically, it tries to estimate the conditional average treatment effect (CATE) of a treatment based on the characteristics of a given person (age, gender, income, ...). While the type of predictive modeling has applications in a wide variety of field, we will focus on a marketing framework. Therefore, when we talk about "treatments" we usually refer to marketing campaigns, ads, etc. <br>
So far, most research in this field has been done on a single treatment case and on cases with a binary outcome. We aim to fill that gap by looking exclusively at cases where multiple possible treatments are compared and focus on continuous outcome variables.
<br>
# Implemented Models
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

