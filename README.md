# Code repository for SAD classification

This is a code repository for the SAD classification project. 

Background:
Seasonal affective disorder (SAD) is characterized by seasonally recurring depression.
At NRU, we have identified several interesting phenotypes of SAD with different imaging modalities. Using fMRI, we have observed lower amygdala response to fear, angry and neutral faces in SAD individuals compared to healthy controls (HC), independent of season (Borgsted et al., in prep.). Using DASB PET, we have seen that SAD patients have higher serotonin transporter binding compared to HC (P = 0.01) in winter (McMahon et al., 2016). Neuropsychological tests of cognition indicate group differences in working memory, cognitive processing speed and motor speed (Hjort et al., in prep.). Finally, there is evidence for significant group differences in neuroticism (Hjort et al. 2, in prep.). However, beyond identifying group differences, determining whether these measures informatively predict group status would support their potential diagnostic value. Current SAD diagnostic criteria include the presence of symptoms for two years, thus identifying markers aiding or facilitating clinical determination could save resources and benefit patients. 

Hypothesis:
We can predict group status (SAD or HC) significantly above chance using informative phenotypic markers. Specific hypotheses:
1.	Can fMRI data/PET data/Neuropsych data individually predict group status?
2.	Can a combination of all the data significantly predict group status?
3.	Which markers (fMRI, PET or neuropsychological) contribute most to prediction accuracy?

Subjects and Methods:
Different subsets of individuals have various measures of interest:
29 SAD and 30 HC individuals have neuropsychological and personality measures at summer/winter.
17 SAD and 23 HC individuals have [11C]DASB PET scans at summer/winter.
14 SAD and 7 HC individuals have [11C]SB207145 PET scans at summer/winter.
17 SAD and 15 HC individuals have MRI scans (including structural MR and fMRI) at summer/winter.

Data analysis:
We will use machine-learning models to predict group status (SAD or HC). This will include logistic regression and randomForest and possibly other model types as necessary. Model performance will be prediction of unseen (test) datasets.

