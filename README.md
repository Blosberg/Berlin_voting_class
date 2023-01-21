# Berlin Voting Classifier

Voting data sets are a convenient platform for data science analysis. Many social science studies rely on sentiment analysis from Twitter or other social media, but are susceptible to bias from bots and 'sock puppet' accounts. 

The security around ballot boxes (and the geographical restrictions applied to ballot collection) ensures much more controlled results.
In Germany, in particular, the compulsory `Wohnmeldungs' allow for very detailed information in each riding. Moreover, in Berlin (as well as at the national level), there is a physical boundary separating formerly eastern and western districts; this historical border continues to have political and cultural significance, and serves as useful landmark for binary classification. 

As such, machines can learn to identify the historical region of such regions (in German 'bezirk's) just by looking at the way people there vote --with over 95% accuracy, given certain allowances. Some regions are better integrated than others, and these regions are overwhelmingly dominated by a particular party. Can you guess which? 

Full results are shown in [this report](Writeup/Berlin_report.pdf). Please note that this project is motivated by personal interest to develop and apply data science tools and associated fields of machine learning. Feel free to comment on remaining questions and features that could still be explored.

![](figures/BMap_Ccoded.png) 

