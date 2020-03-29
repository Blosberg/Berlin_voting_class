# Berlin Voting Classifier

Voting data sets are a convenient platform for data science analysis for several reasons: they are quantifiable measurements of human social psychology, and they exhibit the right balance of noise as well as intuitive trends. Moreover, the data is well controlled. Frequenty social science studies rely on sentiment analysis from Twitter or other social media, and are therefore susceptible to bias from bots and `sock puppet' accounts. Security around ballot boxes ensures each citizen votes no more than once.
In Germany, in particular, the compulsory `Wohnmeldungs' allow for very detailed information in each riding. Moreover, in Berlin (as well as at the national level), there is an unabigious physical boundary separating formerly eastern and western districts; this historical border continues to have political and cultural significance, and serves as useful landmark for binary classification. 

As such, machines can learn to identify the historical region of a voting bezirk just by looking at the way people in that bezirk vote, with over 95% accuracy, given certain allowances. Some regions are better integrated than others, and these regions are overwhelmingly dominated by a particular party. Can you guess which?

Full results are shown in [this report](Writeup/Berlin_report.pdf), though the project is an ongoing work. Feel free to comment.

![](figures/BMap_Ccoded.png) 

