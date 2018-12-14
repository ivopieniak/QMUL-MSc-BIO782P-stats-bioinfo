#Dataset 1
marine = read.delim('part_1_student_1063.tdf', header = T) 
str(marine)
attach(marine)
summary(marine)
#boxplot diversity vs latitude
div_lat = boxplot(UniFracInd ~ latitude, xlab='Latitdue',
                  ylab='Microbial diversity', col=c('blue','orange'))
div_lat$stats         # distribution stats
#boxplot diversity vs season
div_seas = boxplot(UniFracInd ~ season, xlab='Season',
                   ylab='Microbial diversity', col=c('blue','orange'))
div_seas$stats        #distribution stats
#shapiro-wilk test to check if data is normally distributed
norm_tropical = shapiro.test(UniFracInd[latitude=="tropical"])
norm_temp = shapiro.test(UniFracInd[latitude=="temperate"])
norm_jan = shapiro.test(UniFracInd[season=="Jan"])
norm_aug = shapiro.test(UniFracInd[season=="Aug"]) 
# it is approximately normally *norm_aug
# data is good to do the t-test
#
t_latitude = t.test(UniFracInd ~ latitude)
t_season = t.test (UniFracInd ~ season)


# Q3
#Interaction model between season and latitdue
model_lat_seas = lm(UniFracInd ~ season * latitude)
par(mfrow=c(2,2))
plot(model_lat_seas)

#Assesing quality
drop1(model_lat_seas, test = 'F') 
#Interaction plots 
par(mfrow=c(2,1))
interaction.plot(season,latitude,UniFracInd)
interaction.plot(latitude,season,UniFracInd)
############################################################################
#DATASET 2
#DATASET 2

RNA = read.delim('part_2_student_1063.tdf', header = TRUE) 
str(RNA)
attach(RNA)
# Relationship between expression and distance
plot(expression_fold~distance, xlab='Genetic distance', ylab='RNA expression', pch=19, col='blue')
#Creating a model
rna_model = lm(expression_fold~distance)
#ANOVA to determine signifficant difference in variation 
anova(rna_model)
par(mfrow=c(2,2))
#diagnostic plot
plot(rna_model)
par(mfrow=(c(1,1)))
# how does expression change with distance
plot(distance, expression_fold, xlab= "Genetic distance", ylab="Homologue Expression", pch=19
     , col='blue')
abline(rna_model, col='orange')

####################################################################################
#Dataset 3

#Reading the data + checking the data structre
HIV = read.delim('part_3_student_1063.tdf', header=TRUE)
str(HIV)
summary(HIV)
attach(HIV)
par(mfrow=c(2,2))

#Plots to represent the relationship between viral load and cell count,tissue,diversity, and evolutionary distance
plot(VLoad ~ CD4, xlab = 'CD4+ cell count', ylab='Log10 viral load (particles/ml)', names=c('High','Low'),
     col=c('blue','orange'))
plot(VLoad ~ tissue, xlab = 'Tissue', ylab='Log10 viral load (particles/ml)', names=c('Brain','Spinal Cord'),
     col=c('blue','orange'))
plot(VLoad ~ score_shannon, xlab='Shannon population diversity',ylab='Log10 viral load (particles/ml)',
     pch=19, col='blue')
plot(VLoad ~ score_distance, xlab='Relative genetic distance', ylab='Log10 viral load (particles/ml)',
     pch=19, col='blue')
par(mfrow=c(1,1))
#Automated model selection
hiv_model= lm(VLoad ~ CD4 * tissue * score_shannon * score_distance)
anova(hiv_model)     # data quality assesment
# Automated stepwise regression forward/backward, both
backwards = step(hiv_model,direction ='backward')
forwards = step(lm(VLoad ~ 1), scope=(~ score_shannon* CD4 * tissue * score_distance),
                direction = 'forward')
both = step(lm(VLoad ~ score_shannon + score_distance), scope = c(lower=~score_shannon, upper=~score_shannon*score_distance*tissue))
#Final model determined by the lowest AIC score
final_model = lm(VLoad ~ score_shannon + score_distance)
anova(final_model)
