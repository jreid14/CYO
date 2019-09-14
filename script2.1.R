##loading packages
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(broom)) install.packages("broom", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require(ggcorrplot)) install.packages("ggcorrplot", repos = "http://cran.us.r-project.org")
if(!require(factoextra)) install.packages("factoextra", repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(gsubfn)) install.packages("gsubfn", repos = "http://cran.us.r-project.org")
if(!require(ggvis)) install.packages("ggvis", repos = "http://cran.us.r-project.org")
if(!require(knitr)) install.packages("knitr", repos = "http://cran.us.r-project.org")
if(!require(magrittr)) install.packages("magrittr", repos = "http://cran.us.r-project.org")
if(!require(matrixStats)) install.packages("matrixStats", repos = "http://cran.us.r-project.org")
if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
if(!require(RColorBrewer)) install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
if(!require(rlang)) install.packages("rlang", repos = "http://cran.us.r-project.org")
if(!require(rmarkdown)) install.packages("rmarkdown", repos = "http://cran.us.r-project.org")
if(!require(rpart)) install.packages("rpart", repos = "http://cran.us.r-project.org")
if(!require(GGally)) install.packages("GGally", repos = "http://cran.us.r-project.org")
if(!require(stringr)) install.packages("stringr", repos = "http://cran.us.r-project.org")
if(!require(car)) install.packages("car", repos = "http://cran.us.r-project.org")
if(!require(class)) install.packages("class", repos = "http://cran.us.r-project.org")
if(!require(gbm)) install.packages("gbm", repos = "http://cran.us.r-project.org")
if(!require(randomForest)) install.packages("randomForest", repos = "http://cran.us.r-project.org")


##data import
adult <- read.csv("adult.csv")
str(adult)

table(adult$workclass)


### Functions ####

stat_smooth_func = function(mapping = NULL, data = NULL,
                            geom = "smooth", position = "identity",
                            ...,
                            method = "auto",
                            formula = y ~ x,
                            se = TRUE,
                            n = 80,
                            span = 0.75,
                            fullrange = FALSE,
                            level = 0.95,
                            method.args = list(),
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE,
                            xpos = NULL,
                            ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc = ggproto("StatSmooth", Stat,
                         
                         setup_params = function(data, params) {
                           # Figure out what type of smoothing to do: loess for small datasets,
                           # gam with a cubic regression basis for large data
                           # This is based on the size of the _largest_ group.
                           if (identical(params$method, "auto")) {
                             max_group <- max(table(data$group))
                             
                             if (max_group < 1000) {
                               params$method <- "loess"
                             } else {
                               params$method <- "gam"
                               params$formula <- y ~ s(x, bs = "cs")
                             }
                           }
                           if (identical(params$method, "gam")) {
                             params$method <- mgcv::gam
                           }
                           
                           params
                         },
                         
                         compute_group = function(data, scales, method = "auto", formula = y~x,
                                                  se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                  xseq = NULL, level = 0.95, method.args = list(),
                                                  na.rm = FALSE, xpos=NULL, ypos=NULL) {
                           if (length(unique(data$x)) < 2) {
                             # Not enough data to perform fit
                             return(data.frame())
                           }
                           
                           if (is.null(data$weight)) data$weight <- 1
                           
                           if (is.null(xseq)) {
                             if (is.integer(data$x)) {
                               if (fullrange) {
                                 xseq <- scales$x$dimension()
                               } else {
                                 xseq <- sort(unique(data$x))
                               }
                             } else {
                               if (fullrange) {
                                 range <- scales$x$dimension()
                               } else {
                                 range <- range(data$x, na.rm = TRUE)
                               }
                               xseq <- seq(range[1], range[2], length.out = n)
                             }
                           }
                           # Special case span because it's the most commonly used model argument
                           if (identical(method, "loess")) {
                             method.args$span <- span
                           }
                           
                           if (is.character(method)) method <- match.fun(method)
                           
                           base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                           model <- do.call(method, c(base.args, method.args))
                           
                           m = model
                           eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                            list(a = format(coef(m)[1], digits = 3), 
                                                 b = format(coef(m)[2], digits = 3), 
                                                 r2 = format(summary(m)$r.squared, digits = 3)))
                           func_string = as.character(as.expression(eq))
                           
                           if(is.null(xpos)) xpos = min(data$x)*0.9
                           if(is.null(ypos)) ypos = max(data$y)*0.9
                           data.frame(x=xpos, y=ypos, label=func_string)
                           
                         },
                         
                         required_aes = c("x", "y")
)

draw_confusion_matrix <- function(cm, class1, class2, strTitle) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(paste(strTitle, 'CONFUSION MATRIX'), cex.main=2)
  
  color1 = '#88b56c' # correct
  color2 = '#a83939' # miss
  
  # create the matrix 
  rect(150, 430, 240, 370, col = color1)
  text(195, 435, class1, cex=1.2)
  rect(250, 430, 340, 370, col = color2)
  text(295, 435, class2, cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col = color2)
  rect(250, 305, 340, 365, col = color1)
  text(140, 400, class1, cex=1.2, srt=90)
  text(140, 335, class2, cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  


unemp <- function(job){
  job <- as.character(job)
  if(job=='Never-worked' | job=='Without-pay'){
    return('Unemployed')
  }else{
    return(job)
  }
}

adult$workclass <- sapply(adult$workclass,unemp)

table(adult$workclass)


##combine Local and State government jobs into a category called SL-gov and combine all self-employed jobs into a category called self-emp
group_emp <- function(job){
  if(job=='Local-gov' | job=='State-gov'){
    return('SL-gov')
  }else if(job=='Self-emp-inc' | job=='Self-emp-not-inc'){
    return('self-emp')
  }else{
    return(job)
  }
}

adult$workclass <- sapply(adult$workclass,group_emp)
table(adult$workclass)

##marital status
table(adult$marital.status)

group_marital <- function(mar){
  mar <- as.character(mar)
  #Not-Married
  if(mar=='Divorced' | mar=='Separated' | mar=='Widowed'){
    return('Not-Married')
  }else if(mar=='Never-married'){
    return(mar)
  }else{
    return('Married')
  }
}
adult$marital.status <- sapply(adult$marital.status,group_marital)
table(adult$marital.status)

##Country column
adult = adult %>%
  filter(native.country != 'Holand-Netherlanads') ##remove netherlands as it is an outlier
table(adult$native.country)
levels(adult$native.country)

Asia <- c('China','Hong','India','Iran','Cambodia','Japan','Laos','Philippines','Vietnam','Taiwan','Thailand')

North.America <- c('Canada','United-States','Puerto-Rico' )

Europe <- c('England','France','Germany','Greece','Hungary','Ireland','Italy','Poland','Portugal','Scotland','Yugoslavia')

Latin.and.South.America <- c('Columbia','Cuba','Dominican-Republic','Ecuador','El-Salvador','Guatemala','Haiti','Honduras','Mexico','Nicaragua',
                             'Outlying-US(Guam-USVI-etc)','Peru','Jamaica','Trinadad&Tobago')


Other <- c('South')

group_country <- function(ctry){
  if(ctry %in% Asia){
    return('Asia')
  }else if(ctry %in% North.America){
    return('North.America')
  }else if(ctry %in% Europe){
    return('Europe')
  }else if(ctry %in% Latin.and.South.America){
    return('Latin.and.South.America')
  }else{
    return('Other')
  }
}
adult$native.country <- sapply(adult$native.country,group_country)
table(adult$native.country)

##checking the str() of adult again to make sure any of the columns we changed have factor levels with factor()
str(adult)

##converting to factors
adult$workclass <- sapply(adult$workclass,factor)
adult$native.country <- sapply(adult$native.country,factor)
str(adult)
adult$marital.status <- sapply(adult$marital.status,factor)


#convert all the '?' values to a NA value
adult[adult == '?'] <- NA

table(adult$workclass)
adult$workclass <- sapply(adult$workclass,factor)
adult$native.country <- sapply(adult$native.country,factor)
adult$marital.status <- sapply(adult$marital.status,factor)
adult$occupation <- sapply(adult$occupation,factor)
adult <- na.omit(adult)
str(adult)

#Splitting our data in Training and Testing set
train_s = sample(seq_len(nrow(adult)), size = floor(0.75 * nrow(adult)))

train = adult[train_s, ]
test = adult[-train_s, ]

#Training the Model
model <- glm(income ~ ., family=binomial(logit), data=train)
summary(model)


set.seed(1092)

train_s = sample(seq_len(nrow(adult)), size = floor(0.75 * nrow(adult)))

adult.train = adult[train_s, ]
adult.test = adult[-train_s, ]


adult = rbind(adult.test, adult.train)

all_adult = list(adult.test, adult.train, adult)

# partion to just continuous for correlation testing

adult.train.cont = adult[c("age", "education.num", "capital.gain", "capital.loss", "hours.per.week")]

# create custom order for education

edu_order = (adult %>%
               select(education, education.num) %>%
               distinct(education, education.num) %>%
               arrange(education.num))$education


# bar plot education - income bracket 

ggplot(data = adult, aes(x = education.num, fill = adult$income)) + 
  geom_bar() + scale_x_discrete(labels = edu_order) + 
  xlab("Education") + ylab("Count") + 
  scale_fill_manual(values=c("#999999", "#E69F00"))

# bar plot age - Income Bracket

ggplot(data = adult, aes(x = age, fill = income)) + geom_bar() +
  xlab("Age") + ylab("Count") + 
  scale_fill_manual(values=c("#999999", "#E69F00"))


# bar plot sex - Income Bracket

ggplot(data = adult, aes(x = sex, fill = income)) + geom_bar() + 
  scale_fill_manual(values=c("#999999", "#E69F00"))

edu_by_freq = 
  adult %>%
  group_by(education) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n))

edu_order_number = 
  adult %>%
  mutate(edu = factor(education, levels = edu_order)) %>%
  arrange(edu)

set.seed(4321)

logit.form = income ~ . - fnlwgt - education - native.country - education.num + log(x = education.num) - age + poly(age, 2)

#- age + poly(age, 2) - education

# - Marital.Status - Relationship - native.country - Capital_Gain - Capital_Loss - fnlwgt - education

adult.logit = glm(formula = logit.form, 
               family = binomial(link = "logit"), 
               data = adult.train)

# checking nature of the probabilities from logit model

adult.logit.probs = predict(adult.logit, type = "response")

adult.probs = cbind(adult.train, probs = adult.logit.probs)

# probs and age

adult.probs.age = adult.probs %>%
  select(age, probs) %>%
  group_by(age) %>%
  summarize(total = n(), probs = mean(probs))

# ggplot(adult.probs.age, aes(x = age, y = probs)) + geom_point() + geom_smooth()

# plot of age and probability is erratic when ages are high, this is due to small sample size for ages 85+.  
# A threshold for size of n is necessary to make this plot better.

adult.probs.age = adult.probs.age %>%
  filter(total > 20)

ggplot(adult.probs.age, aes(x = age, y = probs)) + geom_point() + geom_smooth(method = "lm")

# probs and education

adult.probs.edu = adult.probs %>%
  select(education.num, probs) %>%
  group_by(education.num) %>%
  summarize(total = n(), probs = mean(probs))

ggplot(adult.probs.edu, aes(x = education.num, y = probs)) + geom_point() + geom_smooth() +
  xlab("education")

# probs and sex

adult.probs.sex = adult.probs %>%
  select(sex, probs) %>%
  group_by(sex) %>%
  summarize(total = n(), probs = mean(probs))

ggplot(adult.probs, aes(x = sex, y = probs)) + 
  geom_boxplot() + 
  stat_summary(fun.y = "mean", geom = "point", shape = 23, size = 3, fill = "green")


# confusion matrix


adult.test.probs = predict(object = adult.logit, newdata = adult.test, type = "response")

adult.test.probs.str = rep("<=50K", length(adult.test.probs))
adult.test.probs.str[adult.test.probs > .7] = ">50K"
atps<-as.factor(adult.test.probs.str)

table(adult.test.probs.str, adult.test$income)
print(paste(round(mean(adult.test.probs.str == trimws(as.character(adult.test$income))), 3), "% Accuracy"))

cm.logit = confusionMatrix(data = atps, adult.test$income)



draw_confusion_matrix(cm.logit, "<=50K", ">50K", "LOGIT")

summary(adult.logit)


### random forests ####

forest_form = income ~ . - education - fnlwgt - native.country

#- Marital.Status - Relationship - native.country - Capital_Gain - Capital_Loss - fnlwgt - education

set.seed(1234)


adult.rf = randomForest(forest_form, data = adult.train)

importance(adult.rf)

varImpPlot(adult.rf)

adult.rf.test.probs = predict(object = adult.rf, newdata = adult.test, type = "response")
artp<-as.factor(adult.rf.test.probs)

cm.logit = confusionMatrix(data = artp, adult.test$income)

draw_confusion_matrix(cm.logit, "<=50K", "<50K", "RAND FOREST")


### boosted forests ####

set.seed(2314)

adult.train.bernoulli = adult.train
adult.test.bernoulli = adult.test

adult.train.bernoulli$income = ifelse(adult.train.bernoulli$income == ">50K", 1, 0)
adult.test.bernoulli$income = ifelse(adult.test.bernoulli$income == ">50K", 1, 0)

# adult.train.bernoulli$income = as.integer(str_replace(adult.train.bernoulli$income, ">50K", 1))
# adult.train.bernoulli$income[is.na(adult.train.bernoulli$income)] = 0

# adult.test.bernoulli$income = as.integer(str_replace(adult.test.bernoulli$income, ">50K", 1))
# adult.test.bernoulli$income[is.na(adult.test.bernoulli$income)] = 0

adult.boost = gbm(forest_form, 
               data = adult.train.bernoulli, 
               distribution = "bernoulli", 
               n.trees = 5000, 
               interaction.depth = 5,
               shrinkage = .05)


# varImp(adult.boost, n.trees = 5000)



adult.boost.test.probs = predict(object = adult.boost, newdata = adult.test.bernoulli, type = "response", n.trees = 5000)
adult.boost.test.probs[adult.boost.test.probs > .5] = 1
adult.boost.test.probs[adult.boost.test.probs <= .5] = 0

abt<-as.factor(adult.boost.test.probs)
abti<-as.factor(adult.test.bernoulli$income)


cm.logit = confusionMatrix(data = abt, abti)

draw_confusion_matrix(cm.logit, "<=50K", ">50K", "BOOSTED FOREST")

