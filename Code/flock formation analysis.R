
# 1) Prepare network data data ----------------------------------------------------

#install.packages("asnipe")
library(asnipe)


# 1.1) load raw network data ----------------------------------------------

net.data.summer <- read.delim("Data/Mill.data.summer.txt", sep=" ", row.names = 1)
net.data.autumn <- read.delim("Data/Mill.data.autumn.txt", sep=",", row.names = 1)
net.data.winter <- read.delim("Data/Mill.data.winter.txt", sep=",", row.names = 1)
net.data.spring <- read.delim("Data/Mill.data.spring.txt", sep=",", row.names = 1)
# Date.Time = yymmddHHMMSS
# Antenna: A for auxiliary, M for main
# PIT = 10 digit alphanumeric code unique to each individual
# location = Mill1-Mill6 (6 feeders around study area)
# week = experimental week within each season (during each week, we collected data for 48 hours). summer has 14 weeks, all other season 3 weeks


# 1.2) Create group by individual matrices --------------------------------

# NOTE: gmm functions are computationally quite intense - we therefore provide the outputs (gmm objects)
# but we explain the procedures here and provide the code necessary to replicate it

# for summer, we built individual gmm objects per experimental week
# here week 1 as an example, but weeks 2-14 were constructed the same way

# sub.net.1 <- subset(net.data.summer, net.data.summer$week==1)
#   gmm.summer.w1 <- gmmevents(
#     time =  sub.net.1$Date.Time,
#     identity =  sub.net.1$PIT,
#     location =  sub.net.1$location,
#     verbose = TRUE,
#     splitGroups = TRUE,
#     global_ids = unique(net.data.summer$PIT)
#   )
# save(gmm.summer.w1, file="gmm.summer.w1.RData")

# loading the gmm objects per week (summer only)
load("Data/gmm.summer.w1.RData")
gmm.summer.w1 <- gmm.summer.w
load("Data/gmm.summer.w2.RData")
gmm.summer.w2 <- gmm.summer.w
load("Data/gmm.summer.w3.RData")
gmm.summer.w3 <- gmm.summer.w
load("Data/gmm.summer.w4.RData")
gmm.summer.w4 <- gmm.summer.w
load("Data/gmm.summer.w5.RData")
gmm.summer.w5 <- gmm.summer.w
load("Data/gmm.summer.w6.RData")
gmm.summer.w6 <- gmm.summer.w # NOTE: due to computational limitations, this data set only includes great tits!
load("Data/gmm.summer.w7.RData")
gmm.summer.w7 <- gmm.summer.w
load("Data/gmm.summer.w8.RData")
gmm.summer.w8 <- gmm.summer.w
load("Data/gmm.summer.w9.RData")
gmm.summer.w9 <- gmm.summer.w
load("Data/gmm.summer.w10.RData")
gmm.summer.w10 <- gmm.summer.w
load("Data/gmm.summer.w11.RData")
gmm.summer.w11 <- gmm.summer.w
load("Data/gmm.summer.w12.RData")
gmm.summer.w12 <- gmm.summer.w
load("Data/gmm.summer.w13.RData")
gmm.summer.w13 <- gmm.summer.w
load("Data/gmm.summer.w14.RData")
gmm.summer.w14 <- gmm.summer.w


# creating the gmm objects for autumn, winter and spring is straight forward, they contain 3 weeks of data collection each

# gmm.autumn <- gmmevents(
#   time = net.data.autumn$Date.Time,
#   identity = net.data.autumn$PIT,
#   location = net.data.autumn$location,
#   verbose = TRUE,
#   splitGroups = TRUE
# )

#save(gmm.autumn, file="/Data/gmm.autumn.RData")
load(file="Data/gmm.autumn.RData")

# gmm.winter <- gmmevents(
#   time = net.data.winter$Date.Time,
#   identity = net.data.winter$PIT,
#   location = net.data.winter$location,
#   verbose = TRUE,
#   splitGroups = TRUE
# )

#save(gmm.winter, file="/Data/gmm.winter.RData")
load(file="Data/gmm.winter.RData")

# gmm.spring <- gmmevents(
#   time = net.data.spring$Date.Time,
#   identity = net.data.spring$PIT,
#   location = net.data.spring$location,
#   verbose = TRUE,
#   splitGroups = TRUE
# )

#save(gmm.spring, file="/Data/gmm.spring.RData")
load(file="Data/gmm.spring.RData")


# 1.3) load breeding data -------------------------------------------------

fledgling_data <- read.delim("Data/fledgling_data.txt",sep="\t")
# this data frame contains all information on the chicks hatched and their parents
# - Box: name of nest box
# - Ring: number of metal leg ring (unique to individual)
# - Tag: 10-digit alphanumeric code of PIT tag (unique to individual). If NA, bird not tagged
# - Who: 'Chick', 'Male' or 'Female'
# - clutch.size: final egg number before incubation (only values for chicks)
# - Chick.weight: weight in g 15 days after hatching
# - Fledged: Date of fledging in number of days after the 1st of April
# - Fledge.order: The order with which chicks fledged (from faceplating data) - within nest boxes
# - Location: 'Castle' or 'Guett' - this analysis focuses on 'Castle' only
# - Coordinates_long: longitude of nest box location
# - Coordinates_lat: latitude of nest box location
# - Family: numeric identifier of each 'family' - i.e. parents and their chicks

# we subset the fledgling_data to the 'castle' site only
fledgling_data <- subset(fledgling_data, fledgling_data$Location=="Castle")

load(file="Data/distance.RData")
# the object 'dist' is a matrix with Euclidean distances between all nestboxes (in m)

# 1.4) Create lists of great tits only & adults/fledglings ----------------

species.age <- read.delim("Data/species_age.txt",sep="\t")
GT.list <- unique(subset(species.age$Pit, species.age$Species=="GRETI"))
adults <- unique(subset(species.age$Pit, species.age$Age_in_2020=='adult'))
GT.adults <- unique(intersect(GT.list, adults))
juveniles <- unique(subset(species.age$Pit, species.age$Age_in_2020=="juvenile"))
# GT.juveniles are IDs of fledglings
GT.juveniles <- unique(intersect(GT.list, juveniles))
# GT.juveniles.sub contains the IDs of fledglings that hatched in the study area
GT.juveniles.sub <- intersect(GT.juveniles, fledgling_data$Tag)

GT.breeding.pairs <- as.vector(na.omit(unique(subset(fledgling_data$Tag, fledgling_data$Who %in% c("Male", "Female")))))

# GT.adults and GT.juveniles contain PIT codes of adult great tits and juveniles, and adult breeding pairs, respectively
# GT.juveniles.sub are the ones that were hatched inside our nest boxes

# 2) Analysis 1: Ontogeny of relationships during transition to independence -------------------------------------------------------
# How do associations between fledglings and others (siblings, parents, other adults, peers) change with age? Do fledglings inherit their parents'social network?

# 2.1) Prepare data frame for analysis ------------------------------------

create.df.per.week <- function(gmm.data, week){
  # due to compuational reasons, the location for week 6 has extra info - we reduce this again to contain the info about which feeder the group was observed at (Mill1-6)
  if(week==6){
    gmm.data$metadata$Location <- substr(gmm.data$metadata$Location, 1,5)
  }
  
  
  # we here calculate the space use for each individual (i.e. how many times they have used which feeder)
  space_use <- as.data.frame(matrix(0, ncol=6, nrow=length(GT.list)))
  rownames(space_use) <- GT.list
  colnames(space_use) <-   c("Mill1", "Mill2", "Mill3", "Mill4", "Mill5", "Mill6")
  for(i in GT.list){
    # extract how many times they were seen at which feeder
    tab.ID1 <- table(gmm.data$metadata$Location[which(gmm.data$gbi[,which(colnames(gmm.data$gbi)%in% i)]==1)])
    for(j in names(tab.ID1)){
      space_use[i,j] <- tab.ID1[j]
    }
      
  }
# this table is needed further down
  

  # subset to those that were hatched in nest boxes (otherwise we don't have age data)
 # GT.list <- intersect(GT.list, fledgling_data$Tag)

    net <- get_network(
    association_data = gmm.data$gbi,
    data_format = "GBI",
    association_index = "SRI",
    times = gmm.data$metadata$Start,
    identities = colnames(gmm.data$gbi),
    which_identities = GT.list
  )
  
  IDs <- rownames(net)
  comb.IDs <- combn(IDs,2, simplify=TRUE)
  
  # note that we double up the data essentially, since we want the data from the perspective of both individuals
  df.net <- rbind.data.frame(cbind(comb.IDs[1,], comb.IDs[2,]),cbind(comb.IDs[2,], comb.IDs[1,]))
  colnames(df.net) <- c("ID1", "ID2")
  
  df.net.comb <- df.net
  df.net.comb$week <- week
  
  for(i in 1:length(df.net.comb[,1])){
    ID1 <- df.net.comb[i, "ID1"]
    ID2 <- df.net.comb[i, "ID2"]
    # extract the number of times ID1 and ID2 have been seen in the same group
    
    # ID1.ID2.together <- length(which(rowSums(gmm.data$gbi[,which(colnames(gmm.data$gbi)%in% c(ID1, ID2))])==2))
    # 
     ID1.ID2.together <- net[ID1, ID2]
    
    df.net.comb[i, "assoc"] <- ID1.ID2.together
    
    # extract how often ID1 and ID2 were seen
    num.sightings.ID1 <- sum(gmm.data$gbi[,ID1])
    num.sightings.ID2 <- sum(gmm.data$gbi[,ID2])    

    df.net.comb[i, "num.sightings.ID1"] <- num.sightings.ID1
    df.net.comb[i, "num.sightings.ID2"] <- num.sightings.ID2

        
    # if same age
    age.ID1 <- unique(subset(species.age$Age_in_2020, species.age$Pit==ID1))
    age.ID2 <- unique(subset(species.age$Age_in_2020, species.age$Pit==ID2))
    
    # plus add the age and species of ID1
    df.net.comb[i, "age.ID1"] <- age.ID1
    
    # relationship (if siblings, peers (non-sibs), parent-offspring, adult-offspring, adult-adult)
    
    if(age.ID1 == "adult" & age.ID2 == "adult"){
      relationship <- "adult_adult"
    } else if(age.ID1 == "adult" & age.ID2 == "juvenile" |age.ID1 == "juvenile" & age.ID2 == "adult") {
      # check if it's parent-offspring or adult-juvenile
      box.ID1 <- subset(fledgling_data$Box, fledgling_data$Tag==ID1)
      box.ID2 <- subset(fledgling_data$Box, fledgling_data$Tag==ID2)
      if(length(box.ID1)==0){
        box.ID1 <- "Nobox"
      }
      if(length(box.ID2)==0){
        box.ID2 <- "Nobox"
      }
      if(length(intersect(box.ID1,box.ID2))>0){
        relationship <- "parent_offspring"
      } else {
        relationship <- "adult_juvenile"
      }
    } else if(age.ID1 == "juvenile" & age.ID2 == "juvenile"){
      box.ID1 <- subset(fledgling_data$Box, fledgling_data$Tag==ID1)
      box.ID2 <- subset(fledgling_data$Box, fledgling_data$Tag==ID2)

      if(length(box.ID1)==0){
        box.ID1 <- "Nobox"
      }
      if(length(box.ID2)==0){
        box.ID2 <- "Nobox"
      }
      
      if(box.ID1==box.ID2){
        relationship <- "siblings"
      } else {
        relationship <- "peers"
      }  
      
      
    } 
    df.net.comb[i, "relationship"] <- relationship
    
    # next we extract the age of ID1 (only if it a fledgling) in days since fledging
    # we take the first recording time of that week as the time
    if(age.ID1=="juvenile" & ID1 %in% GT.juveniles.sub){
      fledge.date.ID1 <- subset(fledgling_data$Fledged, fledgling_data$Tag==ID1)
      fledge.date.ID1.date <- as.Date("200401", format="%y%m%d")+fledge.date.ID1
      days.since.fl.ID1 <- as.numeric(as.Date(substr(gmm.data$metadata$Start[1], 1, 6), format="%y%m%d")-fledge.date.ID1.date)
      df.net.comb[i, "time.since.fl"] <- days.since.fl.ID1
    } else {
      df.net.comb[i, "time.since.fl"] <- NA
    }
    
# finally, we extract the space use similarity between individuals 1 and 2
    space.use.ID1 <- space_use[rownames(space_use) %in% c(ID1),]
    space.use.ID2 <- space_use[rownames(space_use) %in% c(ID2),]
    
    time.spent.min <- NULL
    
    for(k in 1:6){
      # at each feeder location, extract how much time they have each spent there
      feederk.ID1 <- space.use.ID1[k]
      feederk.ID2 <- space.use.ID2[k]
      
      time.spent.min[k] <- min(unlist(c(feederk.ID1, feederk.ID2)))
    }
    
    # divide by the total number of visits by individual ID1
    space.overlap <- sum(time.spent.min) / sum(space.use.ID1) 
    df.net.comb[i, "space_overlap"] <- space.overlap
    }
   
  
  # we subset to fledglings only which hatched in our nest boxes 
  df.net.comb <- subset(df.net.comb, df.net.comb$ID1 %in% GT.juveniles.sub)
  
  # finally, subset to only those seen a minimum of 5 times in this week
  df.net.comb <- subset(df.net.comb, df.net.comb$num.sightings.ID1>=5 & df.net.comb$num.sightings.ID2>=5)
  
return(df.net.comb)
  }


# we run the function on each gmm object (week 1-14)
df.summer.week1 <- create.df.per.week(gmm.data=gmm.summer.w1, week=1)
df.summer.week2 <- create.df.per.week(gmm.data=gmm.summer.w2, week=2)
df.summer.week3 <- create.df.per.week(gmm.data=gmm.summer.w3, week=3)
df.summer.week4 <- create.df.per.week(gmm.data=gmm.summer.w4, week=4)
df.summer.week5 <- create.df.per.week(gmm.data=gmm.summer.w5, week=5)
df.summer.week6 <- create.df.per.week(gmm.data=gmm.summer.w6, week=6)
df.summer.week7 <- create.df.per.week(gmm.data=gmm.summer.w7, week=7)
df.summer.week8 <- create.df.per.week(gmm.data=gmm.summer.w8, week=8)
df.summer.week9 <- create.df.per.week(gmm.data=gmm.summer.w9, week=9)
df.summer.week10 <- create.df.per.week(gmm.data=gmm.summer.w10, week=10)
df.summer.week11 <- create.df.per.week(gmm.data=gmm.summer.w11, week=11)
df.summer.week12 <- create.df.per.week(gmm.data=gmm.summer.w12, week=12)
df.summer.week13 <- create.df.per.week(gmm.data=gmm.summer.w13, week=13)
df.summer.week14 <- create.df.per.week(gmm.data=gmm.summer.w14, week=14)

# combining all data into one data frame
df.summer.all.week <- rbind(df.summer.week1, 
                            df.summer.week2, 
                            df.summer.week3, 
                            df.summer.week4, 
                            df.summer.week5, 
                            df.summer.week6, 
                            df.summer.week7, 
                            df.summer.week8, 
                            df.summer.week9, 
                            df.summer.week10,
                            df.summer.week11,
                            df.summer.week12,
                            df.summer.week13,
                            df.summer.week14)
head(df.summer.all.week)
dim(df.summer.all.week)

#save(df.summer.all.week, file="Data/df.summer.all.week.RDA")
load("Data/df.summer.all.week.RDA")
# this loads the data frame that is generated with the function above

# how many dyads and how many individual birds?
length(df.summer.all.week[,1])
# 17188 dyads
length(which(unique(c(df.summer.all.week$ID1, df.summer.all.week$ID2)) %in% GT.juveniles))
# 68 juveniles
length(which(unique(c(df.summer.all.week$ID1, df.summer.all.week$ID2)) %in% GT.adults))
# 60 adults

# 2.2) Run Bayesian linear regression -------------------------------------


# install the packages needed to run a Bayesian regression analysis
#install.packages("rstan")
#install.packages("brms")
library(brms)
library(rstan)

# test for multi-collinearity among predictor variables:

library(car)

model.vif <- lm( assoc ~ relationship + scale(time.since.fl) + space_overlap, data = df.summer.all.week)

vif(model.vif)
#                         GVIF Df GVIF^(1/(2*Df))
# relationship         1.077081  3        1.012453
# scale(time.since.fl) 1.019991  1        1.009946
# space_overlap        1.074557  1        1.036609

# variance inflation factors are all below 5 - we can use predictors in the same model

# Note that brms automatically centers predictor variables to improve model convergence! 

# Model 1
model.fl.formation <-
  brm(
   assoc ~ relationship + scale(time.since.fl) + space_overlap + relationship:scale(time.since.fl) + (1 |mm(ID1, ID2)) + (1 |ID1),
    df.summer.all.week,
    family = zero_inflated_beta(),
    chains=4,
    iter=4000,
    cores = 6
  )

save(model.fl.formation, file="Output/Analysis 1/model.fl.formation.RDA")
load("Output/Analysis 1/model.fl.formation.RDA")

# 2.3) Check model performance --------------------------------------------

# we first asses the performance of the model (check for mixing, stationarity)
plot(model.fl.formation)
# and posterior predictive checks
pp_check(model.fl.formation, ndraws= 1e2)


# 2.4) Interpret the effects ----------------------------------------------

# Interpret

summary(model.fl.formation)

# Family: zero_inflated_beta 
# Links: mu = logit; phi = identity; zi = identity 
# Formula: assoc ~ relationship + scale(time.since.fl) + space_overlap + relationship:scale(time.since.fl) + (1 | mm(ID1, ID2)) + (1 | ID1) 
# Data: df.summer.all.week (Number of observations: 17188) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Group-Level Effects: 
#   ~ID1 (Number of levels: 65) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.16      0.03     0.11     0.21 1.00     1966     3585
# 
# ~mmID1ID2 (Number of levels: 128) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.50      0.04     0.43     0.60 1.00     2188     3357
# 
# Population-Level Effects: 
#                                                 Estimate Est.Error l-95% CI u-95% CI Rhat
# Intercept                                          -2.63      0.06    -2.74    -2.51 1.00
# relationshipparent_offspring                        0.05      0.06    -0.07     0.17 1.00
# relationshippeers                                   0.10      0.05    -0.01     0.21 1.00
# relationshipsiblings                                0.21      0.06     0.08     0.33 1.00
# scaletime.since.fl                                  0.14      0.01     0.12     0.17 1.00
# space_overlap                                       0.80      0.03     0.74     0.86 1.00
# relationshipparent_offspring:scaletime.since.fl    -0.03      0.05    -0.13     0.07 1.00
# relationshippeers:scaletime.since.fl               -0.08      0.02    -0.11    -0.05 1.00
# relationshipsiblings:scaletime.since.fl            -0.21      0.03    -0.27    -0.14 1.00
# Bulk_ESS Tail_ESS
# Intercept                                           1549     2903
# relationshipparent_offspring                       10765     6161
# relationshippeers                                   1326     2151
# relationshipsiblings                                1771     3335
# scaletime.since.fl                                  5498     6177
# space_overlap                                       5316     5857
# relationshipparent_offspring:scaletime.since.fl     9976     6624
# relationshippeers:scaletime.since.fl                6187     5781
# relationshipsiblings:scaletime.since.fl             9109     5835
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# phi    14.59      0.22    14.16    15.02 1.00     9957     5896
# zi      0.45      0.00     0.45     0.46 1.00    13669     6209
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


# These effect sizes in the summary output are all in log odds. 
# We transform these into odds by exponentiating (we report odds in the manuscript)
exp(fixef(model.fl.formation))

#                                                   Estimate Est.Error       Q2.5      Q97.5
# Intercept                                       0.07237722  1.061627 0.06432874 0.08120151
# relationshipparent_offspring                    1.05167231  1.063914 0.93014157 1.18324304
# relationshippeers                               1.10718195  1.055595 0.99421788 1.23021083
# relationshipsiblings                            1.23009066  1.064227 1.08640779 1.38806808
# scaletime.since.fl                              1.15226704  1.013548 1.12235402 1.18296669
# space_overlap                                   2.22489143  1.032453 2.08733283 2.36803248
# relationshipparent_offspring:scaletime.since.fl 0.97215843  1.053700 0.87647316 1.07710672
# relationshippeers:scaletime.since.fl            0.92102105  1.015923 0.89344942 0.95001137
# relationshipsiblings:scaletime.since.fl         0.81101010  1.033986 0.76002347 0.86588021


# we extract and plot effect sizes for the global model

mod.pred_flock.summer <- plot(conditional_effects(model.fl.formation, prob = 0.95))

library(ggpubr)

# Plot Figure:

ggarrange(
  mod.pred_flock.summer$relationship +
    labs(x= "Relationship", y= "association strength (SRI)") +
    ylim(c(0.02,0.10))+
    theme_bw()+
    scale_x_discrete(labels=c("adult_juvenile" = "adults", "parent_offspring" = "parents",
                                        "peers" = "peers", "siblings" = "siblings")),

    mod.pred_flock.summer$time.since.fl +
    labs(x= "Time since fledging [days]", y= "association strength (SRI)") +
    labs(y= "") +
    ylim(c(0.02,0.10))+
    theme_bw()+
    geom_line(color="black", lwd=0.8),

  mod.pred_flock.summer$space_overlap +
    labs(x= "Space use overlap") +
    labs(y= "association strength (SRI)") +
    ylim(c(0.02,0.10))+
    theme_bw()+
    geom_line(color="black", lwd=0.8),
  
  mod.pred_flock.summer$`time.since.fl:relationship` +
    labs(x= "Time since fledging [days]", y= "association strength (SRI)") +
    labs(y= "") +
    ylim(c(0.02,0.10))+
    theme_bw()+
    theme(legend.background=element_blank())+
    theme(legend.title=element_blank(), legend.position = c(0.78, 0.2))+
    scale_fill_brewer(palette="RdYlBu", breaks=c("adult_juvenile", "parent_offspring", "peers", "siblings"), labels=c("adults", "parents", "peers", "siblings"))+
    scale_color_brewer(palette="RdYlBu", breaks=c("adult_juvenile", "parent_offspring", "peers", "siblings"), labels=c("adults", "parents", "peers", "siblings")),
  ncol=2,
  nrow=2,
#legend="right",
  common.legend = FALSE,
  labels= c("a", "b", "c", "d"))

ggsave("Output/Figures/Figure_Analysis1.tiff", units="in", width=8, height=7, dpi=300, compression = 'lzw')


# 3) Analysis 2: Flock formation - phenotype (fledglings only) -----------------------------------------

# Which factors predict associations among non-sibling fledglings after they have reached independence

# for this analysis, we use the summer data set

# 3.0) --------------------------------------------------------------------

# In prep for the analysis on inheritance of the parents' social network, we prepare data to calculate association strengths among breeding adults only.

breeders.gmm <- function(gmm, week){
  # we create one empty group by indiviudal matrix for the breeding adults
  gbi.summer.comb <- matrix(0, ncol=length(unique(c(GT.adults, GT.breeding.pairs)))+1)
  colnames(gbi.summer.comb) <- c("week", unique(c(GT.breeding.pairs, GT.adults)))
  # for each group in the group by individual matrix
  for(i in 1:nrow(gmm$gbi)){
    # add a row to the matrix
    gbi.summer.comb <- rbind(gbi.summer.comb, rep(0, length(unique(c(GT.breeding.pairs, GT.adults)))+1))
    gbi.summer.comb[nrow(gbi.summer.comb),] <- 0
    # we save information on which week it was
    gbi.summer.comb[nrow(gbi.summer.comb),"week"] <- week
    # we extract the IDs of who was present in that group
    who.present <- names(which(gmm$gbi[i,]==1))
    who.present <- subset(who.present, who.present %in% unique(c(GT.breeding.pairs, GT.adults)))
    # and add an entry of 1 to the same IDs in our newly created group by individual matrix
    if(length(who.present)>0){
      gbi.summer.comb[nrow(gbi.summer.comb),who.present] <- 1
      
    }
  }
  gbi.summer.comb <- subset(gbi.summer.comb, gbi.summer.comb[,1]!=0)
  return(gbi.summer.comb)
}

# we run this function across all weeks
breeders.w1 <- breeders.gmm(gmm=gmm.summer.w1, week=1)
breeders.w2 <- breeders.gmm(gmm=gmm.summer.w2, week=2)
breeders.w3 <- breeders.gmm(gmm=gmm.summer.w3, week=3)
breeders.w4 <- breeders.gmm(gmm=gmm.summer.w4, week=4)
breeders.w5 <- breeders.gmm(gmm=gmm.summer.w5, week=5)
breeders.w6 <- breeders.gmm(gmm=gmm.summer.w6, week=6)
breeders.w7 <- breeders.gmm(gmm=gmm.summer.w7, week=7)
breeders.w8 <- breeders.gmm(gmm=gmm.summer.w8, week=8)
breeders.w9 <- breeders.gmm(gmm=gmm.summer.w9, week=9)
breeders.w10 <- breeders.gmm(gmm=gmm.summer.w10, week=10)
breeders.w11 <- breeders.gmm(gmm=gmm.summer.w11, week=11)
breeders.w12 <- breeders.gmm(gmm=gmm.summer.w12, week=12)
breeders.w13 <- breeders.gmm(gmm=gmm.summer.w13, week=13)
breeders.w14 <- breeders.gmm(gmm=gmm.summer.w14, week=14)

# and combine it into one data frame
gbi.breeders <- rbind(breeders.w1,
                      breeders.w2,
                      breeders.w3,
                      breeders.w4,
                      breeders.w5,
                      breeders.w6,
                      breeders.w7,
                      breeders.w8,
                      breeders.w9,
                      breeders.w10,
                      breeders.w11,
                      breeders.w12,
                      breeders.w13,
                      breeders.w14)

# this now contains a group by individual matrix for adults only

# 3.1) prepare gmm data ---------------------------------------------------
# this function extracts associaiton strengths and individual-level covariates (age similarity, weight similarity, association among parents) for each fledgling dyad

create.df.per.week.fledglings <- function(gmm.data, week){
  # due to compuational reasons, the location for week 6 has extra info - we reduce this again to contain the info about which feeder the group was observed at (Mill1-6)
  if(week==6){
    gmm.data$metadata$Location <- substr(gmm.data$metadata$Location, 1,5)
  }
  
  # create network
  IDs.10sightings <- names(which(colSums(gmm.data$gbi)>=5))

  # we here calculate the space use for each individual (i.e. how many times they have used which feeder)
  space_use <- as.data.frame(matrix(0, ncol=6, nrow=length(IDs.10sightings)))
  rownames(space_use) <- IDs.10sightings
  colnames(space_use) <-   c("Mill1", "Mill2", "Mill3", "Mill4", "Mill5", "Mill6")
  for(i in IDs.10sightings){
    # extract how many times they were seen at which feeder
    tab.ID1 <- table(gmm.data$metadata$Location[which(gmm.data$gbi[,which(colnames(gmm.data$gbi)%in% i)]==1)])
    for(j in names(tab.ID1)){
      space_use[i,j] <- tab.ID1[j]
    }
    
  }
  # this table is needed further down
  
  net <- get_network(
    association_data = gmm.data$gbi,
    data_format = "GBI",
    association_index = "SRI",
    times = gmm.data$metadata$Start,
    identities = colnames(gmm.data$gbi),
    which_identities = IDs.10sightings
  )
  
  # we create a network among breeding adults (testing for inheritance of social network)
  # we use the data from the start of data collection up until the respective week (e.g. for week 4 we use weeks 1-4)
  breeders.sub <- subset(gbi.breeders, gbi.breeders[,1] %in% c(1:week))
  # remove the first column (week)
  breeders.sub <- breeders.sub[,2:ncol(breeders.sub)]
  # create the social network
  net.breeders <-   get_network(
    association_data = breeders.sub,
    data_format = "GBI",
    association_index = "SRI",
#    times = gmm.data$metadata$Start,
    identities = colnames(breeders.sub),
    which_identities = intersect(c(GT.breeding.pairs, GT.adults),net.data.summer$PIT)
  )
  
  IDs <- rownames(net)
  comb.IDs <- combn(IDs,2, simplify=TRUE)
  
  # note that we double up the data essentially, since we want the data from the perspective of both individuals
  df.net <- rbind.data.frame(cbind(comb.IDs[1,], comb.IDs[2,]),cbind(comb.IDs[2,], comb.IDs[1,]))
  colnames(df.net) <- c("ID1", "ID2")
  
  df.net.comb <- df.net
  df.net.comb$week <- week
  
  # subset to juveniles as focal individuals
  
  df.net.comb <- subset(df.net.comb, df.net.comb$ID1 %in% GT.juveniles.sub)
  
  for(i in 1:length(df.net.comb[,1])){
    ID1 <- df.net.comb[i, "ID1"]
    ID2 <- df.net.comb[i, "ID2"]
    
    
    # extract the association strength of ID1 and ID2
    
    ID1.ID2.together <- net[ID1, ID2]
    
    df.net.comb[i, "assoc"] <- ID1.ID2.together
    
    
    # we extract the space use similarity between individuals 1 and 2
    space.use.ID1 <- space_use[rownames(space_use) %in% c(ID1),]
    space.use.ID2 <- space_use[rownames(space_use) %in% c(ID2),]
    
    time.spent.min <- NULL
    
    for(k in 1:6){
      # at each feeder location, extract how much time they have each spent there
      feederk.ID1 <- space.use.ID1[k]
      feederk.ID2 <- space.use.ID2[k]
      
      time.spent.min[k] <- min(unlist(c(feederk.ID1, feederk.ID2)))
    }
    
    # divide by the total number of visits by individual ID1
    space.overlap <- sum(time.spent.min) / sum(space.use.ID1) 
    df.net.comb[i, "space_overlap"] <- space.overlap
    
    # extract ID1's age in days in fledging
    
    # next we extract the age of ID1 in days since fledging
    # we take the first recording time of that week as the time
    fledge.date.ID1 <- subset(fledgling_data$Fledged, fledgling_data$Tag==ID1)
    fledge.date.ID1.date <- as.Date("200401", format="%y%m%d")+fledge.date.ID1
    days.since.fl.ID1 <- as.numeric(as.Date(substr(gmm.data$metadata$Start[1], 1, 6), format="%y%m%d")-fledge.date.ID1.date)
    df.net.comb[i, "time.since.fl.ID1"] <- days.since.fl.ID1
    
    # extract the parents of ID1
    box.ID1 <- subset(fledgling_data$Box, fledgling_data$Tag==ID1)
    parents.ID1 <- subset(fledgling_data$Tag, fledgling_data$Box==box.ID1 & fledgling_data$Who %in% c("Male", "Female"))
    
    # if both individuals are juveniles
    # we extract similarities in phenotypes
    # and the association rates of their parents
    if(ID2 %in% GT.juveniles.sub){
      box.ID1 <- subset(fledgling_data$Box, fledgling_data$Tag==ID1)
      box.ID2 <- subset(fledgling_data$Box, fledgling_data$Tag==ID2)
      
      # extract the parents of both IDs
      parents.ID1 <- subset(fledgling_data$Tag, fledgling_data$Box==box.ID1 & fledgling_data$Who %in% c("Male", "Female"))
      parents.ID2 <- subset(fledgling_data$Tag, fledgling_data$Box==box.ID2 & fledgling_data$Who %in% c("Male", "Female"))
      
      
      sum.assoc.among.parents <- sum(net.breeders[rownames(net.breeders) %in% c(parents.ID1), colnames(net.breeders) %in% c( parents.ID2)])
      
      df.net.comb[i, "sum.assoc.parents"] <- sum.assoc.among.parents
      df.net.comb[i, "type"] <- "peers"
      # extract age 
      age.diff <- abs(subset(fledgling_data$Fledged, fledgling_data$Tag==ID1)-subset(fledgling_data$Fledged, fledgling_data$Tag==ID2))
      df.net.comb[i,"age.diff"] <- age.diff
      
      # extract sibling status
      if(box.ID1==box.ID2){
        df.net.comb[i,"siblings"] <- "yes"
      } else {
        df.net.comb[i,"siblings"] <- "no"
      }
      
      # extract fledge weight of ID1 and ID2
      weight.ID1 <- subset(fledgling_data$Chick.weight, fledgling_data$Tag==ID1)
      weight.ID2 <- subset(fledgling_data$Chick.weight, fledgling_data$Tag==ID2)
      
      weight.diff <- abs(weight.ID1-weight.ID2)
      
      df.net.comb[i, "weight.diff"] <- weight.diff
      
      # and add the weight of ID1
      df.net.comb[i, "weight.ID1"] <- weight.ID1
      
      # extract relative fledge order - for those without fledge data, we just take 0
      
      fledge.order.ID1 <- subset(fledgling_data$Fledge.order, fledgling_data$Tag==ID1)
      fledge.order.ID2 <- subset(fledgling_data$Fledge.order, fledgling_data$Tag==ID2)
      
      # extract max weight of the their clutch
      max.fledge.order.box.ID1 <- max(na.omit(subset(fledgling_data$Fledge.order, fledgling_data$Box==box.ID1)))
      max.fledge.order.box.ID2 <- max(na.omit(subset(fledgling_data$Fledge.order, fledgling_data$Box==box.ID2)))
      
      if(is.infinite(max.fledge.order.box.ID1)){
        max.fledge.order.box.ID1 <- NA
      }
      
      if(is.infinite(max.fledge.order.box.ID2)){
        max.fledge.order.box.ID2 <- NA
      }
      
      # now extract the relative fledge order
      
      rel.fledge.order.diff <- abs(fledge.order.ID1/max.fledge.order.box.ID1-fledge.order.ID2/max.fledge.order.box.ID2)
      
      
      df.net.comb[i, "rel.fledge.order.diff"] <- rel.fledge.order.diff
      
      # add the relative fledge order of ID1
      df.net.comb[i, "rel.fledge.order.ID1"] <- fledge.order.ID1/max.fledge.order.box.ID1
      
      # this is only if ID2 is an adult
    #  df.net.comb[i, "sum.assoc.with.parents"] <- NA
    
      
    } else if(ID2 %in% c(GT.adults, GT.breeding.pairs) & !(ID2 %in% parents.ID1)){ 
      # if only ID1 is a juvenile and the other an adult (and not ID1s parents)
   #   df.net.comb[i, "sum.assoc.among.parents"] <- NA
      df.net.comb[i, "age.diff"] <- NA
      df.net.comb[i, "siblings"] <- NA
      df.net.comb[i, "weight.diff"] <- NA
      df.net.comb[i, "weight.ID1"] <- NA
      df.net.comb[i, "rel.fledge.order.diff"] <- NA
      df.net.comb[i, "rel.fledge.order.ID1"] <- NA
      
      # extract the association betetween ID2 and the parents of ID1
      
      # extract the parents of both IDs
      box.ID1 <- subset(fledgling_data$Box, fledgling_data$Tag==ID1)
      parents.ID1 <- subset(fledgling_data$Tag, fledgling_data$Box==box.ID1 & fledgling_data$Who %in% c("Male", "Female"))
    
      sum.assoc.w.parents <- sum(net.breeders[rownames(net.breeders) %in% c(parents.ID1), colnames(net.breeders) %in% ID2])
      
      
      df.net.comb[i, "sum.assoc.parents"] <- sum.assoc.w.parents
      df.net.comb[i, "type"] <- "adult"
      
      # relationships with parents
    } else if(ID2 %in% c(GT.adults, GT.breeding.pairs) & ID2 %in% parents.ID1){
      
      df.net.comb[i, "sum.assoc.parents"] <- NA
      df.net.comb[i, "type"] <- NA
      df.net.comb[i, "age.diff"] <- NA
      df.net.comb[i, "siblings"] <- NA
      df.net.comb[i, "weight.diff"] <- NA
      df.net.comb[i, "weight.ID1"] <- NA
      df.net.comb[i, "rel.fledge.order.diff"] <- NA
      df.net.comb[i, "rel.fledge.order.ID1"] <- NA
   #   df.net.comb[i, "sum.assoc.with.parents"] <- NA
    }

 

  }
  return(df.net.comb)
}

# we run the function on each gmm object (week 1-14)


# 3.2) Prepare data set ---------------------------------------------------



# NOTE: there are no fledglings present in weeks 1-3, hence we exclude those weeks
#df.summer.week1 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w1, week=1))
#df.summer.week2 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w2, week=2))
#df.summer.week3 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w3, week=3))
df.summer.week4 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w4, week=4))
df.summer.week5 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w5, week=5))
df.summer.week6 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w6, week=6))
df.summer.week7 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w7, week=7))
df.summer.week8 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w8, week=8))
df.summer.week9 <-  suppressWarnings(create.df.per.week.fledglings(gmm.data=gmm.summer.w9, week=9))
df.summer.week10 <- suppressWarnings( create.df.per.week.fledglings(gmm.data=gmm.summer.w10, week=10))
df.summer.week11 <- suppressWarnings( create.df.per.week.fledglings(gmm.data=gmm.summer.w11, week=11))
df.summer.week12 <- suppressWarnings( create.df.per.week.fledglings(gmm.data=gmm.summer.w12, week=12))
df.summer.week13 <- suppressWarnings( create.df.per.week.fledglings(gmm.data=gmm.summer.w13, week=13))
df.summer.week14 <- suppressWarnings( create.df.per.week.fledglings(gmm.data=gmm.summer.w14, week=14))

# combining all data into one data frame
df.summer.all.week.fledgies <- rbind(#df.summer.week1, 
                            #df.summer.week2, 
                            #df.summer.week3, 
                            df.summer.week4, 
                            df.summer.week5, 
                            df.summer.week6, 
                            df.summer.week7, 
                            df.summer.week8, 
                            df.summer.week9, 
                            df.summer.week10,
                            df.summer.week11,
                            df.summer.week12,
                            df.summer.week13,
                            df.summer.week14)

#save(df.summer.all.week.fledgies, file="Data/df.summer.all.week.fledgies.RData")
load("Data/df.summer.all.week.fledgies.RData")

# create a df for non-sibling juvenile relationships only removing siblings
df.summer.all.week.fledgies <- subset(df.summer.all.week.fledgies, df.summer.all.week.fledgies$siblings =="no" ) # will be for model 2

head(df.summer.all.week.fledgies)

length(df.summer.all.week.fledgies$assoc)
length(unique(df.summer.all.week.fledgies$ID1))
# 8906 dyads among 65 non-sibling fledgies

# test for multi-collinearity among predictor variables:

library(car)

model.vif <- lm( assoc ~ age.diff + space_overlap + weight.diff , data = df.summer.all.week.fledgies)

vif(model.vif)

# age.diff space_overlap   weight.diff 
# 1.003159      1.004423      1.003654 

# 3.3) Run brm model ------------------------------------------------------

library(brms)
library(rstan)


# we test whether weight and age difference predict associations, while controlling for space use

# Model 2
model.phenotype <-
  brms::brm(
    assoc ~  scale(age.diff) + 
      space_overlap +
      scale(weight.diff) + 
      scale(age.diff):scale(time.since.fl.ID1) +
      scale(weight.diff):scale(time.since.fl.ID1) +
      (1 |mm(ID1, ID2)) + (1 |ID1) ,
    df.summer.all.week.fledgies,
    family = zero_inflated_beta(),
    chains=4,
    iter=6000,
    cores = 6
  )

#save(model.phenotype, file="Output/Analysis 2/model.fl.formation.phenotype.RDA")
load("Output/Analysis 2/model.fl.formation.phenotype.RDA")

# 3.4) Look at model performance ------------------------------------------

plot(model.phenotype)

# and posterior predictive checks
pp_check(model.phenotype, ndraws= 1e2)


# 3.5) Look at effect sizes -----------------------------------------------

summary(model.phenotype)
# Family: zero_inflated_beta 
# Links: mu = logit; phi = identity; zi = identity 
# Formula: assoc ~ scale(age.diff) + space_overlap + scale(weight.diff) + scale(age.diff):scale(time.since.fl.ID1) + scale(weight.diff):scale(time.since.fl.ID1) + (1 | mm(ID1, ID2)) + (1 | ID1) 
# Data: df.summer.all.week.fledgies (Number of observations: 8906) 
# Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 12000
# 
# Group-Level Effects: 
#   ~ID1 (Number of levels: 65) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.19      0.04     0.12     0.27 1.00     3502     6172
# 
# ~mmID1ID2 (Number of levels: 65) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.60      0.07     0.48     0.75 1.00     4921     7612
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# Intercept                                  -2.55      0.09    -2.73    -2.38 1.00     2962
# scaleage.diff                               0.01      0.02    -0.03     0.04 1.00    13406
# space_overlap                               0.83      0.05     0.74     0.91 1.00     8623
# scaleweight.diff                            0.01      0.01    -0.01     0.04 1.00    20558
# scaleage.diff:scaletime.since.fl.ID1       -0.05      0.01    -0.07    -0.02 1.00     7548
# scaleweight.diff:scaletime.since.fl.ID1    -0.01      0.01    -0.03     0.01 1.00    20929
# Tail_ESS
# Intercept                                   5497
# scaleage.diff                               9518
# space_overlap                               9083
# scaleweight.diff                            9506
# scaleage.diff:scaletime.since.fl.ID1        9517
# scaleweight.diff:scaletime.since.fl.ID1     8801
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# phi    13.01      0.26    12.51    13.52 1.00    20226     8336
# zi      0.39      0.01     0.38     0.40 1.00    22492     8457
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# We transform these into odds ratios
exp(fixef(model.phenotype))

# Estimate Est.Error       Q2.5      Q97.5
# Intercept                               0.07788493  1.092704 0.06550054 0.09249411
# scaleage.diff                           1.00925954  1.017571 0.97509184 1.04449842
# space_overlap                           2.28248452  1.046832 2.08658581 2.49219030
# scaleweight.diff                        1.01334401  1.012876 0.98825883 1.03890777
# scaleage.diff:scaletime.since.fl.ID1    0.95240908  1.012863 0.92818883 0.97581980
# scaleweight.diff:scaletime.since.fl.ID1 0.98893812  1.009474 0.97092892 1.00734307


# 3.6) Create a figure ----------------------------------------------------

# setting the conditions for plotting
int.conditions <- list(
  time.since.fl.ID1 = setNames(c(10,30,60), c(10,30,60)),
  age.diff = setNames(c(0,20,40), c(0,20,40)),
  weight.diff = setNames(c(1,3,5), c(1,3,5))
)


pl.phenotype  <- plot(conditional_effects(model.phenotype, c("age.diff", "weight.diff", "space_overlap",
                                                             
                                                             "time.since.fl.ID1:age.diff",
                                                             "time.since.fl.ID1:weight.diff"
                                                             ), 
                                          int_conditions = int.conditions))

library(ggpubr)


ggarrange(

  #  age difference
  pl.phenotype$age.diff +
    labs(x= "Age difference [days]", y= "") +
    ylim(c(0.0,0.106))+
    theme_bw()+
    ylab("association strength (SRI)")+
    geom_line(colour="black", lwd=1),
  
  # weight difference
  pl.phenotype$weight.diff +
    labs(x= "Fledge weight difference [g]", y= "") +
    ylab("")+
    ylim(c(0.0,0.106))+
    theme_bw()+
    geom_line(colour="black", lwd=1),
  

  # age similarity over time
  pl.phenotype$`time.since.fl.ID1:age.diff` +
    labs(x= "Time since fledging [days]", y= "association strength (SRI)") +
    ylim(c(0.02,0.106))+
    theme_bw()+
 #   theme( legend.background=element_blank())+
    theme(legend.position = c(0.73, 0.25), legend.text = element_text(size=8))+
    scale_fill_brewer(palette="RdYlBu", breaks=c(0, 20, 40), labels=c(0, 20, 40), name = "age diff [days]")+
    scale_color_brewer(palette="RdYlBu", breaks=c(0, 20, 40), labels=c(0, 20, 40), name = "age diff [days]")+
    theme(legend.title = element_text(size=10),legend.key = element_rect(fill = "white")),
  
# age similarity over time
pl.phenotype$`time.since.fl.ID1:weight.diff` +
  labs(x= "Time since fledging [days]", y= "") +
  labs(y= "") +
  ylim(c(0.02,0.106))+
  theme_bw()+
#  theme( legend.background=element_blank())+
  theme(legend.position = c(0.73, 0.25), legend.text = element_text(size=8))+
  scale_fill_brewer(palette="RdYlBu", breaks=c(1, 3, 5), labels=c(1, 3, 5), name = "weight diff [g]")+
  scale_color_brewer(palette="RdYlBu", breaks=c(1, 3, 5), labels=c(1, 3, 5), name = "weight diff [g]")+
  theme(legend.title = element_text(size=10),legend.key = element_rect(fill = "white")),
 
  labels= c("a", "b", "c", "d")
  

)

ggsave("Output/Figures/Figure_Analysis2.tiff", units="in", width=6, height=6, dpi=300, compression = 'lzw')



# 4 ) Inheritance of parental associations -----------------------------
# Analysis 3
# load the data frame (same as for model 2)

# have to reload, as we subset it to only juveniles in the previous analysis
load("Data/df.summer.all.week.fledgies.RData")

df.summer.all.week.inh <- subset(df.summer.all.week.fledgies, !is.na(df.summer.all.week.fledgies$sum.assoc.parents))


# 4.1. Direct inheritance (association with adults) -----------------------

# subset the data frame to only juvenile/adult associations
df.summer.all.week.inh.ad <- subset(df.summer.all.week.fledgies, df.summer.all.week.fledgies$type=="adult")

length(df.summer.all.week.inh.ad$assoc)
# 7415 dyads
length(unique(df.summer.all.week.inh.ad$ID1))
# 65 juvneiles
length(unique(df.summer.all.week.inh.ad$ID2))
# 60 adults

# check for multi-collinearity
model.vif <- lm( assoc ~ space_overlap + sum.assoc.parents + time.since.fl.ID1 , data = df.summer.all.week.inh.ad)
library(car)
vif(model.vif)
# space_overlap sum.assoc.parents time.since.fl.ID1 
# 1.287721          1.285966          1.001651 



# Run Bayesian linear regression
library(brms)
library(rstan)


# we test whether associations between adult and parents predict association between adult and juvenile, while controlling for space use

# Model 3a
model.inh.ad <-
  brms::brm(
    assoc ~  sum.assoc.parents + 
      space_overlap +
      sum.assoc.parents:scale(time.since.fl.ID1) +
      (1 |mm(ID1, ID2)) + (1 |ID1) ,
    df.summer.all.week.inh.ad,
    family = zero_inflated_beta(),
#    prior=   c(prior(normal(0, 2), class= Intercept),
#               prior(normal(0,  5), class= b)),
    control = list(adapt_delta = .99), # set higher because of divergent transitions
    chains=4,
    iter=6000,
    cores = 6
  )

#save(model.inh.ad, file="Output/Analysis 3/model.fl.formation.inh.ad.RDA")
load("Output/Analysis 3/model.fl.formation.inh.ad.RDA")


plot(model.inh.ad)

# and posterior predictive checks
pp_check(model.inh.ad, ndraws= 1e2)

summary(model.inh.ad)
# Family: zero_inflated_beta 
# Links: mu = logit; phi = identity; zi = identity 
# Formula: assoc ~ sum.assoc.with.parents + space_overlap + sum.assoc.with.parents:scale(time.since.fl.ID1) + (1 | mm(ID1, ID2)) + (1 | ID1) 
# Data: df.summer.all.week.fledgies.inh.adults (Number of observations: 7415) 
# Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 12000
# 
# Group-Level Effects: 
#   ~ID1 (Number of levels: 65) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.09      0.06     0.00     0.22 1.00      784     2273
# 
# ~mmID1ID2 (Number of levels: 125) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.44      0.05     0.34     0.55 1.00     2216     4148
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# Intercept                                        -2.59      0.06    -2.71    -2.47 1.00     5014
# sum.assoc.with.parents                           -0.03      0.20    -0.43     0.38 1.00    11686
# space_overlap                                     0.69      0.05     0.60     0.78 1.00    14294
# sum.assoc.with.parents:scaletime.since.fl.ID1     0.83      0.10     0.64     1.02 1.00    20797
# Tail_ESS
# Intercept                                         7169
# sum.assoc.with.parents                            9768
# space_overlap                                    10047
# sum.assoc.with.parents:scaletime.since.fl.ID1     9035
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# phi    18.21      0.47    17.32    19.14 1.00    19804     8458
# zi      0.55      0.01     0.54     0.56 1.00    24907     8474
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


exp(fixef(model.inh.ad))
# Estimate Est.Error       Q2.5     Q97.5
# Intercept                                     0.07534555  1.060789 0.06683017 0.0842368
# sum.assoc.with.parents                        0.97440694  1.226645 0.64913232 1.4573657
# space_overlap                                 2.00131563  1.047380 1.82878498 2.1894189
# sum.assoc.with.parents:scaletime.since.fl.ID1 2.28738414  1.101939 1.89040212 2.7611296


# create a plot
# setting the conditions for plotting
int.conditions <- list(
  time.since.fl.ID1 = setNames(c(10,30,60), c(10,30,60)),
  sum.assoc.with.parents = setNames(c(0,0.25, 0.5), c(0,0.25, 0.5))
)


pl.inh.ad <- plot(conditional_effects(model.inh.ad, c("time.since.fl.ID1:sum.assoc.with.parents"), int_conditions = int.conditions))



# 4.2.) Inheritance of social networks among juveniles ----------------------------------------------------
# Model 3b

# subset the data frame to only juvenile/adult associations
df.summer.all.week.inh.fl <- subset(df.summer.all.week.fledgies, df.summer.all.week.fledgies$type=="peers")

length(df.summer.all.week.inh.fl$assoc)
# 9420 dyads
length(unique(df.summer.all.week.inh.fl$ID1))
# 65 juvneiles
length(unique(df.summer.all.week.inh.fl$ID2))
# 65 juveniles



library(car)

model.vif <- lm( assoc ~ space_overlap + sum.assoc.parents + time.since.fl.ID1, data = df.summer.all.week.inh.fl)

vif(model.vif)

# space_overlap sum.assoc.parents time.since.fl.ID1 
# 1.072248          1.054027          1.023844 


library(brms)
library(rstan)


# we test whether associaitons among parents predict association juveniles, while controlling for space use

# Model 3b
model.inh.juv <-
  brms::brm(
    assoc ~  sum.assoc.among.parents + 
      space_overlap +
      sum.assoc.among.parents:scale(time.since.fl.ID1) +
      (1 |mm(ID1, ID2)) + (1 |ID1) ,
    df.summer.all.week.fledgies.inh,
    family = zero_inflated_beta(),
    #    prior=   c(prior(normal(0, 2), class= Intercept),
    #               prior(normal(0,  5), class= b)),
   # control = list(adapt_delta = .99), # set higher because of divergent transitions
    chains=4,
    iter=6000,
    cores = 6
  )

save(model.inh.juv, file="Output/Analysis 3/model.fl.formation.inh.juv.RDA")
load("Output/Analysis 3/model.fl.formation.inh.juv.RDA")


plot(model.inh.juv)

# and posterior predictive checks
pp_check(model.inh.juv, ndraws= 1e2)

summary(model.inh.juv)
# Family: zero_inflated_beta 
# Links: mu = logit; phi = identity; zi = identity 
# Formula: assoc ~ sum.assoc.among.parents + space_overlap + sum.assoc.among.parents:scale(time.since.fl.ID1) + (1 | mm(ID1, ID2)) + (1 | ID1) 
# Data: df.summer.all.week.fledgies.inh (Number of observations: 8906) 
# Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 12000
# 
# Group-Level Effects: 
#   ~ID1 (Number of levels: 65) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.14      0.03     0.09     0.20 1.00     3780     6392
# 
# ~mmID1ID2 (Number of levels: 65) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.63      0.07     0.50     0.78 1.00     4220     7055
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
# Intercept                                         -2.53      0.09    -2.71    -2.36 1.00     2367
# sum.assoc.among.parents                            0.12      0.11    -0.09     0.33 1.00    14619
# space_overlap                                      0.79      0.05     0.70     0.88 1.00     7924
# sum.assoc.among.parents:scaletime.since.fl.ID1     0.21      0.06     0.09     0.33 1.00    16854
# Tail_ESS
# Intercept                                          4223
# sum.assoc.among.parents                            9426
# space_overlap                                      8545
# sum.assoc.among.parents:scaletime.since.fl.ID1     9615
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# phi    12.95      0.25    12.46    13.44 1.00    17849     9148
# zi      0.39      0.01     0.38     0.40 1.00    20773     8632
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

exp(fixef(model.inh.juv))
# Estimate Est.Error       Q2.5      Q97.5
# Intercept                                      0.07960202  1.093582 0.06659405 0.09466547
# sum.assoc.among.parents                        1.12539572  1.111518 0.91337654 1.38869146
# space_overlap                                  2.21281347  1.046505 2.02000332 2.41483046
# sum.assoc.among.parents:scaletime.since.fl.ID1 1.23756493  1.062250 1.09851113 1.39044023

# create a plot
# setting the conditions for plotting
int.conditions <- list(
  time.since.fl.ID1 = setNames(c(10,30,60), c(10,30,60)),
  sum.assoc.among.parents = setNames(c(0,0.25, 0.5), c(0,0.25, 0.5))
)


pl.inh.juv <- plot(conditional_effects(model.inh.juv, c("time.since.fl.ID1:sum.assoc.among.parents"), int_conditions = int.conditions))

# 5) Analysis 4: Assortment by age across seasons -------------------------------------------------


# 5.1) Prepare network data set --------------------------------------------

gmm.summer.3weeks <- NULL
gmm.summer.3weeks$gbi <- rbind(gmm.summer.w12$gbi,
                               gmm.summer.w13$gbi,
                               gmm.summer.w14$gbi)

gmm.summer.3weeks$metadata <- rbind(gmm.summer.w12$metadata,
                                    gmm.summer.w13$metadata,
                                    gmm.summer.w14$metadata)


create.net.gen <- function(gmm.data, season){
  object <- NULL
  age.vec <- NULL
  # create network
  IDs.10sightings <- names(which(colSums(gmm.data$gbi)>=5))
  # subset to GTs
  IDs.10sightings <- IDs.10sightings[IDs.10sightings %in% GT.list]
  
  net <- get_network(
    association_data = gmm.data$gbi,
    data_format = "GBI",
    association_index = "SRI",
    times = gmm.data$metadata$Start,
    identities = colnames(gmm.data$gbi),
    which_identities = IDs.10sightings
  )
  
  IDs <- rownames(net)
  
  # extract the ages of the all IDs
   for(i in IDs){
     if(i %in% GT.juveniles){
       age.ID <- "fledgling"  
     } else {
       age.ID <- "adult"
     }
     age.vec <- c(age.vec, age.ID)
     
   }
  
    
  object$net <- net
  object$age <- age.vec
  object$season <- season
  
  return(object)
}


summer.net <- create.net.gen(gmm.data = gmm.summer.3weeks, season = "summer")
autumn.net <- create.net.gen(gmm.data = gmm.autumn, season = "autumn")
winter.net <- create.net.gen(gmm.data = gmm.winter, season = "winter")
spring.net <- create.net.gen(gmm.data = gmm.spring, season = "spring")

# extract the total number of birds in each season
dim(summer.net$net)
dim(autumn.net$net)
dim(winter.net$net)
dim(spring.net$net)

# how many first-year and adult birds
length(summer.net$age[summer.net$age=="fledgling"])
length(summer.net$age[summer.net$age=="adult"])

length(autumn.net$age[autumn.net$age=="fledgling"])
length(autumn.net$age[autumn.net$age=="adult"])

length(winter.net$age[winter.net$age=="fledgling"])
length(winter.net$age[winter.net$age=="adult"])

length(spring.net$age[spring.net$age=="fledgling"])
length(spring.net$age[spring.net$age=="adult"])

# 5.2) Calculate assortment (analysis 4) --------------------------------------------


library(assortnet)


# calculate significance by doing node-permutation

set.seed(5)

assortment.function <- function(network){
    vec.rand <- NULL
    assort <- assortment.discrete(graph=network$net, types = network$age, weighted = TRUE)
    object <- NULL
  for(i in 1:1000){
    # we use node based permutation
    rand.phenotype <- sample(network$age)
    r.rand <- assortment.discrete(graph = network$net, types = rand.phenotype)$r
    vec.rand[i] <- r.rand
    
  }
    # extract where the real r falls among the computed r (which corresponds to our p value)
  p <- length(which(vec.rand > assort$r))/1000
  object$p <- p
  object$r <- assort$r
  return(object)
}

# summer (Model 4a)
assortment.function(network = summer.net)

# $p
# [1] 0.099
# $r
# [1] 0.00558

# autumn (Model 4b)
set.seed(10)
assortment.function(network = autumn.net)
# $p
# [1] 0.032
# 
# $r
# [1] 0.09698163

# winter (Model 4c)
set.seed(8)
assortment.function(network=winter.net)
# $p
# [1] 0.233
# 
# $r
# [1] -0.005050664

# spring (Model 4d)
set.seed(6)
assortment.function(network=spring.net)
# $p
# [1] 0.267
# 
# $r
# [1] -0.01198448


# 6) Analysis 5: Network stability across seasons -------------------------------------
gmm.summer.3weeks <- NULL
gmm.summer.3weeks$gbi <- rbind(gmm.summer.w12$gbi,
                               gmm.summer.w13$gbi,
                               gmm.summer.w14$gbi)

gmm.summer.3weeks$metadata <- rbind(gmm.summer.w12$metadata,
                                    gmm.summer.w13$metadata,
                                    gmm.summer.w14$metadata)

# 6.1) Prepare dyadic data ------------------------------------------------


create.df.gen <- function(gmm.data, season, prev.net=NULL){
  
  # create network
  IDs.10sightings <- names(which(colSums(gmm.data$gbi)>=5))
  # subset to GTs
  IDs.10sightings <- IDs.10sightings[IDs.10sightings %in% GT.list]
  
  net <- get_network(
    association_data = gmm.data$gbi,
    data_format = "GBI",
    association_index = "SRI",
    times = gmm.data$metadata$Start,
    identities = colnames(gmm.data$gbi),
    which_identities = IDs.10sightings
  )
  
  IDs <- rownames(net)
  comb.IDs <- combn(IDs,2, simplify=TRUE)
  
  # we again double up so that we can repeat the analyses from the perspective of the fledging and the adults
  df.net <- rbind.data.frame(cbind.data.frame(comb.IDs[1,], comb.IDs[2,]),
                             cbind.data.frame(comb.IDs[2,], comb.IDs[1,]))

  colnames(df.net) <- c("ID1", "ID2")
  
  df.net.comb <- df.net
  df.net.comb$season <- season
  
  
  # we here calculate the space use for each individual (i.e. how many times they have used which feeder)
  space_use <- as.data.frame(matrix(0, ncol=6, nrow=length(IDs.10sightings)))
  rownames(space_use) <- IDs.10sightings
  colnames(space_use) <-   c("Mill1", "Mill2", "Mill3", "Mill4", "Mill5", "Mill6")
  for(i in IDs.10sightings){
    # extract how many times they were seen at which feeder
    tab.ID1 <- table(gmm.data$metadata$Location[which(gmm.data$gbi[,which(colnames(gmm.data$gbi)%in% i)]==1)])
    for(j in names(tab.ID1)){
      space_use[i,j] <- tab.ID1[j]
    }
    
  }
  # this table is needed further down
  
  for(i in 1:length(df.net.comb[,1])){
    ID1 <- df.net.comb[i, "ID1"]
    ID2 <- df.net.comb[i, "ID2"]
    # extract the association strength of ID1 and ID2
    
    ID1.ID2.together <- net[ID1, ID2]
    
    df.net.comb[i, "assoc"] <- ID1.ID2.together
    
    if(ID1 %in% GT.juveniles){
      age.ID1 <- "fledgling"  
    } else {
      age.ID1 <- "adult"
    }
    
    if(ID2 %in% GT.juveniles){
      age.ID2 <- "fledgling"  
    } else {
      age.ID2 <- "adult"
    }
    
    df.net.comb[i, "age.ID1"] <- age.ID1
    df.net.comb[i, "age.ID2"] <- age.ID2    
    
    if(age.ID1==age.ID2){
      df.net.comb[i, "same.age"] <- 1
      
    } else {
      df.net.comb[i, "same.age"] <- 0
    }
      
    

    # if there is a previous season, we extract the dyad's previous network strength 
    if(!(is.null(prev.net))){
      # we don't know if ID1 and ID2 were the same in the prev network, so we have to consider both ways
      v1 <- intersect(which(prev.net$ID2==ID1), which(prev.net$ID1==ID2))
      v2 <- intersect(which(prev.net$ID1==ID1), which(prev.net$ID2==ID2))
      prev.assoc <- prev.net[unique(c(v1, v2)[1]), "assoc"]
      if(length(prev.assoc)==0){
        prev.assoc <- 0
      }
      df.net.comb[i, "prev.assoc"] <- prev.assoc
      
    }
    
    # finally, we extract the space use similarity between individuals 1 and 2
    space.use.ID1 <- space_use[rownames(space_use) %in% c(ID1),]
    space.use.ID2 <- space_use[rownames(space_use) %in% c(ID2),]
    
    time.spent.min <- NULL
    
    for(k in 1:6){
      # at each feeder location, extract how much time they have each spent there
      feederk.ID1 <- space.use.ID1[k]
      feederk.ID2 <- space.use.ID2[k]
      
      time.spent.min[k] <- min(unlist(c(feederk.ID1, feederk.ID2)))
    }
    
    # divide by the total number of visits by individual ID1
    space.overlap <- sum(time.spent.min) / sum(space.use.ID1) 
    df.net.comb[i, "space_overlap"] <- space.overlap
    
  }

  df.net.comb$same.age <- as.factor(df.net.comb$same.age)
  
  return(df.net.comb)
}

# run function - for autumn, winter and spring, we always include the association during the previous season as a predictor
df.fl.summer.gen <- create.df.gen(gmm.data=gmm.summer.3weeks, season="summer")
df.fl.autumn.gen <- create.df.gen(gmm.data=gmm.autumn, season="autumn", prev.net = df.fl.summer.gen)
df.fl.winter.gen <- create.df.gen(gmm.data=gmm.winter, season="winter", prev.net = df.fl.autumn.gen)
df.fl.spring.gen <- create.df.gen(gmm.data=gmm.spring, season="spring", prev.net = df.fl.winter.gen)


# extract the number of dyads, first-year and adults in each season
# for dyads with data across consecutive seasons

length(df.fl.autumn.gen[!is.na(df.fl.autumn.gen$prev.assoc),"ID1"])
length(unique(subset(df.fl.autumn.gen$ID1, df.fl.autumn.gen$ID1 %in% GT.juveniles & !is.na(df.fl.autumn.gen$prev.assoc))))
length(unique(subset(df.fl.autumn.gen$ID1, df.fl.autumn.gen$ID1 %in% GT.adults& !is.na(df.fl.autumn.gen$prev.assoc))))

length(df.fl.winter.gen[!is.na(df.fl.winter.gen$prev.assoc),"ID1"])
length(unique(subset(df.fl.winter.gen$ID1, df.fl.winter.gen$ID1 %in% GT.juveniles & !is.na(df.fl.winter.gen$prev.assoc))))
length(unique(subset(df.fl.winter.gen$ID1, df.fl.winter.gen$ID1 %in% GT.adults& !is.na(df.fl.winter.gen$prev.assoc))))

length(df.fl.spring.gen[!is.na(df.fl.spring.gen$prev.assoc),"ID1"])
length(unique(subset(df.fl.spring.gen$ID1, df.fl.spring.gen$ID1 %in% GT.juveniles & !is.na(df.fl.spring.gen$prev.assoc))))
length(unique(subset(df.fl.spring.gen$ID1, df.fl.spring.gen$ID1 %in% GT.adults& !is.na(df.fl.spring.gen$prev.assoc))))


# 6.2. Run brm models -----------------------------------------------------

library(brms)
library(rstan)

# now we run models 9-11 for previous associations predicting association strength
# this reduces the data set to dyads that were observed over two consecutive seasons

# Model 5a
model.autumn.prev <-
  brms::brm(
    assoc~  prev.assoc*age.ID1 + space_overlap + (1|ID1) + (1 |mm(ID1, ID2)),
    df.fl.autumn.gen,
    family = zero_inflated_beta(),
    chains=4,
    iter=6000,
    cores = 6
  )
# save the object
#save(model.autumn.prev, file="Output/Analysis 5/model.autumn.prev.RDA")
load("Output/Analysis 5/model.autumn.prev.RDA")


# Model 5b
model.winter.prev <-
  brms::brm(
    assoc~  prev.assoc*age.ID1+ space_overlap + (1|ID1) + (1 |mm(ID1, ID2)),
    df.fl.winter.gen,
    family = zero_inflated_beta(),
    chains=4,
    iter=6000,
    cores = 6
  )
# save the object
#save(model.winter.prev, file="Output/Analysis 5/model.winter.prev.RDA")
load("Output/Analysis 5/model.winter.prev.RDA")

# Model 5c
model.spring.prev <-
  brms::brm(
    assoc~ prev.assoc*age.ID1+ space_overlap + (1|ID1) + (1 |mm(ID1, ID2)),
    df.fl.spring.gen,
    family = zero_inflated_beta(),
    chains=4,
    iter=6000,
    cores = 6
  )
# save the object
#save(model.spring.prev, file="Output/Analysis 5/model.spring.prev.RDA")
load("Output/Analysis 5/model.spring.prev.RDA")

# 6.3) Model checks -------------------------------------------------------

plot(model.autumn.prev) # Model 5a
plot(model.winter.prev) # Model 5b
plot(model.spring.prev) # Model 5c

# chains seems to have mixed well

# and posterior predictive checks
pp_check(model.autumn.prev, ndraws= 1e2) # Model 5a
pp_check(model.winter.prev, ndraws= 1e2) # Model 5b
pp_check(model.spring.prev, ndraws= 1e2) # Model 5c


# 6.5) Look at effects ----------------------------------------------------

# Model 5a
summary(model.autumn.prev)
# Family: zero_inflated_beta 
# Links: mu = logit; phi = identity; zi = identity 
# Formula: assoc ~ prev.assoc * age.ID1 + space_overlap + (1 | ID1) + (1 | mm(ID1, ID2)) 
# Data: df.fl.autumn.gen (Number of observations: 870) 
# Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 12000
# 
# Group-Level Effects: 
#   ~ID1 (Number of levels: 29) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.39      0.11     0.21     0.63 1.00     3424     4598
# 
# ~mmID1ID2 (Number of levels: 30) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.49      0.12     0.29     0.74 1.00     3847     6026
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                      -2.80      0.21    -3.23    -2.40 1.00     5700     7290
# prev.assoc                      2.52      0.77     0.98     4.00 1.00    10578     9473
# age.ID1fledgling                0.16      0.23    -0.30     0.60 1.00     6831     7814
# space_overlap                   0.99      0.19     0.63     1.37 1.00     5486     8275
# prev.assoc:age.ID1fledgling    -2.88      1.23    -5.29    -0.45 1.00     9546     9596
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# phi    31.20      2.94    25.72    37.20 1.00     7586     8457
# zi      0.65      0.02     0.61     0.68 1.00    26336     8336
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# get odds for model 5a (autumn)
exp(fixef(model.autumn.prev))
# Estimate Est.Error        Q2.5       Q97.5
# Intercept                    0.06072894  1.233013 0.039592936  0.09092318
# prev.assoc                  12.37673902  2.170491 2.669310902 54.64587211
# age.ID1fledgling             1.17529685  1.259784 0.741557716  1.82372704
# space_overlap                2.69661174  1.206447 1.877215405  3.91599022
# prev.assoc:age.ID1fledgling  0.05638129  3.421654 0.005024984  0.63524070

# extract probability
plogis(fixef(model.autumn.prev))

#                               Estimate Est.Error       Q2.5      Q97.5
# Intercept                   0.05725208 0.5521746 0.03808504 0.08334517
# prev.assoc                  0.92524336 0.6845914 0.72746926 0.98202922
# age.ID1fledgling            0.54029263 0.5574799 0.42580140 0.64585812
# space_overlap               0.72948200 0.5467827 0.65244173 0.79658218
# prev.assoc:age.ID1fledgling 0.05337210 0.7738403 0.00499986 0.38846923

# Model 5b
 summary(model.winter.prev)

 # Family: zero_inflated_beta 
 # Links: mu = logit; phi = identity; zi = identity 
 # Formula: assoc ~ prev.assoc * age.ID1 + space_overlap + (1 | ID1) + (1 | mm(ID1, ID2)) 
 # Data: df.fl.winter.gen (Number of observations: 1260) 
 # Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 1;
 # total post-warmup draws = 12000
 # 
 # Group-Level Effects: 
 #   ~ID1 (Number of levels: 35) 
 # Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
 # sd(Intercept)     0.65      0.10     0.47     0.88 1.00     3499     6459
 # 
 # ~mmID1ID2 (Number of levels: 36) 
 # Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
 # sd(Intercept)     0.63      0.11     0.45     0.86 1.00     3513     5959
 # 
 # Population-Level Effects: 
 #   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
 # Intercept                      -3.41      0.20    -3.80    -3.01 1.00     3257     5437
 # prev.assoc                      1.12      0.28     0.58     1.68 1.00    14021    10191
 # age.ID1fledgling               -0.08      0.24    -0.55     0.40 1.00     3053     4898
 # space_overlap                   2.64      0.09     2.46     2.83 1.00     9786     9521
 # prev.assoc:age.ID1fledgling    -0.12      0.42    -0.96     0.69 1.00    14251     9514
 # 
 # Family Specific Parameters: 
 #   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
 # phi    54.59      3.13    48.67    60.89 1.00    13018     8649
 # zi      0.45      0.01     0.42     0.47 1.00    29359     8194
 # 
 # Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
 # and Tail_ESS are effective sample size measures, and Rhat is the potential
 # scale reduction factor on split chains (at convergence, Rhat = 1).
 
 # get odds for model 5b (winter)
 exp(fixef(model.winter.prev))
 
 # Estimate Est.Error        Q2.5       Q97.5
 # Intercept                    0.03320661  1.223236  0.02228763  0.04944905
 # prev.assoc                   3.06545536  1.324828  1.78405495  5.33945693
 # age.ID1fledgling             0.92108834  1.272434  0.57555686  1.49177993
 # space_overlap               14.04674428  1.098203 11.70471122 16.89923217
 # prev.assoc:age.ID1fledgling  0.88358770  1.519478  0.38370531  2.00128572

plogis(fixef(model.winter.prev))
# Estimate Est.Error       Q2.5      Q97.5
# Intercept                   0.03213937 0.5502053 0.02180173 0.04711906
# prev.assoc                  0.75402509 0.5698607 0.64081169 0.84225778
# age.ID1fledgling            0.47946173 0.5599433 0.36530377 0.59868045
# space_overlap               0.93354044 0.5234016 0.92128904 0.94413168
# prev.assoc:age.ID1fledgling 0.46909825 0.6030924 0.27730277 0.66680946
 

# Model 5c
summary(model.spring.prev)

# Family: zero_inflated_beta 
# Links: mu = logit; phi = identity; zi = identity 
# Formula: assoc ~ prev.assoc * age.ID1 + space_overlap + (1 | ID1) + (1 | mm(ID1, ID2)) 
# Data: df.fl.spring.gen (Number of observations: 2256) 
# Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 12000
# 
# Group-Level Effects: 
#   ~ID1 (Number of levels: 47) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.93      0.12     0.73     1.18 1.00     2374     4419
# 
# ~mmID1ID2 (Number of levels: 48) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.75      0.10     0.57     0.98 1.00     2791     4282
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                      -3.49      0.24    -3.95    -3.03 1.00     1905     3913
# prev.assoc                      0.81      0.26     0.31     1.33 1.00     9277     8956
# age.ID1fledgling                0.23      0.31    -0.38     0.83 1.00     2042     3645
# space_overlap                   2.16      0.09     1.98     2.34 1.00     9285     8568
# prev.assoc:age.ID1fledgling     0.77      0.39    -0.01     1.54 1.00    11491     9458
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# phi    43.67      2.05    39.70    47.73 1.00    12027     8695
# zi      0.54      0.01     0.52     0.56 1.00    15228     8679
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# get odds for model 5c (spring)
exp(fixef(model.spring.prev))

#                            Estimate Est.Error       Q2.5       Q97.5
# Intercept                   0.030459  1.265735 0.01916581  0.04839989
# prev.assoc                  2.257571  1.293650 1.36770143  3.76482721
# age.ID1fledgling            1.264373  1.357279 0.68614033  2.30329053
# space_overlap               8.669631  1.097243 7.22849492 10.34912345
# prev.assoc:age.ID1fledgling 2.152064  1.481931 0.98977716  4.64999069


plogis(fixef(model.spring.prev))
# Estimate Est.Error       Q2.5      Q97.5
# Intercept                   0.02955867 0.5586421 0.01880539 0.04616549
# prev.assoc                  0.69302278 0.5640137 0.57764945 0.79012880
# age.ID1fledgling            0.55837673 0.5757821 0.40692955 0.69727156
# space_overlap               0.89658343 0.5231835 0.87847109 0.91188747
# prev.assoc:age.ID1fledgling 0.68274760 0.5970879 0.49743116 0.82300856

# 6.6) Make plot ----------------------------------------------------------

int.conditions <- list(
  prev.assoc = setNames(c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)), 
  space_overlap = setNames(c(0.5), c(0.5))
)


pl.autumn.prev <- plot(conditional_effects(model.autumn.prev, 
                                           c("prev.assoc:age.ID1",
                                             "prev.assoc", 
                                             "age.ID1"), int_conditions = int.conditions))

pl.winter.prev <- plot(conditional_effects(model.winter.prev, 
                                           c("prev.assoc:age.ID1",
                                             "prev.assoc", 
                                             "age.ID1"), int_conditions = int.conditions))

pl.spring.prev <- plot(conditional_effects(model.spring.prev, 
                                           c("prev.assoc:age.ID1",
                                             "prev.assoc", 
                                             "age.ID1"), int_conditions = int.conditions))


library(ggpubr)

mycol2 <- c("#b23714", "#00177e")


 ggarrange(
  
  # assoc previous season

  pl.autumn.prev$`prev.assoc:age.ID1`+
    labs(x= "association prev. season", y= "association strength (SRI)")+
    ylim(c(0,0.14))+
    theme_bw()+
    theme(legend.position='none')+
    scale_colour_manual(name = "age class", values=mycol2, labels=c("adult", "first year"))+
    scale_fill_manual(name = "age class", values=mycol2, labels=c("adult", "first year"))+
    ggtitle("autumn")+ theme(plot.title = element_text(hjust = 0.5)),
  
  pl.winter.prev$`prev.assoc:age.ID1`+
    labs(x= "association prev. season", y= "")+
    ylim(c(0,0.14))+
    theme_bw()+
    theme(legend.position='none')+
    scale_colour_manual(name = "age class", values=mycol2, labels=c("adult", "first year"))+
    scale_fill_manual(name = "age class", values=mycol2, labels=c("adult", "first year"))+
    ggtitle("winter")+ theme(plot.title = element_text(hjust = 0.5)),
  
  pl.spring.prev$`prev.assoc:age.ID1`+
    labs(x= "association prev. season", y= "")+
    ylim(c(0,0.14))+
    theme_bw()+
    theme(legend.position='none')+
    scale_colour_manual(name = "age class", values=mycol2, labels=c("adult", "first year"))+
    scale_fill_manual(name = "age class", values=mycol2, labels=c("adult", "first year"))+
    ggtitle("spring")+ theme(plot.title = element_text(hjust = 0.5)),
  
  ncol=3,
  nrow=1,
 # widths = c(0.2, 0.2, 0.2, 0.3),
  common.legend = TRUE,
 legend = "right"
# labels= c("b", "c", "d")
)

# save as tiff
 ggsave("Output/Figures/Figure_Analysis4.tiff", units="in", width=8, height=3, dpi=300, compression = 'lzw')
 


# 7) Plot networks --------------------------------------------------------


#install.packages("GGally")
library(GGally)
library(RColorBrewer)

#install.packages("network")
library(network)
#install.packages("sna")
library(sna)
library(ggplot2)
library(igraph)
 library(qgraph)

plot.network <- function(gmm.data, season){
  IDs.to.include <- names(which(colSums(gmm.data$gbi)>=5))
  # subset to gretis
  IDs.to.include <- intersect(IDs.to.include, GT.list)
  
  net <- get_network(
    association_data = gmm.data$gbi,
    data_format = "GBI",
    association_index = "SRI",
    times = gmm.data$metadata$Start,
    identities = colnames(gmm.data$gbi),
    which_identities = IDs.to.include
  )
  
  IDs <- rownames(net)

  
  # for better visibility, we set the edge weights under 0.01 to 0
  
  net[net<0.05] <- 0
  
  age <- NULL
  for(i in rownames(net)){
    age.i <- unique(subset(species.age$Age_in_2020, species.age$Pit==i))
    age[which(rownames(net)==i)] <- age.i
  }
  
  
  
  shape <- age
  g.net <- graph_from_adjacency_matrix(net, mode = "undirected",
                                       weighted = TRUE, diag = TRUE)
  V(g.net)$age <- age
  
  E(g.net)$width <- E(g.net)$weight
 
  age[age=="juvenile"] <- "#00177e" # grey
  age[age=="adult"] <- "#b23714" # white
  
  shape[shape=="juvenile"] <- "circle" # grey
  shape[shape=="adult"] <- "square" # white
  
  E(g.net)$width <- E(g.net)$weight

  V(g.net)$colour <- age
  
  V(g.net)$shape <- shape

  j <- which(V(g.net)$age == "juvenile")
  a <- which(V(g.net)$age == "adult")

  
  E(g.net)$color <- "#8c8989"
  E(g.net)[j %--% j]$color <- "#00177e"
  E(g.net)[a %--% a]$color <- "#b23714"

  e <- get.edgelist(g.net,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g.net),
                                         area=8*(vcount(g.net)^2),repulse.rad=(vcount(g.net)^2.9))
  

  
 p <-  igraph::plot.igraph( g.net,
                       vertex.size = 6,
                       edge.curved = 0.25,
                       edge.color = adjustcolor(E(g.net)$color, alpha = 0.3) ,
                       vertex.color = V(g.net)$colour,
                       vertex.shape = V(g.net)$shape,
                       vertex.label = NA,
                       vertex.frame.colour = "black",
                       edge.width = E(g.net)$width*20,
                       frame = FALSE,
                       layout=l,
                       asp = 1,
                       #     mark.groups = list.D,
                       #    mark.border =NA,
                       #  mark.col=col.adj
                       edge.width = E(g.net)$weight,
                       main=season
                       
                       
  )
  
  

return(p)
  

}


tiff("Output/Figures/Figure_Networks_summer.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
plot.summer <- plot.network(gmm.data =gmm.summer.3weeks, season="Summer")
dev.off()
tiff("Output/Figures/Figure_Networks_autumn.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
plot.autumn <- plot.network(gmm.data =gmm.autumn, season="Autumn")
dev.off()
tiff("Output/Figures/Figure_Networks_winter.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
plot.winter <- plot.network(gmm.data =gmm.winter, season="Winter")
dev.off()
tiff("Output/Figures/Figure_Networks_spring.tiff", units="in", width=5, height=5, res=300, compression = 'lzw')
plot.spring <- plot.network(gmm.data =gmm.spring, season="Spring")
dev.off()

