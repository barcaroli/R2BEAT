## ----setup, include=FALSE-----------------------------------------------------------------------------------------
options(width=400)
knitr::opts_chunk$set(echo = TRUE, fig.width=6, fig.height=4)
# Install ReGenesees if not already installed
if (requireNamespace("ReGenesees", quietly = TRUE)) {
  svystat    <- ReGenesees::svystat
} else {
  stop("The package ReGenesees is needed. \nInstall it by executing the following: \ndevtools::install_github('DiegoZardetto/ReGenesees'")
}
library(ReGenesees)
library(R2BEAT)
library(plyr)
library(sampling)
options(warn=-1)
options(scipen=9999)


## -----------------------------------------------------------------------------------------------------------------
load("R2BEAT_ReGenesees.RData")   # ReGenesees design object


## -----------------------------------------------------------------------------------------------------------------
des


## -----------------------------------------------------------------------------------------------------------------
cal


## -----------------------------------------------------------------------------------------------------------------
# Control the presence of strata with less than two units
ls <- find.lon.strata(des)


## -----------------------------------------------------------------------------------------------------------------
str(des$variables)


## -----------------------------------------------------------------------------------------------------------------
summary(des$variables$income_hh)


## -----------------------------------------------------------------------------------------------------------------
table(des$variables$work)


## -----------------------------------------------------------------------------------------------------------------
table(des$variables$unemployed)


## ---- eval=FALSE--------------------------------------------------------------------------------------------------
## des<-des.addvars(des,active=factor(ifelse(work==1,1,0)))
## des<-des.addvars(des,inactive=factor(ifelse(work==2,1,0)))
## cal<-des.addvars(cal,active=factor(ifelse(work==1,1,0)))
## cal<-des.addvars(cal,inactive=factor(ifelse(work==2,1,0)))


## -----------------------------------------------------------------------------------------------------------------
table(cal$variables$active)


## -----------------------------------------------------------------------------------------------------------------
table(cal$variables$inactive)



## -----------------------------------------------------------------------------------------------------------------
table(cal$variables$unemployed)


## -----------------------------------------------------------------------------------------------------------------
RGdes <- des                           # ReGenesees design object
RGcal <- cal                           # ReGenesees calibrated object

strata_vars <- c("stratum")            # variables of stratification
target_vars <- c("income_hh",
                 "active",
                 "inactive",
                 "unemployed")         # target variables
deff_vars <- "stratum"                 # stratification variables to be used when calculating deff and effst 
                                       #    (n.b: must coincide or be a subset of variables of stratification)
id_PSU <- c("municipality")            # identification variable of PSUs
id_SSU <- c("id_hh")                   # identification variable of SSUs
domain_vars <- c("region")             # domain variables
inp1 <- input_to_beat.2st_1(RGdes,
                            RGcal,
                            id_PSU,
                            id_SSU,
                            strata_vars,
                            target_vars,
                            deff_vars,
                            domain_vars)


## -----------------------------------------------------------------------------------------------------------------
head(inp1$strata)


## -----------------------------------------------------------------------------------------------------------------
head(inp1$deff)


## -----------------------------------------------------------------------------------------------------------------
head(inp1$effst)


## -----------------------------------------------------------------------------------------------------------------
head(inp1$rho)


## -----------------------------------------------------------------------------------------------------------------
# psu <- read.csv2("psu.csv") # Read the external file containing PSU information
head(psu)
psu_id="municipality"        # Identifier of the PSU
stratum_var="stratum"        # Identifier of the stratum
mos_var="ind"                # Variable to be used as 'measure of size'
delta=1                      # Average number of SSUs for each selection unit
minimum <- 50                # Minimum number of SSUs to be selected in each PSU
inp2 <- input_to_beat.2st_2(psu,
                            psu_id,
                            stratum_var,
                            mos_var,
                            delta,
                            minimum)
head(inp2$psu_file)
head(inp2$des_file)


## -----------------------------------------------------------------------------------------------------------------
cv <- as.data.frame(list(DOM=c("DOM1","DOM2"),
                         CV1=c(0.03,0.04),
                         CV2=c(0.06,0.08),
                         CV3=c(0.06,0.08),
                         CV4=c(0.06,0.08)))
cv


## -----------------------------------------------------------------------------------------------------------------
stratif = inp1$strata 
errors = cv 
des_file = inp2$des_file 
psu_file = inp2$psu_file 
rho = inp1$rho 
effst = inp1$effst

alloc <- beat.2st(stratif, 
                  errors, 
                  des_file, 
                  psu_file, 
                  rho, 
                  deft_start = NULL, 
                  effst,
                  epsilon1 = 5, 
                  mmdiff_deft = 1,maxi = 15, 
                  epsilon = 10^(-11), minnumstrat = 2, maxiter = 200, maxiter1 = 25)


## -----------------------------------------------------------------------------------------------------------------
alloc$sensitivity


## -----------------------------------------------------------------------------------------------------------------
alloc$expected


## -----------------------------------------------------------------------------------------------------------------
allocat <- alloc$alloc[-nrow(alloc$alloc),]
sample_2st <- StratSel(dataPop= inp2$psu_file,
                       idpsu= ~ PSU_ID, 
                       dom= ~ STRATUM, 
                       final_pop= ~ PSU_MOS, 
                       size= ~ PSU_MOS, 
                       PSUsamplestratum= 1, 
                       min_sample= minimum, 
                       min_sample_index= FALSE, 
                       dataAll=allocat,
                       domAll= ~ factor(STRATUM), 
                       f_sample= ~ ALLOC, 
                       planned_min_sample= NULL, 
                       launch= F)


## -----------------------------------------------------------------------------------------------------------------
sample_2st[[2]]


## ---- echo=TRUE---------------------------------------------------------------------------------------------------
des <- sample_2st[[2]]
des <- des[1:(nrow(des)-1),]
strat <- c(as.character(as.numeric(des$Domain[1:(nrow(des)-1)])),"Tot")
barplot(t(des[1:(nrow(des)),2:3]), names=strat,
        col=c("darkblue","red"), las=2, xlab = "Stratum", cex.axis=0.7, cex.names=0.7)
legend("topleft", 
       legend = c("Self Representative","Non Self Representative"),
       fill = c("darkblue", "red"))
title("Distribution of allocated PSUs by domain")



## ---- echo=TRUE---------------------------------------------------------------------------------------------------
barplot(t(des[1:(nrow(des)),5:6]), names=strat,
        col=c("darkblue","red"), las=2, xlab = "Stratum", cex.axis=0.7, 
        cex.names=0.7)
legend("topleft", 
       legend = c("Self Representative","Non Self Representative"),
       fill = c("darkblue", "red"))
title("Distribution of allocated SSUs by domain")


## -----------------------------------------------------------------------------------------------------------------
selected_PSU <- sample_2st[[4]]
selected_PSU <- selected_PSU[selected_PSU$PSU_final_sample_unit > 0,]
write.table(sample_2st[[4]],"Selected_PSUs.csv",sep=";",row.names=F,quote=F)
head(selected_PSU)

