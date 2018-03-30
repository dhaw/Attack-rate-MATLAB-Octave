<<<<<<< HEAD
#Actually subtype, not strain
library(reshape2)
library(ggplot2)
library(plyr)

prepStrain <- function(df, 
                       dfID, 
                       zones, 
                       proportion=FALSE){ #df=dfNet
assign("proportion", proportion, .GlobalEnv)
nh <- zones$ISO3[zones$hemis=="northern"]
sh <- zones$ISO3[zones$hemis=="southern"]
#tr <- zones$ISO3[zones$hemis=="tropics"]

dfx <- data.frame(df$ISO3, df$ISO_YEAR, df$ISO_WEEK, df$VIRUS_SUB_CODE, df$ValueNumeric)
colnames(dfx) <- c("code", "year", "week", "subtype", "value") #"Year" of season start

################################

if (proportion==FALSE){
  dfid <- dfID[which(dfID$MEASURE_CODE=="ILI_CASES" & dfID$AGEGROUP_CODE=="All"), ]
  dfid <- data.frame(dfid$ISO3, dfid$ISO_YEAR, dfid$ISO_WEEK, dfid$ValueNumeric)
  colnames(dfid) <- c("code", "year", "week", "ili")
  #Combines with dfx here:
  dfx <- merge(dfx, dfid, by=c("code", "year", "week"))
  dfx$ili[is.na(dfx$value)] <- 0
  dfx$ili[is.na(dfx$ili)] <- 0
  dfx$value <- dfx$value*dfx$ili
  dfx$ili <- NULL

#ILI time series additional output:
dfn <- dfid[dfid$code %in% nh, ]
dfs <- dfid[dfid$code %in% sh, ]
lateSeason <- which(as.numeric(dfn$week) <= 39)
dfn$year[lateSeason] <- dfn$year[lateSeason]-1
dfid <- rbind.data.frame(dfn, dfs)
dfid$week <- NULL

dfid <- ddply(dfid, c("code", "year"), function(dfid)sum(dfid$ili))
colnames(dfid)[3] <- "ili"
dfid <- reshape(dfid, idvar="code", timevar=c("year"), direction="wide")
write.csv(dfid, "timeSerisILI.csv")
}


################################

dfn <- dfx[dfx$code %in% nh, ]
dfs <- dfx[dfx$code %in% sh, ]
#dft <- dfx[dfx$code %in% tr, ]

lateSeason <- which(as.numeric(dfn$week) <= 39)
dfn$year[lateSeason] <- dfn$year[lateSeason]-1

subtypes <- levels(dfn$subtype)
subtypes <- subtypes[!subtypes %in% c("", "TOTAL")]

dfn <- ddply(dfn, c("code", "year", "subtype"), function(dfn)sum(dfn$value))
dfs <- ddply(dfs, c("code", "year", "subtype"), function(dfs)sum(dfs$value))
#dft <- ddply(dft, c("code", "year", "subtype"), function(dft)sum(dft$value))

dfout <- rbind.data.frame(dfn, dfs)#, dft)
dfout <- reshape(dfout, idvar=c("code", "year"), timevar=c("subtype"), direction="wide")
dfout <- dfout[, c("code", "year", "V1.AH1", "V1.AH1N12009", "V1.AH3", "V1.AH5", 
               "V1.ANOTSUBTYPED", "V1.AOTHER_SUBTYPE", 
               "V1.BNOTDETERMINED", "V1.BVICTORIA", "V1.BYAMAGATA", 
               "V1.AH7N9")]
#To date, 2x H1N1, 0x non-subtypable/other subtype ********

write.csv(dfout, "strainYear.csv")
}

################################

extractCountry <- function(dfall, 
                           iso3="USA", 
                           splitOther=FALSE)
  
country <- dfall[dfall$code==iso3, ]
rownames(country) <- country$year
country <- country[, c(seq(3,12))]#(4,12),15)]
colnames(country) <- c("H1", "H1_2009", "H3", "H5", "A_NA", "A_other", "B_NA", "B_Vic", "B_Yam", "H7")#DH: matlab version excludes H7
country[is.na(country)] <- 0
#Add H5 & H7 to A other:
country$A_other <- country$A_other+country$H5+country$H7
#Add pandemic to H1:
country$H1 <- country$H1+country$H1_2009
country <- country[, c("H1", "H3", "A_other", "A_NA", "B_Vic", "B_Yam", "B_NA")]
#country$H1, country$H3, country$A_other, country$A_NA, country$B_Vic, country$B_Yam, country$B_other)]

if (splitOther==TRUE){
  country$H1 <- country$H1+.5*country$A_other
  country$H3 <-country$H3+.5*country$A_other
  country <- data.frame(country$H1, country$H3, country$A_NA, country$B_Vic, country$B_Yam, country$B_other) #Just delete column!
  palette <- c("#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c")
} else {
  palette <- c("#a6cee3",
               "#1f78b4",
               "#b2df8a",
               "#33a02c",
               "#fb9a99",
               "#e31a1c",
               "#fdbf6f")
}
nVirus <- ncol(country)

dat <- data.frame(t(country))
dat[is.na(dat)] <- 0
if (proportion==TRUE){
  dat <- sweep(dat, 2, colSums(dat), FUN="/")
}
dat$row <- seq_len(nrow(dat))
dat2 <- melt(dat, id.vars = "row")

ggplot(dat2, aes(x=variable, y=value, fill=as.factor(row))) + 
  geom_bar(stat="identity") +
  xlab("\nYear of season start") +
  ylab("Proportion of samples\n") +
  scale_fill_manual(values=palette)
}
=======
#Actually subtype, not strain
library(reshape2)
library(ggplot2)
library(plyr)

prepStrain <- function(df, 
                       dfID, 
                       zones, 
                       proportion=FALSE){ #df=dfNet
assign("proportion", proportion, .GlobalEnv)
nh <- zones$ISO3[zones$hemis=="northern"]
sh <- zones$ISO3[zones$hemis=="southern"]
#tr <- zones$ISO3[zones$hemis=="tropics"]

dfx <- data.frame(df$ISO3, df$ISO_YEAR, df$ISO_WEEK, df$VIRUS_SUB_CODE, df$ValueNumeric)
colnames(dfx) <- c("code", "year", "week", "subtype", "value") #"Year" of season start

################################

if (proportion==FALSE){
  dfid <- dfID[which(dfID$MEASURE_CODE=="ILI_CASES" & dfID$AGEGROUP_CODE=="All"), ]
  dfid <- data.frame(dfid$ISO3, dfid$ISO_YEAR, dfid$ISO_WEEK, dfid$ValueNumeric)
  colnames(dfid) <- c("code", "year", "week", "ili")
  #Combines with dfx here:
  dfx <- merge(dfx, dfid, by=c("code", "year", "week"))
  dfx$ili[is.na(dfx$value)] <- 0
  dfx$ili[is.na(dfx$ili)] <- 0
  dfx$value <- dfx$value*dfx$ili
  dfx$ili <- NULL

#ILI time series additional output:
dfn <- dfid[dfid$code %in% nh, ]
dfs <- dfid[dfid$code %in% sh, ]
lateSeason <- which(as.numeric(dfn$week) <= 39)
dfn$year[lateSeason] <- dfn$year[lateSeason]-1
dfid <- rbind.data.frame(dfn, dfs)
dfid$week <- NULL

dfid <- ddply(dfid, c("code", "year"), function(dfid)sum(dfid$ili))
colnames(dfid)[3] <- "ili"
dfid <- reshape(dfid, idvar="code", timevar=c("year"), direction="wide")
write.csv(dfid, "timeSerisILI.csv")
}


################################

dfn <- dfx[dfx$code %in% nh, ]
dfs <- dfx[dfx$code %in% sh, ]
#dft <- dfx[dfx$code %in% tr, ]

lateSeason <- which(as.numeric(dfn$week) <= 39)
dfn$year[lateSeason] <- dfn$year[lateSeason]-1

subtypes <- levels(dfn$subtype)
subtypes <- subtypes[!subtypes %in% c("", "TOTAL")]

dfn <- ddply(dfn, c("code", "year", "subtype"), function(dfn)sum(dfn$value))
dfs <- ddply(dfs, c("code", "year", "subtype"), function(dfs)sum(dfs$value))
#dft <- ddply(dft, c("code", "year", "subtype"), function(dft)sum(dft$value))

dfout <- rbind.data.frame(dfn, dfs)#, dft)
dfout <- reshape(dfout, idvar=c("code", "year"), timevar=c("subtype"), direction="wide")
dfout <- dfout[, c("code", "year", "V1.AH1", "V1.AH1N12009", "V1.AH3", "V1.AH5", 
               "V1.ANOTSUBTYPED", "V1.AOTHER_SUBTYPE", 
               "V1.BNOTDETERMINED", "V1.BVICTORIA", "V1.BYAMAGATA", 
               "V1.AH7N9")]
#To date, 2x H1N1, 0x non-subtypable/other subtype ********

write.csv(dfout, "strainYear.csv")
}

################################

extractCountry <- function(dfall, 
                           iso3="USA", 
                           splitOther=FALSE)
  
country <- dfall[dfall$code==iso3, ]
rownames(country) <- country$year
country <- country[, c(seq(3,12))]#(4,12),15)]
colnames(country) <- c("H1", "H1_2009", "H3", "H5", "A_NA", "A_other", "B_NA", "B_Vic", "B_Yam", "H7")#DH: matlab version excludes H7
country[is.na(country)] <- 0
#Add H5 & H7 to A other:
country$A_other <- country$A_other+country$H5+country$H7
#Add pandemic to H1:
country$H1 <- country$H1+country$H1_2009
country <- country[, c("H1", "H3", "A_other", "A_NA", "B_Vic", "B_Yam", "B_NA")]
#country$H1, country$H3, country$A_other, country$A_NA, country$B_Vic, country$B_Yam, country$B_other)]

if (splitOther==TRUE){
  country$H1 <- country$H1+.5*country$A_other
  country$H3 <-country$H3+.5*country$A_other
  country <- data.frame(country$H1, country$H3, country$A_NA, country$B_Vic, country$B_Yam, country$B_other) #Just delete column!
  palette <- c("#a6cee3",
    "#1f78b4",
    "#b2df8a",
    "#33a02c",
    "#fb9a99",
    "#e31a1c")
} else {
  palette <- c("#a6cee3",
               "#1f78b4",
               "#b2df8a",
               "#33a02c",
               "#fb9a99",
               "#e31a1c",
               "#fdbf6f")
}
nVirus <- ncol(country)

dat <- data.frame(t(country))
dat[is.na(dat)] <- 0
if (proportion==TRUE){
  dat <- sweep(dat, 2, colSums(dat), FUN="/")
}
dat$row <- seq_len(nrow(dat))
dat2 <- melt(dat, id.vars = "row")

ggplot(dat2, aes(x=variable, y=value, fill=as.factor(row))) + 
  geom_bar(stat="identity") +
  xlab("\nYear of season start") +
  ylab("Proportion of samples\n") +
  scale_fill_manual(values=palette)
}
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
