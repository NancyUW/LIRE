################## Libraries ########################
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(plyr)
library(scales)
library(Hmisc)
library(data.table)
library(zoo)
library(icd9)

################## Lire Data  ########################
# There are differences between the ICD9 codes in this data set and in the 
# quan publication. Why?

lire = read.csv("~/Documents/LIRE/lire_idx_2015_09_08.csv")
lire$indexDate = as.Date(as.character(lire$indexDate), '%m/%d/%y')
lire = lire[lire$indexDate >= '2013-10-01' & lire$indexDate <= '2015-03-31',]


lire_co = read.csv("~/Documents/LIRE/lire_comorbidity_test_2015_09_17.csv")
names(lire_co) = c("visitID", "isInpatientVisit", "visitCode", "visitDate",
                   "isPrimaryDx", "icd9")
lire_co$visitID = as.factor(lire_co$visitID)

dim(lire_co)
dim(lire)

lire_co = lire_co[which(lire_co$visitID %in% lire$studyCode),]

################## Dates ###################################
lire_co$visitDate = as.Date(as.character(lire_co$visitDate), '%m/%d/%y')

################## ICD-9 IsReal and Isvalid #################
# These are patients with unreal ICD-9 codes, ie
# the icd9 does not correspond to a real ICD-9 code
lire_co$real = icd9IsReal(lire_co$icd9)

# These are the patients with syntactically invalid ICD-9 codes. 
lire_co$valid = icd9IsValid(lire_co$icd9, isShort = TRUE) 

lire_co$explain = ifelse(lire_co$real == TRUE, icd9Explain(lire_co$icd9), NA)

################## Comorbidity MAPPINGS ######################
ahrq_Comorbid = icd9Comorbid(lire_co, ahrqComorbid)
elix_Comorbid = icd9Comorbid(lire_co, elixComorbid)
quandeyo_Comorbid = icd9Comorbid(lire_co, quanDeyoComorbid, return.df = TRUE)
quanelix_Comorbid = icd9Comorbid(lire_co, quanElixComorbid)

################## Comorbidity INDICES #######################
charlson = icd9Charlson(lire_co, isShort = FALSE, return.df = TRUE)

################## Joining Charlson and Quan/Deyo  ###########
lire_quandeyo = join(lire_co, quandeyo_Comorbid)
lire_charlson = join(lire_quandeyo, charlson)

################## Exploratory Plotting ######################
# ggplot() + 
#   geom_point(data = lire_charlson, aes(x = visitID, y = Charlson)) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   xlab("Visit ID")

this = melt(quandeyo_Comorbid, id = "visitID")
this = this[this$value == 1,]
this$value = as.numeric(this$value)
this$variable_ordered = with(this, reorder(variable, variable, function(x) -length(x)))
ggplot() + 
  geom_bar(data = this, aes(x = variable_ordered), stat = "bin") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################## QuanDeyo Info #####################################
# The Quon paper is MORE inclusive than the R package and the codes that Chris (from email)
# generated. The goal is to figure out which codes are NOT in the R/Chris, and 
# potentially do a literature review of *why* there are differences. 

# First: Verify that the R ICD9 codes for quanDeyoComorbid$_______ are the same as the
# codes from Chris. If these are the same, cross-check with the Quon codes. 

# names(quanDeyoComorbid)
# "MI"          "CHF"         "PVD"         "Stroke"      "Dementia"    
# "Pulmonary"   "Rheumatic"   "PUD"         "LiverMild"   "DM"          
# "DMcx"        "Paralysis"   "Renal"       "Cancer"      "LiverSevere" 
# "Mets"        "HIV"  

################## Chris's ICD9 codes #######################
# Read in a long data frame of comorbidity mappings from Chris. 
df_chris = read.csv("~/Documents/LIRE/chris_comorbidities.csv")
df_chris$code = icd9ShortToDecimal(df_chris$code)

# Wrap the data frame into a list
chris_list = split(df_chris$code, list(df_chris$comorbidity))

################## Exhaustive ICD9s #############################
# This is an exhaustive list of ICD9 codes
all = read.csv("/Users/organn/Documents/LIRE/icd9_dx_codes.csv", colClasses = rep("character", 3))
all$dxCode = icd9ShortToDecimal(all$dxCode)

################# {icd9} quanDeyoComorbid and quanElixComorbid in decimals ###########
myquandeyo = lapply(quanDeyoComorbid, icd9ShortToDecimal)
myquanelix = lapply(quanElixComorbid, icd9ShortToDecimal)

################# Quan paper Enhanced ICD-9-CM ###################
quan_paper = read.csv("~/Documents/LIRE/quan_mappings.csv")
names(quan_paper) = c("code", "comorbidity")

quan_list = split(quan_paper$comorbidity, list(quan_paper$code))

################# Comparison Functions ##########################
compare_icd9 = function(comparison, reference){
  diff_list = list()
  for(name in names(comparison)){
    # Which codes are in comparison that are not in reference?
    diff_list[[name]] = setdiff(comparison[[name]], reference[[name]])
  }
  return(diff_list[diff_list != "character(0)"])
}

expand = function(small, large){
  expanded_list = list()
  for(i in 1:length(small)){
    # For every element in small list, grep() against large list
    expanded_list[[i]] = unique(grep(paste(
      paste("^", small[[i]], ".?", sep = ""),
      collapse="|"), large, value=TRUE))
  }
  names(expanded_list) = names(small)
  return(expanded_list)
}

################# Comparisons ##########################
chris_expanded = expand(small = chris_list, large = melt(myquandeyo)$value)
quan_expanded = expand(small = quan_list, large = melt(myquandeyo)$value)

## What codes does Chris have that are not in {icd9}?
chris_not_in_r = compare_icd9(comparison = chris_expanded, reference = myquandeyo)

## What codes does Chris have that aren't in the Quan paper?
chris_not_in_quan = compare_icd9(comparison = chris_expanded, reference = quan_expanded)

## What codes are in {icd9} that are not in Chris's codes?
r_not_in_chris = compare_icd9(comparison = myquandeyo, reference =  chris_expanded)

## What codes are in {icd9} that are not in the Quan paper?
r_not_in_quan = compare_icd9(comparison = myquandeyo, reference = quan_expanded)

## What codes are in the Quan paper that are not in {icd9}?
quan_not_in_r = compare_icd9(comparison = quan_expanded, reference = myquandeyo)

## What codes are in the Quan paper that Chris does not have?
quan_not_in_chris = compare_icd9(comparison = quan_expanded, reference = chris_expanded)


###################################################################
########### Pattern matching against R's mappings #################

# Expand out Chris's codes by pattern matching against the exhaustive list 


# Expand chris and quan with R instead of the comprehensive list
# in icd9 but not in chris
# in icd9 but not in quan paper
# 
# chris not in quan paper
# quan paper not in chris

# What version of Elixhauser are in the R package and how they map to 
# one or more of the Elixhauser mappings in the Quan paper
