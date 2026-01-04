######---------------------------------------------------------------------#####
### 01: Model fitting ###
######---------------------------------------------------------------------#####

# Creator: Abigail M. Weber
# Correspondence: amw8076@psu.edu

##################
# load libraries #
##################

library(performance)
library(AICcmodavg)
library(DHARMa)
library(gee)
library(lmtest)
library(sandwich)

################################################################################
# Step 1: Run global model for landscape perspective to test for collinearity #
################################################################################

# Read in modeling data
data <- read.csv("data/data_for_HSF_models.csv")
data <- data[ ,-1]

# Subset to be only fawn and random points
land_perspective <- subset(data, species2 == "fawn" | species2 == "random")

# Obtain complete cases across all buffers
names(land_perspective)
land_perspective <- land_perspective[complete.cases(land_perspective[ ,10:189]), ] # Removed 1662 rows

# 50m model
land_50 <- glm(case ~  scale(water_wetland_prop_50m) + scale(exurban_prop_50m) + scale(other_prop_50m) +
              scale(crop_prop_50m) + scale(edge_prop_50m) + scale(slope_50m) +
              log_road_dist_50m + log_water_dist_50m + scale(elev_50m) + 
              scale(canopy_height_50m) + scale(TRI_50m), 
           family = binomial, data = land_perspective)

check_collinearity(land_50) # Issues with slope/TRI

# 100m buffer
land_100 <- glm(case ~  scale(water_wetland_prop_100m) + scale(exurban_prop_100m) + scale(other_prop_100m) +
               scale(crop_prop_100m) + scale(edge_prop_100m) + scale(slope_100m) +
               log_road_dist_100m + log_water_dist_100m + scale(elev_100m) + scale(canopy_height_100m) +
               scale(TRI_100m), 
            family = binomial, data = land_perspective)

check_collinearity(land_100) # Slope and TRI

# 150m buffer
land_150 <- glm(case ~  scale(water_wetland_prop_150m) + scale(exurban_prop_150m) +
               scale(other_prop_150m) +scale(crop_prop_150m) + scale(edge_prop_150m) +
               scale(slope_150m) + log_road_dist_150m + log_water_dist_150m + scale(elev_150m) +
               scale(canopy_height_150m) + scale(TRI_150m), 
            family = binomial, data = land_perspective)

check_collinearity(land_150) # Slope and TRI

# 200m buffer
land_200 <- glm(case ~  scale(water_wetland_prop_200m) + scale(exurban_prop_200m) +
               scale(other_prop_200m) +scale(crop_prop_200m) + scale(edge_prop_200m) +
               scale(slope_200m) + log_road_dist_200m + log_water_dist_200m + scale(elev_200m) +
               scale(canopy_height_200m) +
               scale(TRI_200m), 
             family = binomial, data = land_perspective)

check_collinearity(land_200) # Slope and TRI

# 250m buffer
land_250 <- glm(case ~  scale(water_wetland_prop_250m) + scale(exurban_prop_250m) +
               scale(other_prop_250m) +scale(crop_prop_250m) + scale(edge_prop_250m) +
               scale(slope_250m) + log_road_dist_250m + log_water_dist_250m + scale(elev_250m) + 
               scale(canopy_height_250m) + scale(TRI_250m), 
             family = binomial, data = land_perspective)

check_collinearity(land_250) # Slope and TRI

# 300m buffer
land_300 <- glm(case ~  scale(water_wetland_prop_300m) + scale(exurban_prop_300m) +
               scale(other_prop_300m) +scale(crop_prop_300m) + scale(edge_prop_300m) +
               scale(slope_300m) + log_road_dist_300m + log_water_dist_300m + scale(elev_300m) + 
               scale(canopy_height_300m) + scale(TRI_300m), 
             family = binomial, data = land_perspective)

check_collinearity(land_300) # Slope and TRI

# 350m buffer
land_350 <- glm(case ~  scale(water_wetland_prop_350m) + scale(exurban_prop_350m) +
               scale(other_prop_350m) +scale(crop_prop_350m) + scale(edge_prop_350m) +
               scale(slope_350m) + log_road_dist_350m + log_water_dist_350m + scale(elev_350m) +
               scale(canopy_height_350m) + scale(TRI_350m), 
             family = binomial, data = land_perspective)

check_collinearity(land_350) # Slope and TRI

# 400m buffer
land_400 <- glm(case ~  scale(water_wetland_prop_400m) + scale(exurban_prop_400m) +
              scale(other_prop_400m) +scale(crop_prop_400m) + scale(edge_prop_400m) +
              scale(slope_400m) + log_road_dist_400m + log_water_dist_400m + scale(elev_400m) +
              scale(canopy_height_400m) + scale(TRI_400m), 
             family = binomial, data = land_perspective)

check_collinearity(land_400) # Slope and TRI

# 450m buffer
land_450 <- glm(case ~  scale(water_wetland_prop_450m) + scale(exurban_prop_450m) +
               scale(other_prop_450m) +scale(crop_prop_450m) + scale(edge_prop_450m) +
               scale(slope_450m) + log_road_dist_450m + log_water_dist_450m + scale(elev_450m) +
               scale(canopy_height_450m) + scale(TRI_450m), 
             family = binomial, data = land_perspective)

check_collinearity(land_450) # Slope and TRI

# 500m buffer
land_500 <- glm(case ~  scale(water_wetland_prop_500m) + scale(exurban_prop_500m) +
               scale(other_prop_500m) +scale(crop_prop_500m) + scale(edge_prop_500m) +
               scale(slope_500m) + log_road_dist_500m + log_water_dist_500m + scale(elev_500m) +
               scale(canopy_height_500m) + scale(TRI_500m), 
             family = binomial, data = land_perspective)

check_collinearity(land_500) # Slope, TRI, and now edge (VIF = 5.18)
# Keep edge for now since correlation doesn't appear until broadest scale
  # Reassess if 500m is the best buffer size

# RESULT: Remove slope from all subsequent analyses

################################################################################
# Step 2: Determine most supported buffer size for the landscape perspective #
################################################################################

# Remove everything from environment except 'data'
remove(land_perspective, land_50, land_100, land_150, land_200, land_250,
       land_300, land_350, land_400, land_450, land_500)

# Subset data for the landscape perspective: fawn and random locations
land_perspective <- subset(data, species2 == "fawn" | species2 == "random")

# Obtain complete cases across all buffers
names(land_perspective)
land_perspective <- land_perspective[complete.cases(land_perspective[ ,10:189]), ] # Removed 1662 rows

# Check that the case is right for each model comparison: fawn = 1 while random = 0
table(land_perspective$species2, land_perspective$case) # fawn = 1, random = 0

# Conduct models across buffer size
# 50m model
land_50 <- glm(case ~  scale(water_wetland_prop_50m) + scale(exurban_prop_50m) +
              scale(other_prop_50m) + scale(crop_prop_50m) + scale(edge_prop_50m) +
              log_road_dist_50m + log_water_dist_50m + scale(elev_50m)+
              scale(canopy_height_50m) + scale(TRI_50m),
            family = binomial, data = land_perspective)

# 100m buffer
land_100 <- glm(case ~  scale(water_wetland_prop_100m) + scale(exurban_prop_100m) +
               scale(other_prop_100m) + scale(crop_prop_100m) + scale(edge_prop_100m) +
               log_road_dist_100m + log_water_dist_100m + scale(elev_100m)+
               scale(canopy_height_100m) + scale(TRI_100m),
             family = binomial, data = land_perspective)

# 150m buffer
land_150 <- glm(case ~  scale(water_wetland_prop_150m) + scale(exurban_prop_150m) +
               scale(other_prop_150m) + scale(crop_prop_150m) + scale(edge_prop_150m) +
               log_road_dist_150m + log_water_dist_150m + scale(elev_150m)+
               scale(canopy_height_150m) + scale(TRI_150m),
             family = binomial, data = land_perspective)

# 200m buffer
land_200 <- glm(case ~  scale(water_wetland_prop_200m) + scale(exurban_prop_200m) +
               scale(other_prop_200m) + scale(crop_prop_200m) + scale(edge_prop_200m) +
               log_road_dist_200m + log_water_dist_200m + scale(elev_200m)+
               scale(canopy_height_200m) + scale(TRI_200m),
             family = binomial, data = land_perspective)

# 250m buffer
land_250 <- glm(case ~  scale(water_wetland_prop_250m) + scale(exurban_prop_250m) +
               scale(other_prop_250m) + scale(crop_prop_250m) + scale(edge_prop_250m) +
               log_road_dist_250m + log_water_dist_250m + scale(elev_250m)+
               scale(canopy_height_250m) + scale(TRI_250m),
             family = binomial, data = land_perspective)

# 300m buffer
land_300 <- glm(case ~  scale(water_wetland_prop_300m) + scale(exurban_prop_300m) +
               scale(other_prop_300m) + scale(crop_prop_300m) + scale(edge_prop_300m) +
               log_road_dist_300m + log_water_dist_300m + scale(elev_300m)+
               scale(canopy_height_300m) + scale(TRI_300m),
             family = binomial, data = land_perspective)

# 350m buffer
land_350 <- glm(case ~  scale(water_wetland_prop_350m) + scale(exurban_prop_350m) +
               scale(other_prop_350m) + scale(crop_prop_350m) + scale(edge_prop_350m) +
               log_road_dist_350m + log_water_dist_350m + scale(elev_350m)+
               scale(canopy_height_350m) + scale(TRI_350m),
             family = binomial, data = land_perspective)

# 400m buffer
land_400 <- glm(case ~  scale(water_wetland_prop_400m) + scale(exurban_prop_400m) +
               scale(other_prop_400m) + scale(crop_prop_400m) + scale(edge_prop_400m) +
               log_road_dist_400m + log_water_dist_400m + scale(elev_400m)+
               scale(canopy_height_400m) + scale(TRI_400m),
             family = binomial, data = land_perspective)

# 450m buffer
land_450 <- glm(case ~  scale(water_wetland_prop_450m) + scale(exurban_prop_450m) +
               scale(other_prop_450m) + scale(crop_prop_450m) + scale(edge_prop_450m) +
               log_road_dist_450m + log_water_dist_450m + scale(elev_450m)+
               scale(canopy_height_450m) + scale(TRI_450m),
             family = binomial, data = land_perspective)

# 500m buffer
land_500 <- glm(case ~  scale(water_wetland_prop_500m) + scale(exurban_prop_500m) +
               scale(other_prop_500m) + scale(crop_prop_500m) + scale(edge_prop_500m) +
               log_road_dist_500m + log_water_dist_500m + scale(elev_500m)+
               scale(canopy_height_500m) + scale(TRI_500m),
             family = binomial, data = land_perspective)

# AIC score across landawn RSland models
Cand_mods1 <- list("land_50" = land_50,
                   "land_100" = land_100,
                   "land_150" = land_150,
                   "land_200" = land_200,
                   "land_250" = land_250,
                   "land_300" = land_300,
                   "land_350" = land_350,
                   "land_400" = land_400,
                   "land_450" = land_450,
                   "land_500" = land_500)

# Compute AIC table
AICc_table_land <- aictab(cand.set = Cand_mods1)
AICc_table_land

################################################################################
# Step 3: Determine most supported buffer size for the prey perspective  #
################################################################################

# Subset data for the landscape perspective: fawn and doe locations
prey_perspective <- subset(data, species2 == "fawn" | species2 == "doe")
# fp <- subset(data, species2 == "fawn" | species2 == "predator")

# Obtain complete cases across all buffers
names(prey_perspective)
prey_perspective <- prey_perspective[complete.cases(prey_perspective[ ,10:189]), ] # Removed 0 rows
# fp_comp <- fp[complete.cases(fp[ ,14:193]), ] # Removed 396 rows

# Check that the case is right for each model comparison: fawn = 1 while doe = 0
table(prey_perspective$species2, prey_perspective$case) #fawn = 1, doe = 0
# table(fp$species2, fp$case) #fawn = 1, predator = 0

# Conduct models across buffer size
# 50m model
prey_50 <- glm(case ~  scale(water_wetland_prop_50m) + scale(exurban_prop_50m) +
              scale(other_prop_50m) + scale(crop_prop_50m) + scale(edge_prop_50m) +
              log_road_dist_50m + log_water_dist_50m + scale(elev_50m)+
              scale(canopy_height_50m) + scale(TRI_50m),
            family = binomial, data = prey_perspective)

# 100m buffer
prey_100 <- glm(case ~  scale(water_wetland_prop_100m) + scale(exurban_prop_100m) +
               scale(other_prop_100m) + scale(crop_prop_100m) + scale(edge_prop_100m) +
               log_road_dist_100m + log_water_dist_100m + scale(elev_100m)+
               scale(canopy_height_100m) + scale(TRI_100m),
             family = binomial, data = prey_perspective)

# 150m buffer
prey_150 <- glm(case ~  scale(water_wetland_prop_150m) + scale(exurban_prop_150m) +
               scale(other_prop_150m) + scale(crop_prop_150m) + scale(edge_prop_150m) +
               log_road_dist_150m + log_water_dist_150m + scale(elev_150m)+
               scale(canopy_height_150m) + scale(TRI_150m),
             family = binomial, data = prey_perspective)

# 200m buffer
prey_200 <- glm(case ~  scale(water_wetland_prop_200m) + scale(exurban_prop_200m) +
               scale(other_prop_200m) + scale(crop_prop_200m) + scale(edge_prop_200m) +
               log_road_dist_200m + log_water_dist_200m + scale(elev_200m)+
               scale(canopy_height_200m) + scale(TRI_200m),
             family = binomial, data = prey_perspective)

# 250m buffer
prey_250 <- glm(case ~  scale(water_wetland_prop_250m) + scale(exurban_prop_250m) +
               scale(other_prop_250m) + scale(crop_prop_250m) + scale(edge_prop_250m) +
               log_road_dist_250m + log_water_dist_250m + scale(elev_250m)+
               scale(canopy_height_250m) + scale(TRI_250m),
             family = binomial, data = prey_perspective)

# 300m buffer
prey_300 <- glm(case ~  scale(water_wetland_prop_300m) + scale(exurban_prop_300m) +
               scale(other_prop_300m) + scale(crop_prop_300m) + scale(edge_prop_300m) +
               log_road_dist_300m + log_water_dist_300m + scale(elev_300m)+
               scale(canopy_height_300m) + scale(TRI_300m),
             family = binomial, data = prey_perspective)

# 350m buffer
prey_350 <- glm(case ~  scale(water_wetland_prop_350m) + scale(exurban_prop_350m) +
               scale(other_prop_350m) + scale(crop_prop_350m) + scale(edge_prop_350m) +
               log_road_dist_350m + log_water_dist_350m + scale(elev_350m)+
               scale(canopy_height_350m) + scale(TRI_350m),
             family = binomial, data = prey_perspective)

# 400m buffer
prey_400 <- glm(case ~  scale(water_wetland_prop_400m) + scale(exurban_prop_400m) +
               scale(other_prop_400m) + scale(crop_prop_400m) + scale(edge_prop_400m) +
               log_road_dist_400m + log_water_dist_400m + scale(elev_400m)+
               scale(canopy_height_400m) + scale(TRI_400m),
             family = binomial, data = prey_perspective)

# 450m buffer
prey_450 <- glm(case ~  scale(water_wetland_prop_450m) + scale(exurban_prop_450m) +
               scale(other_prop_450m) + scale(crop_prop_450m) + scale(edge_prop_450m) +
               log_road_dist_450m + log_water_dist_450m + scale(elev_450m)+
               scale(canopy_height_450m) + scale(TRI_450m),
             family = binomial, data = prey_perspective)

# 500m buffer
prey_500 <- glm(case ~  scale(water_wetland_prop_500m) + scale(exurban_prop_500m) +
               scale(other_prop_500m) + scale(crop_prop_500m) + scale(edge_prop_500m) +
               log_road_dist_500m + log_water_dist_500m + scale(elev_500m)+
               scale(canopy_height_500m) + scale(TRI_500m),
             family = binomial, data = prey_perspective)

# get the AIC score across fawn-doe LSD models
Cand_mods2 <- list("prey_50" = prey_50,
                   "prey_100" = prey_100,
                   "prey_150" = prey_150,
                   "prey_200" = prey_200,
                   "prey_250" = prey_250,
                   "prey_300" = prey_300,
                   "prey_350" = prey_350,
                   "prey_400" = prey_400,
                   "prey_450" = prey_450,
                   "prey_500" = prey_500)

# Compute AIC table
AICc_table_prey <- aictab(cand.set = Cand_mods2)
AICc_table_prey

################################################################################
# Step 4: Determine most supported buffer size for the predator perspective  #
################################################################################

# Subset data for the landscape perspective: fawn and doe locations
pred_perspective <- subset(data, species2 == "fawn" | species2 == "predator")

# Obtain complete cases across all buffers
names(pred_perspective)
pred_perspective <- pred_perspective[complete.cases(pred_perspective[ ,10:189]), ] # Removed 396 rows

# Check that the case is right for each model comparison: fawn = 1 while pred = 0
table(pred_perspective$species2, pred_perspective$case) #fawn = 1, doe = 0

# Conduct models across buffer size
# 50m model
pred_50 <- glm(case ~  scale(water_wetland_prop_50m) + scale(exurban_prop_50m) +
               scale(other_prop_50m) + scale(crop_prop_50m) + scale(edge_prop_50m) +
               log_road_dist_50m + log_water_dist_50m + scale(elev_50m)+
               scale(canopy_height_50m) + scale(TRI_50m),
             family = binomial, data = pred_perspective)

# 100m buffer
pred_100 <- glm(case ~  scale(water_wetland_prop_100m) + scale(exurban_prop_100m) +
                scale(other_prop_100m) + scale(crop_prop_100m) + scale(edge_prop_100m) +
                log_road_dist_100m + log_water_dist_100m + scale(elev_100m)+
                scale(canopy_height_100m) + scale(TRI_100m),
              family = binomial, data = pred_perspective)

# 150m buffer
pred_150 <- glm(case ~  scale(water_wetland_prop_150m) + scale(exurban_prop_150m) +
                scale(other_prop_150m) + scale(crop_prop_150m) + scale(edge_prop_150m) +
                log_road_dist_150m + log_water_dist_150m + scale(elev_150m)+
                scale(canopy_height_150m) + scale(TRI_150m),
              family = binomial, data = pred_perspective)

# 200m buffer
pred_200 <- glm(case ~  scale(water_wetland_prop_200m) + scale(exurban_prop_200m) +
                scale(other_prop_200m) + scale(crop_prop_200m) + scale(edge_prop_200m) +
                log_road_dist_200m + log_water_dist_200m + scale(elev_200m)+
                scale(canopy_height_200m) + scale(TRI_200m),
              family = binomial, data = pred_perspective)

# 250m buffer
pred_250 <- glm(case ~  scale(water_wetland_prop_250m) + scale(exurban_prop_250m) +
                scale(other_prop_250m) + scale(crop_prop_250m) + scale(edge_prop_250m) +
                log_road_dist_250m + log_water_dist_250m + scale(elev_250m)+
                scale(canopy_height_250m) + scale(TRI_250m),
              family = binomial, data = pred_perspective)

# 300m buffer
pred_300 <- glm(case ~  scale(water_wetland_prop_300m) + scale(exurban_prop_300m) +
                scale(other_prop_300m) + scale(crop_prop_300m) + scale(edge_prop_300m) +
                log_road_dist_300m + log_water_dist_300m + scale(elev_300m)+
                scale(canopy_height_300m) + scale(TRI_300m),
              family = binomial, data = pred_perspective)

# 350m buffer
pred_350 <- glm(case ~  scale(water_wetland_prop_350m) + scale(exurban_prop_350m) +
                scale(other_prop_350m) + scale(crop_prop_350m) + scale(edge_prop_350m) +
                log_road_dist_350m + log_water_dist_350m + scale(elev_350m)+
                scale(canopy_height_350m) + scale(TRI_350m),
              family = binomial, data = pred_perspective)

# 400m buffer
pred_400 <- glm(case ~  scale(water_wetland_prop_400m) + scale(exurban_prop_400m) +
                scale(other_prop_400m) + scale(crop_prop_400m) + scale(edge_prop_400m) +
                log_road_dist_400m + log_water_dist_400m + scale(elev_400m)+
                scale(canopy_height_400m) + scale(TRI_400m),
              family = binomial, data = pred_perspective)

# 450m buffer
pred_450 <- glm(case ~  scale(water_wetland_prop_450m) + scale(exurban_prop_450m) +
                scale(other_prop_450m) + scale(crop_prop_450m) + scale(edge_prop_450m) +
                log_road_dist_450m + log_water_dist_450m + scale(elev_450m)+
                scale(canopy_height_450m) + scale(TRI_450m),
              family = binomial, data = pred_perspective)

# 500m buffer
pred_500 <- glm(case ~  scale(water_wetland_prop_500m) + scale(exurban_prop_500m) +
                scale(other_prop_500m) + scale(crop_prop_500m) + scale(edge_prop_500m) +
                log_road_dist_500m + log_water_dist_500m + scale(elev_500m)+
                scale(canopy_height_500m) + scale(TRI_500m),
              family = binomial, data = pred_perspective)

# Get the AIC score across fawn-pred LSD models
Cand_mods3 <- list("pred_50" = pred_50,
                   "pred_100" = pred_100,
                   "pred_150" = pred_150,
                   "pred_200" = pred_200,
                   "pred_250" = pred_250,
                   "pred_300" = pred_300,
                   "pred_350" = pred_350,
                   "pred_400" = pred_400,
                   "pred_450" = pred_450,
                   "pred_500" = pred_500)

# Compute AIC table
AICc_table_pred <- aictab(cand.set = Cand_mods3)
AICc_table_pred

################################################################################
# Step 5: Sum log-likelihood and determine top buffer size across perspectives #
################################################################################

# View the 3 AICc tables to sum log-likelihood across
AICc_table_land
AICc_table_prey
AICc_table_pred

# Sum "LL" by column "Modnames"
names(AICc_table_land)

aic_all <- rbind(AICc_table_land, AICc_table_prey, AICc_table_pred)
aic_all
class(aic_all)

aic_df <- as.data.frame(aic_all)
class(aic_df)

# Use 'substr' on Modnames to create new matching name column
aic_df$names <- substr(aic_df$Modnames, start = 4, stop = 11)
aic_df

# Create a table which sums log-likelihood by new names
sum_ll <- aggregate(aic_df$LL ~ aic_df$names, aic_df, FUN = sum)
sum_ll

# Order sum_ll by descending log-likelihood value
sum_ll2 <- sum_ll[order(sum_ll$`aic_df$LL`), ]
max_ll <- round(abs(max(sum_ll2$`aic_df$LL`)),4)
sum_ll2$diff <- round(sum_ll2$`aic_df$LL`, 4) + max_ll
sum_ll2

# Save the sum_ll2 df
write.csv(sum_ll2, "output/sum_loglikelihood_fawn_models.csv")

# RESULT: The top buffer size across the landscape, prey, and predator perspective
  # is the 100m buffer size. However, there are many competitive buffer sizes, but the
  # coeffcients are similar across these buffer sizes. Therefore, we continue with the 100m buffer.

################################################################################
# Step 6: Run 100m model across all perspectives, including predator selection #
################################################################################

# Remove everything from the environment except "data"
remove(aic_all, aic_df, AICc_table_land, AICc_table_pred,
       AICc_table_prey, Cand_mods1, Cand_mods2, Cand_mods3,      
       land_100, land_150, land_200, land_250,       
       land_300, land_350, land_400, land_450,       
       land_50, land_500, land_perspective, max_ll,          
       pred_100, pred_150, pred_200, pred_250,       
       pred_300, pred_350, pred_400, pred_450, 
       pred_50, pred_500, pred_perspective, prey_100,        
       prey_150, prey_200, prey_250, prey_300,
       prey_350, prey_400, prey_450, prey_50,         
       prey_500, prey_perspective, sum_ll, sum_ll2)

# Subset to be only fawn and random points
land_perspective <- subset(data, species2 == "fawn" | species2 == "random")
prey_perspective <- subset(data, species2 == "fawn" | species2 == "doe")
pred_perspective <- subset(data, species2 == "fawn" | species2 == "predator")

# Obtain complete cases for 100m exclusively
names(land_perspective)

land_persp_100 <- land_perspective[complete.cases(land_perspective[ ,c(26:27,29,31,34,36,39,161,171,181)]), ] #removed 1199 points
prey_persp_100 <- prey_perspective[complete.cases(prey_perspective[ ,c(26:27,29,31,34,36,39,161,171,181)]), ] #removed 0 points
pred_persp_100 <- pred_perspective[complete.cases(pred_perspective[ ,c(26:27,29,31,34,36,39,161,171,181)]), ] #removed 221 points

# Check that the case is right for each model comparison
table(land_persp_100$species2, land_persp_100$case) # fawn = 1, random = 0
table(prey_persp_100$species2, prey_persp_100$case) # fawn = 1, doe = 0
table(pred_persp_100$species2, pred_persp_100$case) # fawn = 1, predator = 0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Landscape perspective/Fawn kill site RSF at 100m #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

land_100 <- glm(case ~  scale(water_wetland_prop_100m) + scale(exurban_prop_100m) +
                scale(other_prop_100m) + scale(crop_prop_100m) + scale(edge_prop_100m) +
                log_road_dist_100m + log_water_dist_100m + scale(elev_100m) +
                scale(canopy_height_100m) + scale(TRI_100m),
              family = binomial, data = land_persp_100)

# Check collinearity - good
check_collinearity(land_100)

# Check residuals
out <- simulateResiduals(land_100)
plot(out)

# Model summary
summary(land_100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Prey perspective/Fawn kill site-doe LSDF at 100m #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

prey_100 <- glm(case ~  scale(water_wetland_prop_100m) + scale(exurban_prop_100m) +
                scale(other_prop_100m) + scale(crop_prop_100m) + scale(edge_prop_100m) +
                log_road_dist_100m + log_water_dist_100m + scale(elev_100m) +
                scale(canopy_height_100m) + scale(TRI_100m),
              family = binomial, data = prey_persp_100)

# Check collinearity
check_collinearity(prey_100)

# Check residuals
out <- simulateResiduals(prey_100)
plot(out)

# Model summary
summary(prey_100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Predator perspective/Fawn kill site-predator LSDF at 100m #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pred_100 <- glm(case ~  scale(water_wetland_prop_100m) + scale(exurban_prop_100m) +
                scale(other_prop_100m) + scale(crop_prop_100m) + scale(edge_prop_100m) +
                log_road_dist_100m + log_water_dist_100m + scale(elev_100m) +
                scale(canopy_height_100m) + scale(TRI_100m),
              family = binomial, data = pred_persp_100)

# Check collinearity
check_collinearity(pred_100)

# Check residuals
out <- simulateResiduals(pred_100)
plot(out)

# Model summary
summary(pred_100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Predator selection/Predator RSF at 100m #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Read in predator GEE dataset
pred_gee <- read.csv("data/data_for_pred_GEE_model.csv")
pred_gee <- pred_gee[,-1]

# Check case: pred = 1, random = 0
table(pred_gee$species2, pred_gee$case)

# Subset to complete cases
names(pred_gee)
pred_select_100 <- pred_gee[complete.cases(pred_gee[ ,c(24:25,27,29,32,34,37,159,169,179)]), ] # removed 3498 rows

# Make id a factor 
pred_select_100$id <- as.factor(pred_select_100$id)

# 100m buffer model
pred_sel_100 <- gee(case ~ scale(water_wetland_prop_100m) + scale(exurban_prop_100m) +
                  scale(other_prop_100m) + scale(crop_prop_100m) + scale(edge_prop_100m) +
                  log_road_dist_100m + log_water_dist_100m + scale(elev_100m)+
                  scale(canopy_height_100m) + scale(TRI_100m),
                data = pred_select_100, family = binomial, id = id, corstr = "exchangeable")

# Check collinearity
check_collinearity(pred_sel_100)

# Model summary
summary(pred_sel_100)

################################################################################
# Step 7: Combine coefficients from each of the 4 models, use sandwich estimators #
################################################################################

# coeftest values for LSDF models
m_prey_100_coefs <- coeftest(prey_100, vcov = vcovCL, cluster = ~ id)
m_pred_100_coefs <- coeftest(pred_100, vcov = vcovCL, cluster = ~ id)

# Create a dataframe
m_prey_100_coef_df <- as.data.frame(m_prey_100_coefs[,])
m_pred_100_coef_df <- as.data.frame(m_pred_100_coefs[,])
m_land_100_coef_df <- as.data.frame(summary(land_100)$coefficients)
m_gee_100_coef_df <- as.data.frame(summary(pred_sel_100)$coefficients)

# Calculate the confidence intervals for 2.5% and 97.5% for all dfs, SE * 1.96 = CI
# Along with the 90% CIs aka 5% and 95% for all df's which is SE * 1645

# Fawn risk 100m
for(i in 1:nrow(m_land_100_coef_df)){
  m_land_100_coef_df$lower_95_CI[i] <- (m_land_100_coef_df$Estimate[i] - (m_land_100_coef_df$`Std. Error`[i] * (1.96)))
  m_land_100_coef_df$upper_95_CI[i] <- (m_land_100_coef_df$Estimate[i] + (m_land_100_coef_df$`Std. Error`[i] * (1.96)))
  
  m_land_100_coef_df$lower_90_CI[i] <- (m_land_100_coef_df$Estimate[i] - (m_land_100_coef_df$`Std. Error`[i] * (1.645)))
  m_land_100_coef_df$upper_90_CI[i] <- (m_land_100_coef_df$Estimate[i] + (m_land_100_coef_df$`Std. Error`[i] * (1.645)))
}

# Fawn doe LSDF 100m
for(i in 1:nrow(m_prey_100_coef_df)){
  m_prey_100_coef_df$lower_95_CI[i] <- (m_prey_100_coef_df$Estimate[i] - (m_prey_100_coef_df$`Std. Error`[i] * (1.96)))
  m_prey_100_coef_df$upper_95_CI[i] <- (m_prey_100_coef_df$Estimate[i] + (m_prey_100_coef_df$`Std. Error`[i] * (1.96)))
  
  m_prey_100_coef_df$lower_90_CI[i] <- (m_prey_100_coef_df$Estimate[i] - (m_prey_100_coef_df$`Std. Error`[i] * (1.645)))
  m_prey_100_coef_df$upper_90_CI[i] <- (m_prey_100_coef_df$Estimate[i] + (m_prey_100_coef_df$`Std. Error`[i] * (1.645)))
}

# Fawn pred LSDF 100m 
for(i in 1:nrow(m_pred_100_coef_df)){
  m_pred_100_coef_df$lower_95_CI[i] <- (m_pred_100_coef_df$Estimate[i] - (m_pred_100_coef_df$`Std. Error`[i] * (1.96)))
  m_pred_100_coef_df$upper_95_CI[i] <- (m_pred_100_coef_df$Estimate[i] + (m_pred_100_coef_df$`Std. Error`[i] * (1.96)))
  
  m_pred_100_coef_df$lower_90_CI[i] <- (m_pred_100_coef_df$Estimate[i] - (m_pred_100_coef_df$`Std. Error`[i] * (1.645)))
  m_pred_100_coef_df$upper_90_CI[i] <- (m_pred_100_coef_df$Estimate[i] + (m_pred_100_coef_df$`Std. Error`[i] * (1.645)))
}

# Predator GEE 100m 
for(i in 1:nrow(m_gee_100_coef_df)){
  m_gee_100_coef_df$lower_95_CI[i] <- (m_gee_100_coef_df$Estimate[i] - (m_gee_100_coef_df$`Robust S.E.`[i] * (1.96)))
  m_gee_100_coef_df$upper_95_CI[i] <- (m_gee_100_coef_df$Estimate[i] + (m_gee_100_coef_df$`Robust S.E.`[i] * (1.96)))
  
  m_gee_100_coef_df$lower_90_CI[i] <- (m_gee_100_coef_df$Estimate[i] - (m_gee_100_coef_df$`Robust S.E.`[i] * (1.645)))
  m_gee_100_coef_df$upper_90_CI[i] <- (m_gee_100_coef_df$Estimate[i] + (m_gee_100_coef_df$`Robust S.E.`[i] * (1.645)))
}

# Add columns to ID which model the coefs are coming from
m_land_100_coef_df$model <- "Fawn Kill Site (RSF)"
m_prey_100_coef_df$model <- "Fawn Kill Site - Doe Use (LSDF)"
m_pred_100_coef_df$model <- "Fawn Kill Site - Predator Use (LSDF)"
m_gee_100_coef_df$model <- "Predator Use (RSF)"

# Combine the fawn models into one all coef
all_fawn_coef <- rbind(m_land_100_coef_df, m_prey_100_coef_df, m_pred_100_coef_df)

# Turn the rownames into a column name
all2 <- all_fawn_coef
names <- rownames(all2)
rownames(all2) <- NULL
data2 <- cbind(names, all2)

# Rename the column covariate names
data2$covariate <- rep(c("Intercept",
                         "Water",
                         "Exurban",
                         "Other",
                         "Crop",
                         "Edge",
                         "Distance to Road",
                         "Distance to Water",
                         "Elevation",
                         "Canopy Height",
                         "TRI"))

# Change the names and the order of model
head(data2)
names(data2) <- c("Names", "Estimate", "SE", "Z_value", "P_value", "Lower_95_CI",
                  "Upper_95_CI", "Lower_90_CI", "Upper_90_CI", "Model", "Covariate")

# Compare names of predator gee and all_coef
names(m_gee_100_coef_df)
names(data2)

# Change names of pred gee, add/drop columns as needed
names(m_gee_100_coef_df) <- c("Estimate", "Naive_SE", "Naive_Z", "SE", "Z_value",
                              "Lower_95_CI", "Upper_95_CI", "Lower_90_CI", "Upper_90_CI", "Model")

m_gee_100_coef_df$Covariate <- rep(c("Intercept",
                                     "Water",
                                     "Exurban",
                                     "Other",
                                     "Crop",
                                     "Edge",
                                     "Distance to Road",
                                     "Distance to Water",
                                     "Elevation",
                                     "Canopy Height",
                                     "TRI"))

# Turn the rownames into a column name for gee
m_gee2 <- m_gee_100_coef_df
names <- rownames(m_gee2)
rownames(m_gee2) <- NULL
mgee <- cbind(names, m_gee2)

# Check names again
names(mgee)
names(data2)

# Drop and add columns as needed
mgee$P_value <- NA
mgee2 <- mgee[ ,c(1,2,5:6,13,7:12)]

names(mgee2) <- names(data2)

# Combine coefficient data together
coef_data <- rbind(data2, mgee2)

# Save this for figures
write.csv(coef_data, "output/all_model_coefs_100m_buffer.csv")
