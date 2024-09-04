# ----------------------------------------------------------------------------
# Fluxometrics data analysis
#
# Created: 24/12/2020
# Sarah De Beuckeleer & Hugo Steenberghen
# Last modified: 04/09/2024
# This script accompanies Steenberghen et. al. 2024 
# It is used to analyse fluorescence dynamics after a hyperosmotic and hypoosmotic
# treatment, respectively reducing and increasing fluorescence after the fixed 
# injection times at 10s and 50s.
# ----------------------------------------------------------------------------

# General guidelines:
# 1. Store R file in MAIN folder
# 2. Copy IJ output folders to MAIN folder and rename with format: "Output_R1", separate by replicate/run
# 3. Adapt metadata file and store in MAIN folder
# 4. For plots like plate-wise heat maps, also keep the MW96.csv file in the MAIN folder
# 5. To include ICC data for single-cell correlation etc. store those IJ summary files in the ICC subfolder of MAIN
# A dataset folder and a plot folder will be created automatically, plots are saved as .svg files by default
# datasets are save as .csv files
# PLOTS: there is a fixed general parameter called transp_. which can be used to produce transparent plots. 
# This parameter also helps to reduce standard parameters within the ggplot commands.
#
# METADATA: 
# Each well gets a unique ID that combines plate and well
# Each ROI gets a unique ID that combines plate, well and roi
# "Condition" = 0 marks images without stimulation for photobleaching correction, Condition = 1 are wells used for analysis
# "Treatment" and "Celltype" are parameters that can be used as grouping variables (factors)
# "Included" is a parameter that allows for quick well-wise exclusion from the metadata csv file. 
# Included can be loaded in on a later time in the script as well.

#----Preparations----
# Remove all objects from the environment
rm(list = ls())

# Trigger garbage collection to free up memory
gc()
#### Install and load packages ####

# List of required packages
packages <- c(
  "tidyverse", "Rmisc", "data.table", "gtools", "ggplot2", "ggthemes", "ggpubr", "reshape2", 
  "dplyr", "Hmisc", "RColorBrewer", "car", "inflection", "WRS2", 
  "greekLetters", "plotly", "numDeriv", "gmodels", "broom","progressr", "future","furrr", 
  "factoextra", "GGally", "parallel"
)

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    #install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Loop through the list of packages and install if missing
for (pkg in packages) {
  install_if_missing(pkg)
}


#### Definition of fixed variables and functions ####

dir = dirname(rstudioapi::getSourceEditorContext()$path)
calc.dir = paste0(dir,"/Calcein_Output")
plot.dir = paste0(dir,"/Plots")
datasets.dir = paste0(dir,"/Datasets")
dir.create(plot.dir)
dir.create(datasets.dir)
setwd(dir)

time.interval = 0.137           # sec (7.3 FPS)
reps = 3                        # Number of plates run
# # Assign colors SDB
# col_celltype_light = c("#0073ff","#f03737","#0bdb35") # for bar charts
# col_celltype = c("#038cfc","#ff5c33","#49d444") # darker color scheme -> for ribbons
# col_conditions_light = c("#52c447","#30752e","#6225cc","#7545e6") # for bar charts
# col_conditions = c("#129900","#0307fc","#02adad") # darker color scheme for 3 condition ribbons

# General parameter for transparent plots
transp =theme(
  plot.title = element_text(hjust = 0.5, size = 20),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background = element_blank(),
  legend.title = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))


#### Load treatment condition metadata ####
metadata = fread(paste0(dir,"/metadata.csv"),sep =";",stringsAsFactors = T)
metadata$Well = str_pad(metadata$Well,width = 3, side = "left",pad = "0")
metadata$Unique = paste0("R",metadata$Replicate,"W",metadata$Well)

#### Load and organize image data (+ merge with metadata) ####

im.data = data.table()
im.files = list.files(calc.dir,pattern = '*im_results')
n = length(im.files)
for (i in 1:n){
  print(im.files[i])
  path = paste(calc.dir,im.files[i],sep="/")
  results = fread(path, header = TRUE, sep = "\t", dec = ".")
  plate = substr(im.files[i],str_locate(im.files[i],"_R")[2]+1,str_locate(im.files[i],"_R")[2]+1)
  well = substr(im.files[i],str_locate(im.files[i],"NaCl_")[2]+1 ,str_locate(im.files[i],"NaCl_")[2]+3)
  results$Time = (results$V1) * time.interval
  # results$Norm.Mean = norm.fun(results$Mean1,results$Time)
  results$Plate = plate
  results$Well = well
  results$File = im.files[i]
  im.data = rbind(im.data,results)
}
im.data$Unique = paste0("R",im.data$Plate,"W",im.data$Well)

im.data = im.data %>%
  select(-Well) %>%
  left_join(x=., metadata) %>%
  rename(Im_Area= Area1, Im_Mean = Mean1, Nr = V1)

im.data$File = factor(im.data$File)
im.data$Plate = factor(im.data$Plate)
im.data$Well = factor(im.data$Well)
im.data$Condition = factor(im.data$Condition)

#### Store im.data ####
fwrite(im.data, paste0(datasets.dir,"/im.data.CAM6.csv"),dec= ".")


#### Load and organize Calcein roi data (+ merge with metadata) ####
roi.data = data.table()
roi.files = list.files(calc.dir,pattern = '*roi_results')
n = length(roi.files) # nr of files in directory
for (i in 1:n) # loop over all summary files
{
  path = paste(calc.dir,roi.files[i],sep="/")
  results = fread(path, header = TRUE, sep = "\t", dec = ".")
  print(i)
  print(path)
  instances = dim(results)[1]
  rois = (dim(results)[2]-1) /2 #(2 params are extracted, mean and area)
  melt.results = melt.data.table(results,"V1",variable.factor=T) # merge variables into a single column
  filt.results = melt.results[grep("Mean", melt.results$variable),] # only keep the intensity info
  filt.results_area = melt.results[grep("Area", melt.results$variable),] # only keep the area info
  filt.results = cbind(filt.results, filt.results_area)

  plate = substr(roi.files[i],str_locate(roi.files[i],"_R")[2]+1,str_locate(roi.files[i],"_R")[2]+1)
  well = substr(roi.files[i],str_locate(roi.files[i],"NaCl_")[2]+1 ,str_locate(roi.files[i],"NaCl_")[2]+3)
  filt.results$ROI =  rep(1:rois, each = instances)
  filt.results$Time = (filt.results$V1) * time.interval
  filt.results$Plate = plate
  filt.results$Well = well
  #filt.results$File = roi.files[i]
  filt.results$Unique = paste0("R",filt.results$Plate,"W",filt.results$Well)
  filt.results$V1 = NULL
  roi.data = rbind(roi.data,filt.results)
}

#check number of columns!
colnames(roi.data)[3] = "Nr"
colnames(roi.data)[2] = "Mean_value"
colnames(roi.data)[5] = "Area_value"
roi.data$variable = NULL
roi.data$variable = NULL

roi.data$File = NULL
roi.data$Well = NULL

# Combine with metadata
roi.data = left_join(roi.data,metadata %>% select(-c("Replicate","Well")),by="Unique")
roi.data = roi.data %>% rename(Well = Well_new)

# Assign factor labels
#roi.data$File = factor(roi.data$File)
roi.data$Plate = factor(roi.data$Plate)
#roi.data$Well = factor(roi.data$Well)
roi.data$Condition = factor(roi.data$Condition)
roi.data$Celltype = factor(roi.data$Celltype)

#Cleanup:
filt.results = NULL
filt.results_area = NULL
melt.results = NULL


#### Store roi.data ####
fwrite(roi.data, paste0(datasets.dir,"/roi.data.CAM6_unfiltered.csv"),dec= ".")

#----Data loading --------------------------------------------------------------
#### Read in roi.data ####
roi.data = fread(paste0(datasets.dir,"/roi.data.CAM6_unfiltered.csv"),dec= ".",stringsAsFactors = T, drop="V1")
roi.data$unique = paste0(roi.data$Unique,"T",roi.data$Time)

#### Read in im.data ####
im.data = fread(paste0(datasets.dir,"/im.data.CAM6.csv"),dec= ".",stringsAsFactors = T)
im.data$unique = paste0(im.data$Unique,"T",im.data$Time)


#### Combine ROI and image data ####
data = left_join(roi.data, im.data %>% select(Im_Area,Im_Mean,unique), by = "unique")


#### Save combined dataset ####
fwrite(data, paste0(datasets.dir,"/roi-im_data_unfiltered.csv"))


#### Read in roi.data ####
roi.data = fread(paste0(datasets.dir,"/roi.data.CAM6_unfiltered.csv"),dec= ".",stringsAsFactors = T, drop="V1")

#### Cell and well counts before filtering####
counts = roi.data %>%
  na.omit() %>%
  group_by(Unique, Celltype, ROI, Condition) %>%
  summarise(mean = mean(Mean_value)) %>%
  group_by(Unique,Celltype, Condition) %>%
  summarise(n_cells = n())

counts_well = counts %>%
  group_by(Celltype, Condition) %>%
  summarise(n_cells = sum(n_cells), n_wells = n())
fwrite(counts_well,paste0(datasets.dir,"/roi.data.counts_unfiltered.csv"))
counts = NULL
counnts_well = NULL

#### Filter out wells with technical issues ####
roi.data = roi.data %>%
  filter(Included == 1)

#### Cell and well counts after filtering####
counts = roi.data %>%
  na.omit() %>%
  group_by(Unique, Celltype, ROI, Condition) %>%
  summarise(mean = mean(Mean_value)) %>%
  group_by(Unique,Celltype, Condition) %>%
  summarise(n_cells = n())

counts_well = counts %>%
  group_by(Celltype, Condition) %>%
  summarise(n_cells = sum(n_cells), n_wells = n())
fwrite(counts_well,paste0(datasets.dir,"/roi.data.counts_filtered.csv"))


#### 28/06/2023 plot photobleaching curve +-SD M1M23
plot = roi.data %>%
  filter(Condition == 0, Celltype == "WT",Included == 1, Time <90) %>%
  mutate(norm.mean = Mean_value/mean(Mean_value[Time<1])) %>%
  group_by(Plate, Well_in_name,Time) %>%
  summarise(norm.mean = mean(norm.mean), norm.median = median(norm.mean))

# Calculate final value of photobleaching curve (81,8%)
# plot %>%
#   ungroup() %>%
#   filter(Time >88) %>%
#   summarise(mean = mean(norm.median), sd = sd(norm.median))

#col = c("goldenrod1", "darkred")
col=c("grey30","grey50","grey80")
col = "grey70"
# plot$Celltype = factor(plot$Celltype,levels =c("WT","M1M23"))


## Photobleaching curves

counts_well_filt = counts_well %>%
  filter(Celltype == "WT", Condition == 0)

svg(paste0(plot.dir,"/CAM6_photobleaching_wells_time-curve_norm_well-mean-se_M1M23.svg"),bg="transparent")
ggplot(plot,aes(x=Time,y=norm.mean)) + #, fill = as.factor(Plate)
  stat_summary(fun.data = mean_se, geom="ribbon", fun.args = list(mult = 1), alpha = 0.7) +
  stat_summary(fun = mean, geom="line", size = 1) +
  labs(x = "Time (s)", y = "Raw calcein intensity")+
  scale_fill_manual(values = col )+
  ggtitle("Dynamic range 1321N1 without stimulation") +
  theme_grey(base_size = 18)+
  transp
dev.off()





#### Read in counts datatable ####
counts_well = fread(paste0(datasets.dir, "/roi.data.counts_filtered.csv"))

#### Calculate background intensity ####
colnames(data)
str(roi.data)
well.BG = data %>%
  filter(Time == 0.137) %>%
  ungroup() %>%
  group_by(Plate, Well_in_name, Celltype, Condition) %>%
  summarise(roi.Area = sum(Area_value), roi.Mean = mean(Mean_value), im.Mean = Im_Mean, im.Area = Im_Area) %>%
  distinct() %>%
  mutate(BG.Mean = (im.Mean*im.Area - roi.Mean*roi.Area)/(im.Area-roi.Area)) %>%
  select(Plate, Well_in_name, Celltype, Condition, BG.Mean)

fwrite(well.BG, paste0(datasets.dir,"/Well_BG.csv"))

#### Plot BG intensity over plates ####
plot = well.BG

plot$Celltype = factor(plot$Celltype,levels = c("WT","M1M23"))
# Visualisatie parameters
col=c("Grey70","grey30")


svg(paste0(plot.dir,"/Mean_BG_Int_at_Baseline_filtered.svg"),bg="transparent")
ggplot(plot,aes(x = as.factor(Celltype),y = BG.Mean, fill =  as.factor(Celltype))) +
  stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1) , geom ="errorbar",position=position_dodge(width = 0.9, preserve = "single"),width=.5) +
  stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"),colour="white",width=0.9) +
  scale_fill_manual(values=c(col), guide=guide_legend(nrow=2))+
  facet_wrap(.~Plate)+ #, scales = 'free') , labeller(Plate = "Plate")
  theme_minimal(base_size = 12)+
  ggtitle("Mean background intensity at baseline")+
  labs(x = "Plate", y = "Mean intensity",caption = paste0("Ribbon represents SD (n = +-",round(mean(counts_well$n_wells))," wells / ",reps," plates)"))+
  theme(legend.position = "bottom", panel.spacing = unit(2,"lines"),axis.text.x=element_blank())+
  transp
dev.off()

#----Data normalization steps --------------------------------------------------
#### Normalisation to first 1s for both roi.data and data####
sec1 = roi.data %>% 
  filter(Time<1) %>% 
  group_by(Unique, ROI) %>% 
  summarise(sec1 = mean(Mean_value))

roi.data = roi.data %>% 
  left_join(., sec1, by = c("Unique","ROI")) %>% 
  mutate(Norm.Mean= Mean_value/sec1 ) %>% 
  select(-sec1)
sec1 = NULL

# For photobleaching controls, average out the intensity of all ROIs within one Well
roi.data.pnorm = roi.data %>%
  dplyr::group_by(Plate,Celltype,Condition, Time,  Well_in_name) %>%
  dplyr::filter(Condition==0) %>%
  dplyr::summarise(av.roi.norm = mean(Norm.Mean))

#Summarise photobleaching control per celltype, plate and condition (not well), average out the intensity of every well within celltype/plate/condition
roi.data.pnorm = roi.data.pnorm %>%
  dplyr::group_by(Plate, Celltype, Time) %>%
  dplyr::summarise(roi.pnorm.value = mean(av.roi.norm))

# For im-v-roi dataset:
# roi.data.pnorm %>% 
#   filter(Time >89) %>% 
#   group_by(Plate, Celltype) %>% 
#   summarise(roi.mean = mean(roi.pnorm.value), im.mean=mean(im.pnorm.value)) %>% 
#   mutate(frac = im.mean/roi.mean)

#### Plots M1M23 image vs roi-based intensity trace ####
## im+roi data processing
data =data %>% 
  mutate(unique = paste0(Unique,"R",ROI))

str(data)
data = data %>%
  group_by(Unique, ROI) %>%
  mutate(roi.Norm.Mean = Mean_value / (mean(Mean_value[Time<1])), im.Norm.Mean = Im_Mean / (mean(Im_Mean[Time<1]))) %>% 
  data.table()

# data.pnorm = data %>%
#   dplyr::group_by(Plate,Celltype,Condition, Time,  Well_in_name) %>%
#   dplyr::filter(Condition==0) %>%
#   dplyr::summarise(av.roi.norm = mean(roi.Norm.Mean), av.im.norm = im.Norm.Mean)
# 
# data.pnorm = data.pnorm %>%
#   dplyr::group_by(Plate, Celltype, Time) %>%
#   dplyr::summarise(roi.pnorm.value = mean(av.roi.norm), im.pnorm.value = mean(av.im.norm))

plot = data %>%
  filter(Time < 90) %>% 
  group_by(Plate, Condition, Time,Celltype, Well_in_name)%>%
  summarise(av.roi.norm = mean(roi.Norm.Mean), av.im.norm = im.Norm.Mean)

#col = c("goldenrod1", "darkred")
col=c("grey70","grey30")

# plot$Celltype = factor(plot$Celltype,levels =c("WT","M1M23"))


## Photobleaching curves
plot_long = plot %>% 
  filter(Celltype == "M1M23", Condition==0) %>% 
  pivot_longer(cols= c(av.im.norm,av.roi.norm))

counts_well_filt = counts_well %>% 
  filter(Celltype == "M1M23", Condition == 0)

svg(paste0(plot.dir,"/CAM6_M1M23_im-v-roi_photobleaching_wells_time-curve_norm_well-mean-sd.svg"),bg="transparent")
ggplot(plot_long,aes(x=Time,y=value, fill = name)) +
  stat_summary(fun.data = mean_sdl, geom="ribbon", fun.args = list(mult = 1), alpha = 0.95) +
  stat_summary(fun = mean, geom="line", size = 1) +
  labs(x = "Time (s)", y = "I/I0",caption = paste0("Ribbon represents inter-well SD \n3 repeats \n",counts_well_filt$n_wells," wells \n",counts_well_filt$n_cells," cells"))+
  scale_fill_manual(labels=c("Image mean", "Cell mean"),values = col )+
  ggtitle("Dynamic range 1321N1 without stimulation") +
  theme_grey(base_size = 15)+
  transp
dev.off()

## Stimulation curves
plot_long = plot %>% 
  filter(Celltype == "M1M23", Condition==3) %>% 
  pivot_longer(cols= c(av.im.norm,av.roi.norm))

counts_well_filt = counts_well %>% 
  filter(Celltype == "M1M23", Condition == 3)

svg(paste0(plot.dir,"/CAM6_M1M23_condition3_im-v-roi_stimulation_wells_time-curve_norm_well-mean-sd.svg"),bg="transparent")
ggplot(plot_long,aes(x=Time,y=value, fill = name)) +
  stat_summary(fun.data = mean_sdl, geom="ribbon", fun.args = list(mult = 1), alpha = 0.95) +
  stat_summary(fun = mean, geom="line", size = 1) +
  labs(x = "Time (s)", y = "I/I0")+ #,caption = paste0("Ribbon represents inter-well SD \n3 repeats \n",counts_well_filt$n_wells," wells \n",counts_well_filt$n_cells," cells")
  scale_fill_manual(labels=c("Image mean", "Cell mean"),values = col )+
  ggtitle("Dynamic range 1321N1-M1M23 with stimulation") +
  theme_grey(base_size = 15)+
  transp
dev.off()


ggplot(plot_long %>% filter(name =="av.roi.norm"),aes(x=Time,y=value, col = Well_in_name)) +
  #stat_summary(fun.data = mean_sdl, geom="ribbon", fun.args = list(mult = 1), alpha = 0.8) +
  geom_line(size=1)+
  #stat_summary(fun = mean, geom="line", size = 1) +
  #labs(x = "Time (s)", y = "Mean well intensity",caption = paste0("Ribbon represents inter-well SD \n3 repeats \n",round(mean(counts_well$n_wells))," (+-",sd(counts_well$n_wells), ") wells \n",round(mean(counts_well$n_cells))," (+-",round(sd(counts_well$n_cells)), ") roi"))+
  #scale_fill_manual(values = col )+
  ggtitle("Raw intensity photobleaching wells") +
  facet_grid(Plate~.) +
  theme_grey(base_size = 15)+
  transp

# individual traces of samples cells 29/8/2023

hist(data$ROI)

plot_samples = data%>% 
  filter(Time <0.2, Included == 1, Condition == 3) %>% 
  group_by(Celltype,Plate) %>% 
  slice_sample(n = 50, replace = FALSE) %>% 
  ungroup() %>% 
  dplyr::select(unique)

plot = data %>%
  filter(Time<90, Included==1, unique %in% plot_samples$unique)

plot_long = plot %>% 
  filter(Celltype == "M1M23", Condition==3, Time< 25) %>% 
  pivot_longer(cols= c(im.Norm.Mean,roi.Norm.Mean))

svg(paste0(plot.dir,"/CAM6_M1M23_condition3_stimulation_150_samples_time-curve_norm.svg"),bg="transparent")
ggplot(plot_long %>% filter(name =="roi.Norm.Mean"),aes(x=Time,y=value)) +
  stat_summary(fun = mean, geom="line", color = "red", size = 1) +
  stat_summary(data = plot_long %>% filter(name =="im.Norm.Mean"), fun = mean, geom="line", color = "green", size = 1) +
  geom_point(size=0.5, alpha = 0.5)+
  #stat_summary(fun = mean, geom="line", size = 1) +
  labs(x = "Time (s)", y = "Normalised intensity",caption = paste0("50 sampled cells per plate,\nwith mean value of rois in red \nand mean value of corresponding well image means in green"))+
  #scale_fill_manual(values = col )+
  ggtitle("50 sampled 1321N1-M1M23 with stimulation per plate") +
  facet_grid(Plate~.) +
  theme_grey(base_size = 15)+
  transp
dev.off()




## roi.data processing
# append with normalized values (to first 7 time points (first second))
roi.data = roi.data %>%
  group_by(Unique, ROI) %>%
  mutate(roi.Norm.Mean = Mean_value / mean(Mean_value[Time<1])) %>% 
  data.table()


#### Photobleaching correction on ROI average values (using exponential fit) ####
# Initialize counters
success_count <- 0
failure_count <- 0

# Function to fit model and log success/failure
fit_model <- function(data, plate) {
  temp <- data %>% filter(Plate == plate)
  tryCatch({
    fit.roi <- nls(roi.pnorm.value ~ SSasymp(Time, yf, y0, log_alpha), data = temp)
    temp.roi <- augment(fit.roi)
    success_count <<- success_count + 1
    return(temp.roi)
  }, error = function(e) {
    failure_count <<- failure_count + 1
    return(NULL)
  })
}

# Fit exponential model for each celltype and plate
fit_WT <- roi.data.pnorm %>%
  filter(Celltype == "WT")
list.roi <- data.frame()

for (i in 1:3) {  # Number of plates
  temp.roi <- fit_model(fit_WT, i)
  if (!is.null(temp.roi)) {
    list.roi <- rbind(list.roi, temp.roi)
  }
}
fit_WT <- left_join(fit_WT, list.roi, by = c('roi.pnorm.value', 'Time'))

fit_M1M23 <- roi.data.pnorm %>%
  filter(Celltype == "M1M23")
list.roi <- data.frame()

for (i in 1:3) {  # Number of plates
  temp.roi <- fit_model(fit_M1M23, i)
  if (!is.null(temp.roi)) {
    list.roi <- rbind(list.roi, temp.roi)
  }
}
fit_M1M23 <- left_join(fit_M1M23, list.roi, by = c('roi.pnorm.value', 'Time'))

# Merge all fitted models together and join with data
exp.fit <- rbind(fit_WT, fit_M1M23)
roi.data.pnorm <- left_join(roi.data, exp.fit)

# Calculate success rate
total_attempts <- success_count + failure_count
success_rate <- success_count / total_attempts

# Write success rate to a text file
write(paste0("Fitting success rate: ", success_rate), file = paste0(datasets.dir, "/PB_success_rate.txt"))


#### Calculate pcorr.value ####
roi.data.pcorr = roi.data.pnorm%>%
  mutate(pcorr.value = Norm.Mean/.fitted)

roi.data.pcorr = roi.data.pcorr %>% select(-c('.fitted', '.resid','roi.pnorm.value'))

#### Store roi.data.pcorr ####

fwrite(roi.data.pcorr, paste0(datasets.dir,"/roi.data.pcorr.CAM6.csv"),dec= ".")

roi.data =NULL
roi.data.pnorm = NULL

#----Calculation of decay and recovery rate ------------------------------------
#### Read in roi.data.pcorr ####
roi.data.pcorr = fread(paste0(datasets.dir,"/roi.data.pcorr.CAM6.csv"), header = TRUE, sep = ",", drop = 1) %>% 
  mutate(unique = paste0(Plate,"_",Well,"_",ROI)) %>% 
  filter(Included == 1, Area_value>100,  Area_value<1500, Condition != 0, Condition != 7)


#### Bend point detection (NaCl) w Bisection Extremum Distance Estimator ####

# Filter time range
exp_inflection = roi.data.pcorr %>%
  filter(Time>9 & Time <40) %>%                 # set nr of frame of beginning NaCl administration and end of response
  group_by(Plate, Celltype, Well, ROI) %>%
  mutate(iplast_exp = bede(Time, pcorr.value,index=1)$iplast) %>%
  distinct(iplast_exp) %>% 
  na.omit()%>%
  ungroup()%>% 
  data.table()

svg(paste0(plot.dir,"/CAM7_hist_iplast_exp.svg"))
hist(exp_inflection$iplast_exp)
dev.off()

# add unique identifyer
exp_inflection = exp_inflection %>%
  dplyr::mutate(unique = paste0(Plate,"_",Well,"_",ROI) )
# merge with roi.data.pcorr
roi.data.pcorr = left_join(roi.data.pcorr, exp_inflection %>% select('iplast_exp','unique'), by = "unique")
exp_inflection =NULL

#### Bend point detection (H2O) w Bisection Extremum Distance Estimator ####
# Filter time range
rec_inflection = roi.data.pcorr %>%
  filter(Time>47 & Time <87) %>%                           # set nr of frame of beginning H2O administration and end of response
  group_by(Plate, Celltype,Well, ROI) %>%
  mutate(iplast_rec = bede(Time,pcorr.value,index=0)$iplast) %>%
  distinct(iplast_rec) %>% 
  drop_na() %>% 
  ungroup()

svg(paste0(plot.dir,"/CAM7_hist_iplast_rec.svg"))
hist(rec_inflection$iplast_rec)
dev.off()

rec_inflection = rec_inflection %>%
  dplyr::mutate(unique = paste0(Plate,"_",Well,"_",ROI) )
roi.data.pcorr = left_join(roi.data.pcorr, rec_inflection %>% select('iplast_rec','unique'), by = "unique")
rec_inflection =NULL
#### Model fit on NaCl slope (exponential fit) and H2O slope (sigmoidal fit) ####

# Replace set with roi.data.pcorr to use less RAM
set <- roi.data.pcorr
roi.data.pcorr <- NULL

# # DEBUG with small subset
# first_500_unique <- head(unique(set$unique), 500)
# set <- set[unique %in% first_500_unique]
 
# Set up parallel processing
num_cores <- detectCores()
print(paste("Number of available cores:", num_cores))

options(future.globals.maxSize = 2 * 1024^3)  # Increase to >1.4 GiB
plan(multisession, workers = num_cores-4)

set$unique <- factor(set$unique)
set <- set %>%
  distinct(Time, unique, pcorr.value, iplast_exp, iplast_rec) %>%
  drop_na() %>%
  data.table()

# Columns to keep
cols_to_keep <- head(names(set), -2)

handlers(global = TRUE)

for (j in c(10, 15,20)) {

  
  # Stimulation fit
  exp.fit_list <- vector("list", length(levels(set$unique)))
  gc()
  
  with_progress({
    message(sprintf("Fitting decay rate model for time window %s s", j))
    # Initialize progress bar for stimulation fit
    library(data.table)
    library(broom)
    pb_exp <- progressor(along = levels(set$unique))
    
    exp.fit_list <- future_map(levels(set$unique), function(unique_level) {
      pb_exp()
      temp <- set[Time > 10 & Time < iplast_exp + j & unique == unique_level, .SD, .SDcols = cols_to_keep]
      
      fit_result <- tryCatch({
        fit <- nls(pcorr.value ~ SSasymp(Time, yf, y0, log_alpha), data = temp)
        fit <- tidy(fit)
        data.table(unique = unique_level, log_alpha = fit$estimate[3], std.error = fit$std.error[3])
      }, error = function(e) {
        #message(sprintf("Error fitting model for unique_level %s: %s", unique_level, e$message))
        data.table(unique = unique_level, log_alpha = NA, std.error = NA)
      })
      
    })
  })

  exp.fit <- rbindlist(exp.fit_list)
  fwrite(exp.fit, paste0(datasets.dir, "/exp.fit-on-well_median_iplast_exp_", j, ".csv"), dec = ".")
  gc()

  # Recovery fit
  rec.fit_list <- vector("list", length(levels(set$unique)))

  with_progress({
    message(sprintf("Fitting recovery rate model for time window %s s", j))
    # Initialize progress bar for recovery fit
    pb_rec <- progressor(along = levels(set$unique))

    rec.fit_list <- future_map(levels(set$unique), function(unique_level) {
      library(data.table)
      library(broom)
      pb_rec()
      temp <- set[Time > 50 & Time < iplast_rec + j & unique == unique_level, .SD, .SDcols = cols_to_keep]

      fit_result <- tryCatch({
        fit <- nls(pcorr.value ~ SSasymp(Time, yf, y0, log_alpha), data = temp)
        fit <- tidy(fit)
        data.table(unique = unique_level, rec_alpha = fit$estimate[3], std.error = fit$std.error[3])
      }, error = function(e) {
        #message(sprintf("Error fitting model for unique_level %s: %s", unique_level, e$message))
        data.table(unique = unique_level, rec_alpha = NA, std.error = NA)
      })
    })
  })

  rec.fit <- rbindlist(rec.fit_list)
  fwrite(rec.fit, paste0(datasets.dir, "/rec.fit-on-well_median_iplast_rec_", j, ".csv"), dec = ".")
  gc()
}

#----Calculation of other parameters and aggregation----------------------------
#### Read in roi.data.pcorr ####
roi.data.pcorr = fread(paste0(datasets.dir,"/roi.data.pcorr.CAM6.csv"), header = TRUE, sep = ",", drop = 1) %>% 
  mutate(unique = paste0(Plate,"_",Well,"_",ROI)) %>% 
  filter(Included == 1, Area_value>100,  Area_value<1500, Condition != 0, Condition != 7)

#### Read in fit dataframes####
# Match the name to the timeframe (here 10s) used for the model fit
exp.fit = fread(paste0(datasets.dir,"/exp.fit-on-well_median_iplast_exp_10.csv"), header = TRUE, sep = ",")
rec.fit = fread(paste0(datasets.dir,"/rec.fit-on-well_median_iplast_rec_10.csv"), header = TRUE, sep = ",")

#### Create new combined data table from filtered roi data at roi level, removing time dimension ####
idvars_cell = c("ROI","Plate","Well_in_name", "unique", "Unique","Celltype", "Condition")

comb.data = roi.data.pcorr %>%
  select(idvars_cell) %>%
  distinct()
#### Read in fit data and combine with comb.data ####

exp.fit = fread(paste0(datasets.dir,"/exp.fit-on-well_median_iplast_exp_20.csv"), header = TRUE, sep = ",")
comb.data = left_join(comb.data, exp.fit %>% select(log_alpha, unique), by = "unique")       # verify column number!

rec.fit = fread(paste0(datasets.dir,"/rec.fit-on-well_median_iplast_rec_10.csv"), header = TRUE, sep = ",")
comb.data = left_join(comb.data, rec.fit %>% select(rec_alpha, unique), by = "unique")       # verify column number!

comb.data = comb.data %>% drop_na()

#### slope calculation ####
#filter timeframe of response + calculate slope + add as variable

slope_exp = roi.data.pcorr %>%
  dplyr::group_by(Celltype, Plate, Well, ROI) %>%
  filter(Time > iplast_exp-0.5 & Time < iplast_exp+0.5) %>%
  dplyr::mutate(slope = paste0(lm(pcorr.value ~ Time)$coeff[[2]]))%>%
  #dplyr::mutate(intersect = paste0(lm(pcorr.value ~ Time)$coeff[[1]]))%>%
  distinct(slope)%>%
  dplyr::mutate(unique = paste0(Plate,"_",Well,"_",ROI) )%>%
  na.omit()

slope_rec = roi.data.pcorr %>%
  dplyr::group_by(Celltype,Plate, Well, ROI) %>%
  filter(Time > 50 & Time < iplast_rec+2) %>%
  dplyr::mutate(slope = paste0(lm(pcorr.value ~ Time)$coeff[[2]]))%>%
  #dplyr::mutate(intersect = paste0(lm(pcorr.value ~ Time)$coeff[[1]]))%>%
  distinct(slope)%>%
  dplyr::mutate(unique = paste0(Plate,"_",Well,"_",ROI) )%>%
  na.omit()

comb.data = left_join(comb.data, slope_exp[,c(5,6)], by = "unique")
comb.data = left_join(comb.data, slope_rec[,c(5,6)], by = "unique")
comb.data = comb.data%>%
  rename(slope_exp = slope.x, slope_rec = slope.y)%>%
  na.omit()


#### Minimum calculation ####
roi.data.pcorr = roi.data.pcorr %>%
  mutate(unique = paste0(Plate,"_", Well,"_",ROI))

min= roi.data.pcorr %>%
  filter(Condition != 0 & Condition != 7 & Time >50) %>%
  group_by(Celltype, Plate, Well, ROI, unique) %>%
  summarise(min = pcorr.value[which.min(pcorr.value)])

comb.data = left_join(comb.data, min[,c(5,6)], by = "unique")

comb.data$slope_exp <- as.numeric(comb.data$slope_exp)
comb.data$slope_rec <- as.numeric(comb.data$slope_rec)
comb.data$min <- as.numeric(comb.data$min)


#### Plot: how good is the slope calculation? ####

### plot data with slope
##set ROI for which you want to plot slope
roi.data.pcorr = roi.data.pcorr%>%
  dplyr::mutate(unique = paste0(Plate, "_",Well,"_",ROI) )

# plot.data = roi.data.pcorr %>%
#   filter(unique == "2_013_64")
# plot.data = comb.data %>%
#   filter(unique == "2_013_64")

plot.data = roi.data.pcorr %>%
  filter(well_label == "3_019") %>% group_by(Time) %>% mutate(value = mean(pcorr.value)) %>% distinct(Time,value)

## Calculate slope and intercept for that ROI
plot.data$slope_rec[1]
#plot.data$intersect_exp[1]
plot.data$iplast_exp[1]

## Plot curve for that ROI and add slope
ggplot(plot.data, aes(x = Time, y = value)) +
  ggtitle(paste0("pcorr value vs time for CAM6 "))+
  geom_line(size = 2) +
  coord_cartesian( xlim = c(10, 75), ylim = c(0, 1))+
  geom_vline(xintercept =10, colour="red")+
  geom_vline(xintercept =50, colour="cornflowerblue")+
  geom_vline(xintercept = 10.001)+
  geom_abline(slope = 0.0251167	, intercept = -0.5, colour = "green")+               ###fill in slope and intercept
  transp

##plot correlation between slope and inflection point
data.plot = comb.data%>%
  distinct(slope_exp, slope_rec, iplast_exp, iplast_rec, unique, Celltype)
ggplot(data.plot, aes(slope_rec, iplast_rec, col = Celltype))+
  geom_point()+
  coord_cartesian(xlim = c(0,0.11), ylim = c(49, 61))+
  ggtitle("Comparison of slope_exp vs. iplast_exp")+
  transp


#### Maximum calculation ####


max_well= roi.data.pcorr %>%
  filter(Condition != 0 & Condition != 7 & Time >50) %>%
  group_by(Celltype, Plate, Well, ROI) %>%
  summarise(max = pcorr.value[which.max(pcorr.value)])%>%
  dplyr::mutate(unique=paste0(Plate,"_",Well,"_",ROI))

comb.data = left_join(comb.data, max_well[,c(5,6)],by = 'unique')


#### Recovery calculation ####
hist(temp$mean_start)

plot = roi.data.pcorr %>%
  filter(unique == "2_008_1")

norm.mean <- plot[, list(norm.mean = Mean_value / (sum(Mean_value[Time<1]) / length(Mean_value[Time<1]))), by=key(plot)]
plot =plot %>%
  mutate(norm.mean =  Mean_value / (sum(Mean_value[Time<1]) / length(Mean_value[Time<1])))

plot(plot$Time, plot$Mean_value)
plot(plot$Time, plot$Norm.Mean)
plot(plot$Time, plot$norm.mean)
plot(plot$Time, plot$pcorr.value)
plot(plot$Time, plot$pcorr.value_new)

plot = plot %>%
  mutate(pcorr.value_new = pcorr.value/max(pcorr.value))

temp = roi.data.pcorr%>%
  filter(Time>0 & Time < 9)%>%
  group_by(unique)%>%
  mutate(mean_start = mean(pcorr.value))%>%
  distinct(unique, mean_start)
comb.data = left_join(comb.data, temp, by = "unique")

comb.data = comb.data%>%
  group_by(unique)%>%
  na.omit%>%
  #mutate(rec_baseline = (max / mean_start)*100)%>%
  #mutate(rec_min = (max - min)*100)%>%
  mutate(signal_rec = (max-min)/(mean_start-min))


#### Extra parameters: dynamic range, recovery vs min, overshoot #### 

group_var = "Condition"
id_vars_cell = c("Celltype", group_var,"Unique", "unique")
id_vars = c("Celltype", group_var,"Unique")
num_vars = c("dyn_range","recovery_end","overshoot")
group_levels = c("1","2","3","4")

extra_params = roi.data.pcorr %>%
  filter(Condition != 0 & Time >45) %>%
  group_by_at(c(id_vars_cell, group_var)) %>%
  summarise(tmax = Time[which.max(pcorr.value)],
            tend = mean(tail(Time, 7)),
            min = min(pcorr.value),
            max = max(pcorr.value),
            total_max = mean(ifelse(max>1, max, 1)),
            dyn_range = mean((total_max-min)*100),
            end = mean(tail(pcorr.value, 7)),
            slope_RVD = (end-max)/(tend-tmax),
            overshoot = (max-end)/end*100,
            recovery_end = (end-min)/(1-min)*100) %>% 
  distinct() %>% 
  ungroup()

# summary at well level
data_temp = extra_params %>% 
  select(-tend, -tmax, -max, -total_max,-end,-min) %>% 
  group_by_at(factor(id_vars)) %>% 
  summarise_if(is.numeric, median)


# set to long format for looped printing
data_temp_long = data_temp %>% 
  pivot_longer(-all_of(id_vars), names_to = "variables", values_to = "value") %>% 
  drop_na() 
  

# Plotting vars

transp =
  theme_grey(base_size = 22)+
  theme(
  plot.title = element_text(hjust = 0.5, size = 20),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background = element_blank(),
  legend.title = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA),
  aspect.ratio = 3)

# col = c("grey70","grey30")
col = c("grey30","#0075D4")
col_trt = c("grey20","grey45","grey60","grey75")

plot(as.factor(RVD$Plate),RVD$slope_RVD)

for(var.select in num_vars){

  # Plot both celltypes, treatments
  plot = data_temp_long %>% filter(variables == var.select)
  plot$Celltype = factor(plot$Celltype, levels = c("WT","M1M23"))
  plot$variables = factor(plot$variables, levels = group_levels)
  
  
  svg(paste0(plot.dir,"/CAM6_bar_plot_well_median-se_",var.select,"_bothCelltype.svg"), bg = "transparent")
  p = ggplot(plot,aes(x = as.factor(Condition),y = value, fill =  as.factor(Celltype))) +
    stat_summary(fun.data = mean_se , geom ="errorbar",position=position_dodge(width = 0.9, preserve = "single"),width=.5) +
    stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"),colour="white",width=0.9) +
    scale_fill_manual(values=c(col), guide=guide_legend(nrow=2))+
    xlab("Well mean parameters")+
    ylab("")+
    facet_wrap(Condition~., scales = 'free')+ #,space = "free")+
    theme(legend.position = "bottom")+
    transp
  print(p)
  
  dev.off()
  
  # Plot M1M23 only
  plot = data_temp_long %>% filter(variables == var.select, Celltype=="M1M23")
  plot$Celltype = factor(plot$Celltype, levels = c("WT","M1M23"))
  
  
  svg(paste0(plot.dir,"/CAM6_bar_plot_well_median-se_",var.select,"_M1M23.svg"), bg = "transparent")
  p = ggplot(plot ,aes(x = as.factor(Condition),y = value, fill= as.factor(Condition))) +
    stat_summary(fun.data = mean_se , geom ="errorbar",position=position_dodge(width = 0.9, preserve = "single"),width=.5) +
    stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"),colour="white",width=0.9) +
    scale_fill_manual(values=col_trt)+
    xlab("Well mean parameters")+
    ylab("")+
    facet_wrap(.~variables, scales = 'free')+ #,space = "free")+
    scale_y_continuous(expand = expansion(c(0.2, 0.3))) +
    theme(legend.position = "bottom")+
    transp
  print(p)
  dev.off()

  # ANOVA both celltypes
  sink(paste0(plot.dir,"/CAM6_well_",var.select,"_anova.txt"))
  print(var.select)
  print("---------------Anova------------------")
  test = aov(get(var.select)~as.factor(Condition)*Celltype, data_temp)
  Tukey = TukeyHSD(test)
  print(Tukey)
  sink()
  
  
  # ANOVA M1M23 only
  sink(paste0(plot.dir,"/CAM6_M1M23_well_",var.select,"_anova.txt"))
  print(var.select)
  print("---------------Anova------------------")
  test = aov(get(var.select)~as.factor(Condition), data_temp)
  Tukey = TukeyHSD(test)
  print(Tukey)
  sink()
}

# #### Clean up dataset ####
# 
# comb.data = left_join(comb.data, exp.fit[,c(1,2)])
# comb.data = left_join(comb.data, rec.fit[,c(1,2)])
# #complete.data = comb.data
# # comb.data = comb.data%>%
# #   distinct(ROI, Plate,Well, Celltype, unique, iplast_exp, iplast_rec, slope_exp, slope_rec, min, max, mean_start, rec_baseline, rec_min, signal_rec, log_alpha, rec_alpha)
# # comb.data = comb.data %>% mutate(well_label = paste0(Plate,"_",Well))
# # comb.data = left_join(comb.data, roi.data.pcorr[,c(9,14)])
# # test = comb.data %>% filter(Celltype == "M1M23") %>% distinct(well_label)
# 
# #### write combined datasets to HDD ####
# fwrite(comb.data, paste0(datasets.dir,"/comb.data.CAM6.csv"),dec= ".")
# comb.data = fread(paste0(dir,"/Datasets/comb.data.CAM6.csv"), header = TRUE, sep = ",", drop = 1,stringsAsFactors = T)
# 
# #### Cell and well counts after filtering comb data ####
# counts = comb.data %>%
#   na.omit() %>%
#   group_by(Unique, Celltype, ROI, Condition) %>%
#   summarise(mean = mean(log_alpha)) %>%
#   group_by(Unique,Celltype, Condition) %>%
#   summarise(n_cells = n())
# 
# counts_well = counts %>%
#   group_by(Celltype, Condition) %>%
#   summarise(n_cells = sum(n_cells), n_wells = n())
# fwrite(counts_well,paste0(datasets.dir,"/comb.data.counts_filtered.csv"))
# counts = NULL
# #counts_well = NULL
# 
# 
# 
# 
# 
# 
# # Remove time variable
# comb.data = comb.data %>% distinct(Area_value, ROI, Plate, Well, Condition, Condition_name, Celltype, Included, well_label, unique, iplast_exp, iplast_rec, slope_exp, slope_rec, min, max, mean_start, rec_baseline, rec_min, signal_rec, log_alpha, rec_alpha)
# 
# fwrite(comb.data, paste0(datasets.dir,"/comb.data_notime.CAM6.csv"),dec= ".")
# 
#### Read in combined data ####
comb.data = fread(paste0(dir,"/Datasets/comb.data.CAM6.csv"), header = TRUE, sep = ",", drop = 1,stringsAsFactors = T)
comb.data =comb.data %>% mutate(Decay_rate = -log_alpha, Recovery_rate=-rec_alpha) %>% 
  rename(Baseline = mean_start,Recovery= signal_rec) %>% 
  select(-c(log_alpha, rec_alpha))


#----Statistics and plots ------------------------------------------------------
#### Statistics plots SDB ####
#calcein stabilty

calcein = roi.data.pcorr %>%
    group_by(well_label)%>%
    filter(Condition != 0 & Condition != 7, Time < 2)%>%
    mutate(well_mean = mean(Mean_value))%>%
    group_by(Celltype, Condition)%>%
    mutate(Calcein = mean(well_mean))%>%
    mutate(sd = sd(well_mean))%>%
    distinct(Condition_name, Calcein, sd, Celltype, well_mean)

  positions <- c("WT", "M1M23")
  calcein$Celltype <- factor(calcein$Celltype, levels = c("WT","M1M23"))
  
svg("fig9.svg", bg = "transparent", width = 6, height = 5)
  ggplot(calcein, aes(x = Condition_name, y = Calcein, fill = Celltype))+
    geom_col(stat = "identity", position = 'dodge')+
    scale_fill_manual(values = c("goldenrod1", "darkred"))+
    ggtitle("Calcein stability at baseline")+
    labs(
      x = "Condition",
      y = "Mean raw baseline calcein intensity",
      caption = "Bars indicate standard deviation")+
    geom_errorbar(stat = "identity", position = position_dodge(.9),aes(ymin = Calcein-sd,  ymax = Calcein+sd),color = "#22292F", width = .1) +
    theme_grey(base_size = 25)+
  transp
dev.off()

test = calcein %>% filter(Celltype == 'M1M23', Condition_name == "T3")
str(plate)
svg("qq_plot_calcein_stability.svg", bg = "transparent")
qqnorm(calcein$well_mean, pch = 1, frame = FALSE)
qqline(weighted_av_rec$weighted_av, col = "steelblue", lwd = 2)
dev.off()



# Calculate median of every variable
median= left_join(comb.int.data, AQP4.data[,c(2,20,38)], by = "unique")
median = left_join(median, roi.data.pcorr[,c(15,9,10,14)])
str(median)
temp = median %>% distinct(Celltype, Condition, Condition_name, well_label,Circ.)

plot = temp %>%
  mutate(cond = paste0(Condition_name,"_",Celltype))%>%
  group_by(cond, Condition, Celltype, well_label)%>%
  na.omit()%>%
  mutate(mean_well = mean(Circ.))%>%     ## fill in variable
  group_by(cond, Condition, Celltype)%>%
  mutate(median = mean(mean_well))%>%
  mutate(sd = sd(mean_well))%>%
  distinct(sd, median, Condition_name, Condition, mean_well, cond, Celltype, well_label)

plot_2 = temp %>%
  group_by(well_label)%>% na.omit() %>%
  mutate(mean_well = mean(min))%>%     ## fill in variable
  group_by(Celltype)%>%
  mutate(median = mean(mean_well))%>%
  mutate(sd = sd(mean_well))%>%
  distinct(sd, median, Celltype, mean_well,well_label)

test = plot %>% filter(cond == "T1_M1M23")
qqnorm(test$mean_well, pch = 1, frame = FALSE)
qqline(weighted_av_rec$weighted_av, col = "steelblue", lwd = 2)


test = plot %>% filter(Celltype == "M1M23") %>% distinct(mean_well, Condition, well_label)
test$Condition <- as.factor(test$Condition)
kruskal.test(mean_well ~ Condition, data = test)
pairwise.wilcox.test(test$mean_well, test$cond,
                     p.adjust.method = "BH")


plot(test$Decay_rate)

test = plot %>% dplyr::filter(Celltype == "M1M23") %>% distinct(mean_well, Condition, well_label)
test$Condition <- as.factor(test$Condition)
anova <- aov(mean_well ~ Condition, data = test)
TukeyHSD(anova)

t.test(mean_well~Celltype, data = plot_2, equal.var = F)
wilcox.test(mean_well ~Celltype, data = plot_2, equal.var = F)

stat = plot %>% 
  distinct(sd, median, Celltype, Condition_name)


positions <- c("WT", "M1M23")
stat$Celltype <- factor(stat$Celltype, levels = c("WT","M1M23"))

svg("fig20b.svg", bg = "transparent", width = 6, height = 5)
ggplot(stat, aes(x = Condition_name, y = median, fill = Celltype))+
  geom_col(stat = "identity", position = 'dodge')+
  scale_fill_manual(values = c("goldenrod1", "darkred"))+
  ggtitle("Stability of the cell circularity ")+
  labs(
     x = "Condition",
     y = "Mean cell circularity",
     caption = "Bars indicate standard deviation")+
  geom_errorbar(stat = "identity", position = position_dodge(.9),aes(ymin = median-sd,  ymax = median+sd),color = "#22292F", width = .1) +
  theme_grey(base_size = 25)+
  #scale_y_reverse()+
  transp
dev.off()

svg("fig17.svg", bg = "transparent", width = 6, height = 5)
ggplot(stat, aes(x = Celltype, y = median, fill = Celltype))+
  geom_col(stat = "identity", position = 'dodge')+
  scale_fill_manual(values = c("goldenrod1", "darkred"))+
  ggtitle("Slope-rec value")+
  labs(
    x = "Celltype",
    y = "Mean slope-rec value",
    caption = "Bars indicate standard deviation")+
  geom_errorbar(stat = "identity", position = position_dodge(.9),aes(ymin = median-sd,  ymax = median+sd),color = "#22292F", width = .1) +
  theme_grey(base_size = 25)+
  #scale_y_reverse()+
  transp
dev.off()

test = comb.int.data %>% mutate(comb = paste0(Celltype,"_",Condition_name)) %>% distinct(comb, well_label)
table(test$comb)
# ggplot(stat, aes(x = Condition_name, y = median))+
#    geom_col()+
#    #coord_cartesian(ylim = c(0, 30))+
#    ggtitle("Median slope after NaCl stimulation (WT)")+
#    labs(
#        x = "Condition",
#        y = "-Median slope after NaCl stimulation",
#        caption = "Bars indicate standard deviation")+
#    geom_errorbar(aes(ymin = median-sd,  ymax = median+sd),color = "#22292F", width = .1) +
#    # geom_line(data = p_value_one,
#    #                  aes(x = x, y = y, group = 1)) +
#    # geom_line(data = p_value_two,
#    #                  aes(x = x, y = y, group = 1)) +
#    # geom_line(data = p_value_three,
#    #                  aes(x = x, y = y, group = 1)) +
#    # annotate("text", x = 2, y = 1.25,
#    #          label = "p =0.0293",
#    #          size = 4, color = "#22292F") +
#    # annotate("text", x = 2.5, y = 1.35,
#    #           label = "p = 0.004329",
#    #           size = 4, color = "#22292F") +
#    # annotate("text", x = 3, y = 1.45,
#    #          label = "p = 0.001166",
#    #          size = 4, color = "#22292F") +
#   transp

plot = roi.data.pcorr %>%
  group_by(well_label, Time)%>%
  mutate(mean = mean(pcorr.value))%>%
  group_by(Celltype, Time)%>%
  mutate(low = ci(mean)[2])%>%
  mutate(high = ci(mean)[3])%>%
  mutate(mean = mean(mean))%>%
  distinct(Time, Celltype, mean, low, high)

plot$Celltype <- factor(plot$Celltype, levels = c("WT","M1M23"))

svg("fig21.svg", bg = "transparent", width = 6, height = 5)
ggplot(data = plot, aes(Time, mean, fill = Celltype))+
  geom_line()+
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.7)+
  scale_fill_manual(values = c("goldenrod1", "darkred"))+
  ggtitle("1321N1 WT vs M1M23 average response")+
  labs(caption = "ribbon represents 95% confidence interval")+
  theme_gray(base_size = 25)+
  transp
dev.off()



#### Define numerical columns ####
numvars_cell = c("Time","Mean_value",".fitted", "pcorr.value", "iplast_exp","iplast_rec","log_alpha", "rec_alpha", "min","max","mean_start","rec_baseline","rec_min","signal_rec")
numvars_cell_summary = c("iplast_exp","iplast_rec","log_alpha", "rec_alpha", "min","max","mean_start","rec_baseline","rec_min","signal_rec")
numvars_well = c("well_mean_Time", "well_mean_iplast_exp","well_mean_iplast_rec","well_mean_log_alpha", "well_mean_rec_alpha", "well_mean_min","well_mean_max","well_mean_mean_start","well_mean_rec_baseline","well_mean_rec_min","well_mean_signal_rec")
numvars_well_median_summary = c("well_median_iplast_exp","well_median_iplast_rec","well_median_log_alpha", "well_median_rec_alpha", "well_median_min","well_median_max","well_median_rec_baseline","well_median_rec_min","well_median_signal_rec")
numvars_well_median_summary_v2 = c("Area", "well_median_iplast_exp","well_median_iplast_rec","well_median_log_alpha", "well_median_rec_alpha", "well_median_min","well_median_max","well_median_rec_baseline","well_median_mean_start","well_median_rec_min","well_median_signal_rec")

numvars_well_summary = c("well_mean_iplast_exp","well_mean_iplast_rec","well_mean_log_alpha", "well_mean_rec_alpha", "well_mean_min","well_mean_max","well_mean_mean_start","well_mean_rec_baseline","well_mean_rec_min","well_mean_signal_rec")

idvars_cell = c("ROI","Plate","Well_in_name","Condition_name","Celltype") #Treatment
idvars_well = c("Plate","Well_in_name","Condition_name","Celltype") #Treatment

#### Barplots mean of well median HST ####
# Data in long format zetten
data.temp = comb.data %>% 
  ungroup() %>% 
  select(-ROI, -well_label, -unique, -Unique, -V1) %>% 
  group_by(Plate, Well_in_name, Condition_name, Celltype) %>% 
  summarise_all(median) %>% 
  distinct()
  
data.temp.long <- data.temp  %>%
  ungroup() %>%
  select(c("Celltype", "Concentration", numvars_well_median_summary_v3)) %>%
  pivot_longer(-c("Celltype", "Concentration"), names_to = "variables", values_to = "value")
data.temp.long$variables = sub("well_median_","",data.temp.long$variables)

# Data prep plots
plot = data.temp.long
plot$Celltype = factor(plot$Celltype,levels = c("WT","M1M23"))
plot$Concentration = factor(plot$Concentration,levels = c("PB", "C", "0", "1", "2","3","4"))
# Visualisatie parameters
col=c("Grey30","grey70")

transp =theme(
  plot.title = element_text(hjust = 0.5, size = 20),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background = element_blank(),
  legend.title = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))

# Bar plots
svg(paste0(plot.dir, "/CAM19_CTRL_variable_barplots.svg"), bg = "transparent")

ggplot(plot,aes(x = as.factor(Concentration),y = value, fill =  as.factor(Celltype))) +
  stat_summary(fun.data = mean_se , geom ="errorbar",position=position_dodge(width = 0.9, preserve = "single"),width=.5) +
  stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"),colour="white",width=0.9) +
  scale_fill_manual(values=c(col), guide=guide_legend(nrow=2))+
  xlab("Median well parameters")+
  ylab("")+
  facet_wrap(.~variables, scales = 'free')+ #,space = "free")+
  theme_minimal(base_size = 12)+
  theme(legend.position = "bottom", panel.spacing = unit(2,"lines"))+
  transp
dev.off()






#### Linear Mixed Model ####
# 1. Is LMM necessary and wanted?
# 2. Are the random effects stable over time (repeated sampling)? Test full (all fixed parameters) model with both and compare REML (lower = better)
#     if stable: intercept only: (1|parameter)
#     if not stable: slope: (1+Time|parameter)
#     Covariance of intercept and slope?
library(pbkrtest)
library(lme4)
library(multcomp)
LMM_well1 = comb.int.data %>% 
  ungroup() %>% 
  select(-ROI, -well_label, -unique, -Unique, -V1) %>% 
  group_by(Plate, Well_in_name, Condition_name, Celltype) %>% 
  summarise_all(median) %>% 
  distinct()

lmm1 = lmer(slope_exp~Celltype*Condition_name+(1|Plate),data = LMM_well1)
summary(lmm1)
lmm2 = lmer(slope_exp~Celltype*Condition_name+(1|Plate)+(1|Condition_name),data = LMM_well1)
summary(lmm2)
lmm3 = lmer(slope_exp~Celltype*Condition_name+(1|Plate),data = LMM_well1)
summary(lmm3)

plot
numvars_plot
dep_vars = "log_alpha"

lmm1 = lmer(slope_exp~Celltype+(Well_in_name|Plate),data = plot)

model = lmer(log(var.select+1)~Celltype+
               (Well_in_name|Plate), 
             data = plot, control =lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

drop1(model,test="Chisq") # p-value by likelihood ratio test

# If interaction significant: Make new var and rerun
data.temp$CombVar = factor(paste(data.temp$Condition,data.temp$,data.temp$Dilution,sep="_")) 
model = lmer(var.select~CombVar+(1|Rep/Well), data = data.temp, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glht(model, linfct = mcp(CombVar = "Tukey"))) # or Dunett 

comb.data$Condition = factor(comb.data$Condition)
lmm1 = lmer(Decay_rate~Condition*Celltype + (1|Plate), data = comb.data)
summary(glht(lmm1, linfct = mcp(Condition = "Tukey")))
?glht

data.temp = comb.data %>% filter(Celltype =="M1M23") %>% 
anova = aov(Decay_rate~Condition, data = data.temp)
TukeyHSD(anova)







#### ANOVA on well level  ####
numvars_well = c("Decay_rate","Recovery_rate","min", "max", "Baseline","Recovery")

data_temp = comb.data %>% 
  drop_na() %>% 
  group_by(Plate, Well_in_name, Condition) %>% 
  summarise_if(is.numeric, mean)
data_temp$Condition = factor(data_temp$Condition)

sink(paste0(plot.dir,"/CAM06_M1M23_Conditions_var_anova_wells.txt"))
print("Mean values per condition/plate:")
tmp1 = data_temp %>% group_by(Condition,Plate) %>% 
  summarise_if(is.numeric, mean) %>% select(-ROI)
print(tmp1)
print("Standard deviations:")
tmp2 = data_temp %>% group_by(Condition,Plate) %>% 
  summarise_if(is.numeric, sd) %>% select(-ROI)
print(tmp2)

for(i in numvars_well){
  print("===========================================================")
  print(i)
  print("-------------------------------")
  test = aov(get(i)~Condition, data_temp)
  Tukey = TukeyHSD(test)
  print(Tukey)
}
sink()

# LMM new
LMM_well1 = comb.int.data %>% 
  ungroup() %>% 
  select(-ROI, -well_label, -unique, -Unique, -V1) %>% 
  group_by(Plate, Well_in_name, Condition_name, Celltype) %>% 
  summarise_all(median) %>% 
  distinct()

lmm1 = lmer(slope_exp~Celltype*Condition_name+(1|Plate),data = LMM_well1)
summary(lmm1)








