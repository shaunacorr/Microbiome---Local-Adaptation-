library(tidyr)
library(growthcurver)
library(dplyr)
library(lme4)
library(tidyverse)
library (cowplot)
library(emmeans)
library(ggplot2)
library(FSA)
library(car)

##DD (Dulse isolates in Dulse Media)####
#Read in file DD
All_DD <- read.csv("All_DD_Plates.csv",header = T)
#Make Id COlumn
All_DD$variable <- paste(All_DD$true_rep, All_DD$isolate)
# transform into format needed for growth curver
DD_input <- pivot_wider(
  data = All_DD,
  id_cols = time,
  names_from = variable,
  values_from = value)
#zero all negative values
DD_input <- mutate_all(DD_input, function(x) ifelse(x < 0, 0, x))

##Growth curver won't work for NA values subset for loss of measurment of Dulse 2 Isolate 7 at 6 hrs

DD_input_sub1 <- DD_input[,-43]
DD_input_sub2 <- DD_input[-3, c(1, 43)]

#Complete growth curver analysis
gc_DD_sub1 <- SummarizeGrowthByPlate(DD_input_sub1, bg_correct = 'none')
gc_DD_sub2 <- SummarizeGrowthByPlate(DD_input_sub2, bg_correct = 'none')


##DF (Dulse Isolates in Fucus Media) #####
#Read in file DF
All_DF <- read.csv("All_DF_Plates.csv",header = T)
#Make Id COlumn
All_DF$variable <- paste(All_DF$true_rep, All_DF$isolate)
# transform into format needed for growth curver
DF_input <- pivot_wider(
  data = All_DF,
  id_cols = time,
  names_from = variable,
  values_from = value)
#zero all negative values
DF_input <- mutate_all(DF_input, function(x) ifelse(x < 0, 0, x))

##Growth curver won't work for NA values subset for loss of measurment of Dulse 3 Isolate 11 & 15 at 4 hrs and Dulse 2 Isolate 14 & 18 and Dulse 5 Isolate 24 at 12hrs

DF_input_sub1 <- DF_input[,c(-60,-64,-39,-35,-121)]
DF_input_sub2 <- DF_input[-2, c(1, 60,64)]
DF_input_sub3 <- DF_input[-6, c(1, 39,35,121)]

#Complete growth curver analysis
gc_DF_sub1 <- SummarizeGrowthByPlate(DF_input_sub1, bg_correct = 'none')
gc_DF_sub2 <- SummarizeGrowthByPlate(DF_input_sub2, bg_correct = 'none')
gc_DF_sub3 <- SummarizeGrowthByPlate(DF_input_sub3, bg_correct = 'none')

##FF (Fucus Isolates in Fucus Media) ####
#Read in file FF
All_FF <- read.csv("All_FF_Plates.csv",header = T)
#Make Id COlumn
All_FF$variable <- paste(All_FF$true_rep, All_FF$isolate)
# transform into format needed for growth curver
FF_input <- pivot_wider(
  data = All_FF,
  id_cols = time,
  names_from = variable,
  values_from = value)
# Growth curver won'ts work on negatives, zero all negative values

FF_input <- mutate_all(FF_input, function(x) ifelse(x < 0, 0, x))

#Growth curver won't work for NA values subset for loss of meaurement at 4 hrs for Fucus rep 2 Isolate 9,21,22 & 24 and Fucus rep 3 Isolate 2,3 ,5,6 & 7 and Fucus 4 Isolate 13, 2 & 9 and Fucus 5 isolate 15 & 22 at 6 hrs.
FF_input_sub1 <- FF_input[,c(-45:-47,-49,-52,-53,-56,-57,-78,-81,-88,-108,-119,-55)]
FF_input_sub2 <- FF_input[-2, c(1, 45:47, 49, 52, 53, 55:57)]
FF_input_sub3 <- FF_input[-3, c(1, 78, 81, 88, 108, 119)]
#Complete growth curver analysis
gc_FF_sub1 <- SummarizeGrowthByPlate(FF_input_sub1, bg_correct = 'none')
gc_FF_sub2 <- SummarizeGrowthByPlate(FF_input_sub2, bg_correct = 'none')
gc_FF_sub3 <- SummarizeGrowthByPlate(FF_input_sub3, bg_correct = 'none')

#Read in file FD (Fucus Isolates in Dulse Media)####
All_FD <- read.csv("All_FD_Plates.csv",header = T)
#Make Id COlumn
All_FD$variable <- paste(All_FD$true_rep, All_FD$isolate)
# transform into format needed for growth curver
FD_input <- pivot_wider(
  data = All_FD,
  id_cols = time,
  names_from = variable,
  values_from = value)
#zero all negative values
FD_input <- mutate_all(FD_input, function(x) ifelse(x < 0, 0, x))

#Growth curver won't work for NA values subset for loss of measurement at 4 hrs for Fucus 3 Isolate 1 & 7 and Fucus 5 Isolate 9 and 6 hrs for Fucus 5 Isolate 7 & 18 

FD_input_sub1 <- FD_input[, c(-50,-57,-115,-117,-111)]
FD_input_sub2 <- FD_input[-2, c(1, 50,57,117)]
FD_input_sub3 <- FD_input[-3, c(1, 115,111)]

#Complete growth curver analysis
gc_FD_sub1 <- SummarizeGrowthByPlate(FD_input_sub1, bg_correct = 'none')
gc_FD_sub2 <- SummarizeGrowthByPlate(FD_input_sub2, bg_correct = 'none')
gc_FD_sub3 <- SummarizeGrowthByPlate(FD_input_sub3, bg_correct = 'none')


################################################################
#Look at each growth curve in turn to check fit
# DD - Growth Curve Plot Check ####

# Get the names of the columns to loop through
column_names <- names(DD_input_sub1)

# Create an empty list to store the results
gc_fits <- list()

# Loop through each column and summarize/plot
for (col in column_names) {
  gc_fit <- SummarizeGrowth(DD_input$time, DD_input_sub1[[col]])
  
  # Plot the result
  plot(gc_fit)
  
  gc_fits[[col]] <- gc_fit
}

# Save the plots to separate PDF files

pdf("DD_input_sub1.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

#DD_input_sub2
column_names <- names(DD_input_sub2)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(DD_input_sub2$time, DD_input_sub2[[col]])
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("DD_input_sub2.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
   
  }
}
dev.off()

# DF - GC Plot Check ####

column_names <- names(DF_input_sub1)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(DF_input_sub1$time, DF_input_sub1[[col]])
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("DF_input_sub1.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# DF_input_sub2
column_names <- names(DF_input_sub2)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(DF_input_sub2$time, DF_input_sub2[[col]])
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("DF_input_sub2.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

#DF_input_sub3
column_names <- names(DF_input_sub3)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(DF_input_sub3$time, DF_input_sub3[[col]])
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("DF_input_sub3.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# FF - GC Plot Check ####

column_names <- names(FF_input_sub1)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(FF_input_sub1$time, FF_input_sub1[[col]])
  
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}


pdf("FF_input_sub1.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# FF_input_sub2
column_names <- names(FF_input_sub2)

gc_fits <- list()
for (col in column_names) {
  gc_fit <- SummarizeGrowth(FF_input_sub2$time, FF_input_sub2[[col]])
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("FF_input_sub2.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# FF_input_sub3
column_names <- names(FF_input_sub3)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(FF_input_sub3$time, FF_input_sub3[[col]])
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("FF_input_sub3.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# FD - GC Plot Check ####

column_names <- names(FD_input_sub1)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(FD_input_sub1$time, FD_input_sub1[[col]])
  
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("FD_input_sub1.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# FD Input Sub2
column_names <- names(FD_input_sub2)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(FD_input_sub2$time, FD_input_sub2[[col]])
  plot(gc_fit)
  gc_fits[[col]] <- gc_fit
}

pdf("FD_input_sub2.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# FD_input_sub3
column_names <- names(FD_input_sub3)
gc_fits <- list()

for (col in column_names) {
  gc_fit <- SummarizeGrowth(FD_input_sub3$time, FD_input_sub3[[col]])
  plot(gc_fit)
  
  gc_fits[[col]] <- gc_fit
}

pdf("FD_input_sub3.pdf")
for (col in column_names) {
  if (!is.null(gc_fits[[col]])) {
    plot(gc_fits[[col]], main = paste("Growth Curve for", col))
  }
}
dev.off()

# Some isolates did not grow and curves fitted are erroneous - those that did not grow set to 0. 
gc_DD_sub1[c(41:42,47,58,60,67,83,113,115),2:9] <- 0
gc_DF_sub1[c(25,27,42,30,50,58,63:64,114,129,137),2:9] <- 0
gc_DF_sub2[c(1:2),2:9] <- 0
gc_FD_sub1[c(3,28,44,31,55,56,61,66,69,82,84,85,91,93,95,110,100,124,117,127,128,131,135,138),2:9] <- 0
gc_FD_sub2[c(2),2:9] <- 0
gc_FF_sub1[c(25:27,29,42,32:33,36,38,51:52,55:56,61,63,75,91,94,113,129),2:9] <- 0

gc_DD <- bind_rows(gc_DD_sub1,gc_DD_sub2)
gc_DF <- bind_rows(gc_DF_sub1,gc_DF_sub2,gc_DF_sub3)
gc_FD <- bind_rows(gc_FD_sub1,gc_FD_sub2,gc_FD_sub3)
gc_FF <- bind_rows(gc_FF_sub1,gc_FF_sub2,gc_FF_sub3)

write.csv(gc_DD, file = "All_DD_AUC.csv", row.names = FALSE)
write.csv(gc_DF, file = "All_DF_AUC.csv", row.names = FALSE)
write.csv(gc_FD, file = "All_FD_AUC.csv", row.names = FALSE)
write.csv(gc_FF, file = "All_FF_AUC.csv", row.names = FALSE)

#Mixed Effects Model Analysis####

##read in the data set
mixed <- read.csv("interaction.csv", header = TRUE)

#Give the factors names
mixed$media<- factor(mixed$media,
                     levels = c("1","2"),
                     labels = c("dulse", "fucus"))

mixed$origin<- factor(mixed$origin,
                      levels = c(1,2),
                      labels = c("dulse", "fucus"))

mixed  <- mixed %>% mutate(
  isolate = as.factor(isolate),
  seaweed_rep = as.factor(seaweed_rep))

## remove instances where isolates have zero growth on either media
mixed_II <- mixed %>%
  filter(auc_e != 0)%>%
  mutate(complete = duplicated(isolate) | duplicated(isolate, fromLast = TRUE))%>%
  filter(complete == "TRUE")


## Calculate delta growth for each isolate ##

dulse <- subset(mixed_II, mixed_II$origin == 'dulse')
fit_dulse <- lmList(auc_e ~ media  |name, data = dulse)
fit_dulse <- coef(fit_dulse)
dulse <- unique(dulse[,c(2:3,5)])
dulse <- cbind(dulse,fit_dulse)

fucus <- subset(mixed_II, mixed_II$origin == 'fucus')
fit_fucus <- lmList(auc_e ~ media  |name, data = fucus)
fit_fucus <- coef(fit_fucus)
fucus <- unique(fucus[,c(2:3,5)])
fucus <- cbind(fucus,fit_fucus)

host_response <- rbind(dulse,fucus)
host_response <- host_response %>% rename(growth_response = mediafucus)

## linear mixed effects model ##
## Now we have a single predictor (origin) and we code seaweed replicate as a random effect
response_fit <- lmer(growth_response ~ origin + (1|seaweed_rep), host_response, REML=FALSE)
summary(response_fit)

## check that origin provides some predictive power
response_fit_II <- lmer(growth_response ~ 1 + (1|seaweed_rep), host_response, REML=FALSE)
summary(response_fit_II)
anova(response_fit, response_fit_II)

#investigate contrasts
emmeans(response_fit, pairwise~origin )

## variation in host plot ##

#Dulse
hist(dulse$mediafucus)
shapiro.test(dulse$mediafucus)
kruskal.test(mediafucus ~ seaweed_rep, data=dulse)
#Kruskal-Wallis chi-squared = 25.365, df = 5, p-value = 0.0001184
dunnTest(dulse$mediafucus, dulse$seaweed_rep, method="bh")

#Fucus

hist(fucus$mediafucus)
shapiro.test(fucus$mediafucus)
fucus_anova <- aov(mediafucus ~ seaweed_rep, data = fucus)
summary(fucus_anova)
# Check normality of residuals
shapiro.test(residuals(fucus_anova))

leveneTest(mediafucus ~ seaweed_rep, data = fucus)

# Plot a Q-Q plot
qqnorm(residuals(fucus_anova))
qqline(residuals(fucus_anova), col = "blue")
TukeyHSD(fucus_anova)




#### Plots ####

average_df2 <- mixed_II %>%
  group_by(seaweed_rep, origin, media) %>%
  summarise(mean_auc_e = mean(auc_e), .groups = "drop")


# Prepare summary statistics
summary_stats2 <- mixed_II %>%
  group_by(origin, media) %>%
  summarise(
    mean_measurement = mean(auc_e),
    sd_measurement= sd(auc_e),
    se_measurement = sd(auc_e) / sqrt(n()),
    n = n(),
    margin_of_error = qt(0.975, df = n - 1) * se_measurement,
    lower_bound = mean_measurement - margin_of_error,
    upper_bound = mean_measurement + margin_of_error,
    .groups = "drop"
  )

# Create a merged plot
p2 <- ggplot() +
  # Plotting the summary statistics with ribbons
  geom_point(data = summary_stats2, aes(x = media, y = mean_measurement, shape = origin, color = origin), size = 3) +
  geom_errorbar(data = summary_stats2, aes(x = media, ymin = mean_measurement - se_measurement, ymax = mean_measurement + se_measurement, color = origin), width = 0.1, size =1) +  # Error bars inherit color
  geom_line(data = summary_stats2, aes(x = media, y = mean_measurement, group = origin, color = origin), size = 2) +
  geom_ribbon(data = summary_stats2, aes(x = media, ymin = lower_bound, ymax = upper_bound, fill = origin), alpha = 0.3, colour = NA) +
  
  # Plotting the average AUC for seaweed data
  geom_point(data = average_df2, aes(x = media, y = mean_auc_e, shape = origin, color = origin), size = 2.5) +
  geom_line(data = average_df2, aes(x = media, y = mean_auc_e, group = seaweed_rep, color = origin), linetype = "dotdash", size = 0.5) +
  
  # Customizing the plot
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11.5),  # Set x-axis title size
        axis.title.y = element_text(size = 11.5),  # Set y-axis title size
        axis.text.x = element_text(face = "italic")) + 
  labs(y = "Fitness (Average AUC)", x = "Environment") +
  ylim(0, 1.4) +
  scale_x_discrete(expand = c(0.3, 0), labels = c("Palmaria palmata", "Fucus serratus")) +
  scale_color_manual(values = c("#660033", "#669966")) +  
  scale_fill_manual(values = c("#660033", "#669966"))
p2

response_plot2 <- ggplot(host_response) +
  geom_boxplot(aes(growth_response, seaweed_rep, fill = origin)) +
  geom_vline(xintercept = 0) +
  labs(x = expression(paste(italic("P. palmaria"), " Environment","                   ",italic("F. serratus"), " Environment")), 
       y = 'Host Seaweed Replicate') +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10),
        axis.text.x =  element_text(size = 8), 
        axis.title.y = element_text(size = 11.5)) +  # Adjust y-axis title size
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),  # Set custom breaks
                     labels = c("1", "0.5", "0", "0.5", "1")) +  # Custom labels
  scale_color_manual(values = c("#660033", "#669966")) +  
  scale_fill_manual(values = c("#660033", "#669966"))
response_plot2




p2 <- p2 + 
  annotate("text", x = -Inf, y = Inf, label = " a)", hjust = 0, vjust = 1, size = 5)

response_plot2 <- response_plot2 + 
  annotate("text", x = -Inf, y = Inf, label = " b)", hjust = 0, vjust = 1, size = 5)
grid.arrange(p2,response_plot2, nrow=1)
Fig4 <- arrangeGrob(p2, response_plot2, nrow=1)
png("Figure4.png", width = 3000, height = 1500, res = 300)
grid.draw(Fig4)  # Draw the grob to the file
dev.off()



