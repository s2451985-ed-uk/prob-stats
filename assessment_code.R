library(readxl)
library(writexl)
library(tidyverse)
library(janitor)

# import the biomarkers data
biomarkers <- read_excel("biomarkers.xlsx")

# clean the column names
biomarkers <- clean_names(biomarkers)

# import the covariates data
covariates <- read_excel("covariates.xlsx")

# clean the column names
covariates <- clean_names(covariates)

# compute the patient_id from the biomarker column
biomarkers <- 
  biomarkers %>%
  mutate(patient_id = str_split(biomarker, "-")) %>%
  mutate(patient_id = map(patient_id, first) %>% 
           reduce(c)) %>%
  mutate(patient_id = as.numeric(patient_id))

# compute the weeks from the biomarker column
biomarkers <- 
  biomarkers %>%
  mutate(weeks = str_split(biomarker, "-")) %>%
  mutate(weeks = map(weeks, last) %>% 
           reduce(c))

# merge biomarkers and covariates tables by the common patient ids
all_data <- biomarkers %>%
  full_join(covariates, by = "patient_id")

# obtaining the summary of all_data and covariates by column to check the null values after the merge
summary(all_data)
summary(covariates)
# patient 42 and 51 have null values in Vas-12months column

# checking if weeks column need any data cleaning
unique_values <- unique(all_data$weeks)
print(unique_values)
# 3 main column values as expected, no need to clean this column

# creating the inclusion_data table by selecting the data at the week 0
inclusion_data <- all_data %>%
  filter(weeks == "0weeks")

# checking whether any patient id has a lack of data at week 0 (inclusion data)
patient_unique_values_diff <- setdiff(union(unique(covariates$patient_id), unique(inclusion_data$patient_id)), 
                              intersect(unique(covariates$patient_id), unique(inclusion_data$patient_id)))
print(patient_unique_values_diff)
# patient 40 does not have an inclusion data

# selecting the inclusion data for male
inclusion_data_male <- inclusion_data %>%
  filter(sex_1_male_2_female == 1)

# selecting the inclusion data for female
inclusion_data_female <- inclusion_data %>%
  filter(sex_1_male_2_female == 2)

biomarker_names <- c("il_8", "vegf_a", "opg", "tgf_beta_1", "il_6", "cxcl9", 
                     "cxcl1", "il_18", "csf_1")

# generating the histogram of each biomarker data distribution for males at inclusion
inclusion_data_male %>%
  select(all_of(biomarker_names)) %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value)) + 
  geom_histogram(bins = 30, fill = "blue", color = "black") + 
  facet_wrap(~ Variable, scales = "free")

# generating the histogram of each biomarker data distribution for females at inclusion
inclusion_data_female %>%
  select(all_of(biomarker_names)) %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value)) + 
  geom_histogram(bins = 30, fill = "blue", color = "black") + 
  facet_wrap(~ Variable, scales = "free")

# creating an empty dataframe with two columns; biomarker_names and p_value
p_results_male_female_difference <- data.frame(
  biomarker_names = character(0),
  p_value = numeric(0),
  stringsAsFactors = FALSE 
)


#  loop through each biomarker to compute the p value through the t test and save the results in the dataframe
for (biomarker in biomarker_names) { 
  
  pval <- t.test(inclusion_data_male[[biomarker]], inclusion_data_female[[biomarker]])$p.value
  pval <- signif(pval, digits = 4)
  p_results_male_female_difference <- rbind(p_results_male_female_difference, 
                                            data.frame(biomarker_names = biomarker, p_value = pval))
  
}
# PRINT THE TABLE HERE
print(p_results_male_female_difference)

# Calculating the probability of making at least one type I error assuming that 
# the tests are independent and that all null hypotheses are true
1 - dbinom(0,9,0.05)

# Bonferroni correction
1 - dbinom(0,9,0.05/9)

# obtain the regression training data by selecting the independent and dependent variables
# turning the categorical data into 0 and 1 values and selecting 80% patients
regression_train_data <- inclusion_data %>%
  select("il_8", "vegf_a", "opg", "tgf_beta_1", "il_6", "cxcl9", "cxcl1", "il_18",
         "csf_1", "age", "sex_1_male_2_female", "smoker_1_yes_2_no", "vas_12months") %>%
  mutate(sex_1_male_2_female = if_else(sex_1_male_2_female == 2, 0, sex_1_male_2_female)) %>%
  rename(sex_1_male_0_female = sex_1_male_2_female) %>%
  mutate(smoker_1_yes_2_no = if_else(smoker_1_yes_2_no == 2, 0, smoker_1_yes_2_no)) %>%
  rename(smoker_1_yes_0_no = smoker_1_yes_2_no) %>%
  slice(1:94) 


# generate the scatterplots for each independent variable vs the dependent variable
regression_train_data %>% 
  pivot_longer(cols = biomarker_names, names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value, y = vas_12months)) + 
  geom_point(colour = "darkgreen") + 
  facet_wrap(~ Variable, scales = "free")
  
# obtain the independent variables for the remaining 20% of the data  
regression_test_data_x <- inclusion_data %>%
  select("il_8", "vegf_a", "opg", "tgf_beta_1", "il_6", "cxcl9", "cxcl1", "il_18",
         "csf_1", "age", "sex_1_male_2_female", "smoker_1_yes_2_no") %>%
  mutate(sex_1_male_2_female = if_else(sex_1_male_2_female == 2, 0, sex_1_male_2_female)) %>%
  rename(sex_1_male_0_female = sex_1_male_2_female) %>%
  mutate(smoker_1_yes_2_no = if_else(smoker_1_yes_2_no == 2, 0, smoker_1_yes_2_no)) %>%
  rename(smoker_1_yes_0_no = smoker_1_yes_2_no) %>%
  slice(95:117)

# obtain the dependent variable for the remaining 20% of the data  
regression_test_data_y <- inclusion_data %>%
  select("vas_12months") %>%
  slice(95:117)

# train the linear model with the regression_train_data
myLinear <- lm(vas_12months~il_8+vegf_a+opg+tgf_beta_1+il_6+cxcl9+cxcl1+il_18+
                 csf_1+age+sex_1_male_0_female+smoker_1_yes_0_no, data=regression_train_data)

# get the summary of the trained linear model
summary(myLinear)

# obtain the residuals of the trained model
training_residuals <- residuals(myLinear)

# plot the histogram of the residuals
hist(training_residuals)  

# plot the predictions versus the residuals
plot(predict(myLinear),training_residuals)

# make predictions for the remaining 20% of the data, test data
new_predictions <- predict(myLinear, newdata = regression_test_data_x)

# obtain the residuals of the test data
new_residuals <- regression_test_data_y$vas_12months - new_predictions

plot(new_predictions, regression_test_data_y$vas_12months)
cor(new_predictions,regression_test_data_y$vas_12months)

# compute the rmse for the test data
rmse_test <- sqrt(mean(new_residuals^2))

# compute the rmse for the train data
rmse_train <- sqrt(mean(training_residuals^2))
