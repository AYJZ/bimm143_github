install.packages("tidyverse")
library(tidyverse)
source("http://thegrantlab.org/misc/cdc.R")
ggplot(data = cdc, mapping=aes(x = height, y = weight)) +
  geom_point(shape = 1, size = 3)
cor(cdc$height, cdc$weight)
hist(cdc$weight)
hist(cdc$height)

# Create height.m
height_m <- cdc$height * 0.0254

# Create weight conversion
weight_kg <- cdc$weight * 0.454

#bmi
bmi <- (weight_kg)/(height_m^2)
plot(cdc$height, bmi)

#obese individuals number
sum(bmi >= 30)

#percentage of obese individuals
round( (sum(bmi >= 30)/length(bmi)) * 100, 1)

first_hun_height <- cdc[1:100, "height"]
first_hun_weight <- cdc[1:100, "weight"]
ggplot(mapping = aes(x = first_hun_weight, y = first_hun_height)) +
  geom_point()

#obese male
new_cdc <- cbind(cdc, bmi)
sum((new_cdc$gender == "m") & (new_cdc$bmi >= 30))

#another way
obese_gender <- new_cdc$gender[new_cdc$bmi >= 30]
table(obese_gender)
