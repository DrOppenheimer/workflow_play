my_data <- import_data("~/Documents/GitHub/Kevin_R_scripts/GSE198449/GSE198449.data.txt")
all_my_data <- as.vector(my_data)
hist(all_my_data, breaks =100)
#hist(log10(all_my_data), breaks =100) # is roughly log normal

# log transform
my_data_lt <- log10(my_data)
all_my_data_lt <- as.vector(my_data_lt)
hist(all_my_data_lt, breaks =100)

# center/standardize
#my_data_lt_scaled <- scale(my_data_lt) # This fails so use function below

all_my_data_lt_scaled <- as.vector(my_data_lt_scaled)
hist(all_my_data_lt_scaled, breaks =100)

















# with ggplot
all_my_data <- data.frame(Values = as.vector(my_data))
ggplot(data=all_my_data,
       mapping=aes(x=Values))+
  geom_histogram()


