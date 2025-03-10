# Load necessary library
library(ggplot2)
library(tidyr)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Function to generate softmax-transformed weights
softmax <- function(x) {
  exp_x <- exp(x)
  return(exp_x / sum(exp_x))
}

# Number of samples
n_samples <- 10000

# Generate samples from a normal distribution and apply softmax
samples <- replicate(n_samples, softmax(rnorm(3, 0, 1)))

# Convert to a data frame
df <- as.data.frame(t(samples))
colnames(df) <- c("X1", "X2", "X3")

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

# Plot density distributions
ggplot(df_long, aes(x = Value, fill = Variable)) +
  geom_density(alpha = 0.1) +
  labs(title = "Distribution of Softmax-Transformed Weights",
       x = "Weight",
       y = "Density") +
  theme_minimal()

df_long %>%
  group_by(Variable) %>%
  summarize(var(Value))


