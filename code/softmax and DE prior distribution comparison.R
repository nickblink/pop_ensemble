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
softmax_samples <- replicate(n_samples, softmax(rnorm(3, 0, 1)))

# Convert to a data frame
df_softmax <- as.data.frame(t(softmax_samples))
colnames(df_softmax) <- c("u1", "u2", "u3")

DE_samples <- replicate(n_samples,rnorm(3, 1/3, sd = .1))
df_DE <- as.data.frame(t(DE_samples))

X_samples <- replicate(n_samples, rnorm(3, 100, sd = 10))
X_df <- as.data.frame(t(X_samples))

Y_softmax <- rowSums(df_softmax * X_df)
Y_DE <- rowSums(df_DE * X_df)

sd(Y_softmax)
sd(Y_DE)

## Plot softmax
{
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
}

