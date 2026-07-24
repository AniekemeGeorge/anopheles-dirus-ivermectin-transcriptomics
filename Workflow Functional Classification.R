# Set working directory
setwd("C:/Users/19198826/Downloads/RotationProject2/Sequencing/DSEQ Data")
list.files()
library(ggplot2)
library(dplyr)
library(forcats)
df <- read_excel("Publication_Annotation_Table_ Alive_Tolerance.xlsx")
 
names(df)
category_count <- df %>%
  filter(!is.na(`Functional category`)) %>%
  count(`Functional category`, sort = TRUE)
plot_df <- df %>%
  filter(`Functional category` != "Other / Putative") %>%
  count(`Functional category`, sort = TRUE)
category_count

ggplot(plot_df,
       aes(x = forcats::fct_reorder(`Functional category`, n),
           y = n,
           fill = `Functional category`)) +
  
  geom_col(width = 0.8) +
  
  coord_flip() + 
  labs(
    title = "Functional classification of the top 40 DEGs in Alive mosquitoes",
    x = "",
    y = "Number of mapped orthologs"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none"
  )

#Functional Classification for Paralysed
 
df <- read_excel("Publication_Annotation_Table_Paralysed.xlsx")

names(df)
category_count <- df %>%
  filter(!is.na(`Functional category`)) %>%
  count(`Functional category`, sort = TRUE)
plot_df <- df %>%
  filter(`Functional category` != "Other / Putative") %>%
  count(`Functional category`, sort = TRUE)
category_count

ggplot(plot_df,
       aes(x = forcats::fct_reorder(`Functional category`, n),
           y = n,
           fill = `Functional category`)) +
  
  geom_col(width = 0.8) +
  
  coord_flip() + 
  labs(
    title = "Functional classification of the top 40 DEGs in paralysed mosquitoes",
    x = "",
    y = "Number of mapped orthologs"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none"
  )
ggsave("Paralysed_Functional_Classification.png",
       width = 7,
       height = 3.8,
       dpi = 600)
