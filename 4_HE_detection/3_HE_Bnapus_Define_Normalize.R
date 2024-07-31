library(tidyverse)
library(patchwork)
source("src/MLplot.r")
library(circlize)

### 1. read the accession id file and raw counts files from step2
napusA <- read_tsv("BnapusA_reciprocalMatch.cov",col_names = F)
napusC <- read_tsv("BnapusC_reciprocalMatch.cov",col_names = F)
id <- read_tsv("id.he.bnapus.acc",col_names = F)

df_napusA <- napusA %>% select(1,2,6,3,4,5,13:ncol(napusA))
names(df_napusA) <- c("linkC","linkA","geneA","chrA","startA","endA",str_c(id$X1,"A"))
df_napusC <- napusC %>% select(1,2,6,3,4,5,13:ncol(napusC))
names(df_napusC) <- c("linkC","linkA","geneC","chrC","startC","endC",str_c(id$X1,"C"))
# merge and select
df_large <- df_napusC %>% inner_join(df_napusA,by=c("linkC","linkA"))
df_large_select <- df_large %>% select(geneC,chrC,startC,endC,geneA,chrA,startA,endA,everything(),-linkC,-linkA) 
df_large_name <-  df_large %>% select(geneC,chrC,startC,endC,geneA,chrA,startA,endA)




### 2. Normalize gene length
df_large_select_normal_acc <- df_large_select %>%  mutate(lenA = endA - startA, lenC = endC -startC) %>% 
  select(9:last_col()) %>% 
  mutate(across( ends_with("C"), ~(.)/lenC ), across(ends_with("A"), ~(.)/lenA )) %>% select(-lenA, -lenC) 

id_depth <- df_large_select_normal_acc %>% summarise(across( everything(), mean )) %>% 
  pivot_longer(cols = everything(), names_to = c("acc","AC"),names_sep = -1, values_to = "depth"  ) %>% 
  arrange(acc) %>% 
  group_by(acc) %>% 
  summarise(depth_mean = mean(depth))

df_large_select_new <- bind_cols(df_large_name,df_large_select_normal_acc)
id_new <-  id %>% left_join(id_depth, c("X1" = "acc")) 




### 3. Filter the noise caused by local gene duplications or gene loss
#### Use threshold 1.5 && 0.375
output <- vector("list",nrow(id_new))
for (i in seq_along(id_new[[1]])) {
  i_name <- id_new[[i,1]]
  i_name_int <- str_c(i_name, "int")
  i_depth <- id_new[[i,2]]
  i_a <- str_c(i_name,"A")
  i_c <- str_c(i_name,"C")
  
  output[[i]] <- 
    df_large_select_new %>% select(matches(i_name)) %>% 
    mutate( !!i_name_int :=log2(.data[[i_a]]/.data[[i_c]])) %>% 
    mutate( !!i_name := 
              if_else( abs(.data[[i_name_int]]) >1.5 ,
                       if_else( (.data[[i_a]] < 0.375*i_depth & .data[[i_c]] >1.5*i_depth ) | (.data[[i_a]] >1.5*i_depth  & .data[[i_c]] < 0.375*i_depth  ), .data[[i_name_int]],0 ) ,
                       .data[[i_name_int]])
    ) %>% 
    select(!!i_name) 
}
df_tidy_log2 <- bind_cols(df_large_name,output)

### Get the log transformed HE value for each homoeologous gene pairs
write_csv(df_tidy_log2,file = "df_tidy_log2_new_202403_he15_dep0375_15.csv")
