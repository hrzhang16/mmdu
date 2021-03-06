## code to prepare `DATASET` dataset goes here

library("stringi")
library("stringr")

us_senate = read.table("data-raw/us_senate108.txt",sep = ",")
data = as.character(us_senate$V1)
data = data[-1]

comm = vector()
N = length(data)
for(i in 1:N){
  comm = c(comm,substr(data[[i]],21,21))
}
comm = strtoi(comm)
comm_col = rep(0,N)
comm_col[which(comm==1)] = 4
comm_col[which(comm==2)] = 2
comm_col[which((comm!=1)&(comm!=2))]=3

name = rep(0,100);for(i in 1:100){name[i] = substr(data[[i]],26,36)}
state = rep(0,100);for(i in 1:100){state[i] = substr(data[[i]],13,20)}
abbr_state = as.character(read.csv("data-raw/us_states.csv",header = FALSE)[,3])
for(i in 1:100){state[i] = abbr_state[ceiling(i/2)]}

data = strsplit(data, " ")
reformat = list()
for(i in 1:N){
  ind = which(lapply(data[[i]],nchar)>200)
  reformat = c(reformat,list(data[[i]][ind]))
  temp = strtoi(strsplit(reformat[[i]],"")[[1]])
  reformat[[i]] = temp
}

a = unlist(lapply(reformat,length))
J = as.numeric(names(which.max(table(a))))

for(i in 1:N){
  reformat[[i]] = tail(reformat[[i]],J)
}


response = matrix(unlist(reformat), ncol = J, byrow = TRUE)
Omega = matrix(1,N,J)
Omega[(response==0) | (response==7) | (response==8) | (response==9)] = 0

response[ (response>=1) & (response<=3) ] = 1
response[ (response>=4) & (response<=6) ] = 0
response[Omega == 0] <- NA

comm_word = rep(0,100)
comm_word[comm==1] = "D"
comm_word[comm==2] = "R"
comm_word[comm==3] = "I"
info = rep(0,100);for(i in 1:100){info[i] = str_c(c(comm_word[i],state[i]),collapse = "-")}
senator = data.frame(name,info)
names(senator) = c("Name","State")

# bill <- read.csv("data-raw/senate_bill108.csv", sep = "")
# View(bill)

vote108 <- list("data" = response,
                "senate" = name,
                "party" = comm_word,
                "state" = state)


usethis::use_data(vote108, overwrite = TRUE)
