# Load database of TEE data for great apes
setwd("C:/Users/tskra/Dropbox/energetics_tsimane_hadza/manuscript/submission2_science/revise_and_resubmit")
db <- read.csv("./published data/uploaded/great_ape_TEE_DLW.csv", sep=",", header=T)   

## Chimpanzees
mod1 <- lm(log(TEE_kcal)~log(Mass_kg), data=db[which(db$Genus == "Pan"),])
summary(mod1)

# generate estimate and 95% CI at published adult body mass for males/females
chimp.est <- predict(mod1, newdata = data.frame(Mass_kg = c(41.7, 34.4)), se.fit = T, interval="confidence", level=0.95)
chimp.est <- as.data.frame(chimp.est)
chimp.est$fit.fit <- exp(chimp.est$fit.fit)
chimp.est$fit.lwr <- exp(chimp.est$fit.lwr)
chimp.est$fit.upr <- exp(chimp.est$fit.upr)
chimp.est$species <- "Chimpanzee"
chimp.est$sex <- c("M", "F")


## Gorillas
mod2 <- lm(log(TEE_kcal)~log(Mass_kg), data=db[which(db$Genus == "Gorilla"),])
summary(mod2)

# generate estimate and 95% CI at published adult body mass for males/females
gorilla.est <- predict(mod2, newdata = data.frame(Mass_kg = c(170.4, 70.5)), se.fit = T, interval="confidence", level=0.95)
gorilla.est <- as.data.frame(gorilla.est)
gorilla.est$fit.fit <- exp(gorilla.est$fit.fit)
gorilla.est$fit.lwr <- exp(gorilla.est$fit.lwr)
gorilla.est$fit.upr <- exp(gorilla.est$fit.upr)
gorilla.est$species <- "Gorilla"
gorilla.est$sex <- c("M", "F")

## Orangs
mod3 <- lm(log(TEE_kcal)~log(Mass_kg), data=db[which(db$Genus == "Pongo"),])
summary(mod3)

# generate estimate and 95% CI at published adult body mass for males/females
orang.est <- predict(mod3, newdata = data.frame(Mass_kg = c(78.5, 35.8)), se.fit = T, interval="confidence", level=0.95)
orang.est <- as.data.frame(orang.est)
orang.est$fit.fit <- exp(orang.est$fit.fit)
orang.est$fit.lwr <- exp(orang.est$fit.lwr)
orang.est$fit.upr <- exp(orang.est$fit.upr)
orang.est$species <- "Orangutan"
orang.est$sex <- c("M", "F")


# Combine
rbind(chimp.est, gorilla.est, orang.est)
