PCR <- seq(0.01, 0.99, 0.01)
logit_LM <- logit(PCR)

plot(logit_LM, logit(PCR), type = "l")
lines(logit_LM, logit(PCR), col = "red")

expit(logit_LM - logit(PCR))

1/0.377

LM <- expit(logit_LM)

plot(PCR, LM, xlim = c(0, 1), ylim = c(0, 1), type = "l")


logit(0.9)

expit(-0.5)


PCR <- seq(0.01, 0.99, 0.01) #P2
delta_dash <- -0.5 #delta_dash

term <- (PCR/(1-PCR)) * exp(delta_dash)

LM <- term/(1 + term)

plot(PCR, LM, type = "l")
plot(PCR, LM/PCR, type = "l", xlim = c(0,1), ylim = c(0, 1))

logit_LM <- logit(PCR)



