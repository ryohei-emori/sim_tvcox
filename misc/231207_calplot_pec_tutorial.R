library(prodlim)
library(lava)
library(riskRegression)
library(survival)
library(pec)

#### survival
dlearn <- SimSurv(40)
dval <- SimSurv(100)
f <- coxph(Surv(time,status)~X1+X2,data=dlearn,x=TRUE,y=TRUE)
cf=calPlot(f,time=3,data=dval,type = "risk", add=TRUE)
print(cf)
plot(cf)
g <- coxph(Surv(time,status)~X2,data=dlearn,x=TRUE,y=TRUE)
cf2=calPlot(list("Cox regression X1+X2"=f,"Cox regression X2"=g),
            time=4,
            type="risk",
            data=dval)
print(cf2)

sample_list <- list()
sample_list[["5000"]] <- list()
sample_list[["5000"]]$est_model <- list()

sample_list[["3000"]]$est_model[[1]] <- coxph(Surv(time,status)~X1+X2,data=dlearn,x=TRUE,y=TRUE)
sample_list[["3000"]]$est_model[[2]] <- coxph(Surv(time,status)~X1+X2,data=dlearn,x=TRUE,y=TRUE)
sample_list[["5000"]]$est_model[[1]] <- coxph(Surv(time,status)~X2,data=dlearn,x=TRUE,y=TRUE)
sample_list[["5000"]]$est_model[[2]] <- coxph(Surv(time,status)~X1,data=dlearn,x=TRUE,y=TRUE)

rep_times <- 2
col1 = rgb(0.4, 0.6, 1, 0.9)  # 淡い青色の透過度50%
col2 = rgb(0.4, 1, 0.4, 0.9)  # 淡い緑色の透過度50%

cfg=calPlot(object = c(sample_list[["3000"]]$est_model,
                       sample_list[["5000"]]$est_model),
            time = 4,
            data = dval,
            type = "risk",
            pseudo = TRUE,
            showPseudo = TRUE,
            pseudo.col = "gray",
            bars = FALSE,
            add =FALSE,
            legend = FALSE,
            col = c(rep(col1, rep_times), rep(col2, rep_times)),
            lwd = c(rep(2, rep_times), rep(2, rep_times)),
            lty = c(rep("twodash", rep_times), rep("dashed", rep_times)),
            xlab = "Predicted risk probabilities",
            ylab = "Observed risk frequencies"
            )

plot(cf2)
calPlot(f,time=3,data=dval,type="survival")
#calPlot(f,time=3,data=dval,bars=TRUE,pseudo=FALSE)
#calPlot(f,time=3,data=dval,bars=TRUE,type="risk",pseudo=FALSE)
