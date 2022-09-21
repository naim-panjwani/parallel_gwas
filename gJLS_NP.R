gL_test <-function(model, data, correlation=NULL){
  
  ## check if there is missing data
  if(sum(is.na(data)) > 0)  stop("missing value(s) not allowed")
  
  ## Obtain the p-value scale test
  if(is.null(correlation)) {
    model=terms(model)
    #model0=model[-(1:2)]  # modify here so only the g term is removed
    model0=model[-(1)]     # line added by me (Naim)
    fit0<-lm(model0, data = data)
    fit<-lm(model,data=data)
    aovfit <- anova(fit0,fit)
    gL_F<-aovfit[2,5];numDF<-aovfit[2,3];denDF<-aovfit[2,1];gL_p<-aovfit[2,6]
  } else {
    fit<-gls(model,data=data,correlation=correlation,method="ML",control=lmeControl(opt = "optim"))
    aovfit<-anova(fit,Terms=2:dim(anova(fit))[1])
    gL_F<-aovfit[1,2];numDF<-aovfit[1,1];denDF<-fit$dims$N-fit$dims$p;gL_p<-aovfit[1,3]
  }
  ## return the scale p-value
  return(cbind(gL_F, numDF, denDF, gL_p))
}


#Note that p_S is obtained using a stage 1 median regression (rq function, tau=0.5) where group medians are chosen to be the larger of two middle values when the group size is even.
gS_test <-function(model, data, correlation=NULL){
  
  ## check if there is missing data
  if(sum(is.na(data)) > 0)  stop("missing value(s) not allowed")
  model <- terms(model)
  
  ## Obtain the p-value scale test
  lm1 <- rq(model,tau=.5,data=data)
  data$d1 <- abs(resid(lm1))
  model2=reformulate(attr(model, "term.labels"),response="d1")
  
  if(is.null(correlation)) {
    model2=terms(model2)
    #model0=model2[-(1:2)]  # modify here so only the g term is removed
    model0=model2[-(1)]     # line added by me (Naim)
    fit0<-lm(model0, data = data)
    fit<-lm(model2,data=data)
    aovfit<-anova(fit0,fit)
    gS_F<-aovfit[2,5];numDF<-aovfit[2,3];denDF<-aovfit[2,1];gS_p<-aovfit[2,6]
  } else {
    fit<-gls(model2,data=data,correlation= correlation,method="ML",control=lmeControl(opt = "optim"))
    aovfit<-anova(fit,Terms=2:dim(anova(fit))[1])
    gS_F<-aovfit[1,2];numDF<-aovfit[1,1];denDF<-fit$dims$N-fit$dims$p;gS_p<-aovfit[1,3]
  }
  ## return the scale p-value
  return(cbind(gS_F, numDF, denDF, gS_p))
}


gJLS_test <-function(model.loc,model.scale,data,correlation=NULL){
  
  ## check if there is missing data
  if(sum(is.na(data)) > 0)  stop("missing value(s) not allowed")
  
  
  ## Obtain the results from the individual location and scale test
  gL <- gL_test(model=model.loc,data=data,correlation=correlation)
  gS <- gS_test(model=model.scale,data=data,correlation=correlation)
  
  ## JLS test statistic (and corresponding p-value) using Fisher's combined p-value method
  t_JLS <- -2*log(gL[4])-2*log(gS[4])
  p_JLS <- 1-pchisq(t_JLS,4)
  
  ## return the location, scale and JLS p-values
  dataf=data.frame(rbind(gL,gS,c(t_JLS,4,"NA",p_JLS)))
  names(dataf)=c("Statistic","df1","df2","p-value")
  rownames(dataf)=c("gL","gS","gJLS")
  return(dataf)
}

