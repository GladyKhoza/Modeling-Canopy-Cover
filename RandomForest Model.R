################ Empty Varibales
RMSE_train.rf <-numeric(length=2)
RMSE_test.rf <-numeric(length=2)
R.sq_train.rf <- numeric(length=2)
R.sq_test.rf <- numeric(length=2)
rRMSE_test.rf <- numeric(length=2)
multi_temp_master <- NULL
results_final <- NULL
rf.Impor <- matrix() ##Empty Matrix for RF importance
All_band_names <- NULL

## Importing as individual rasters 
CC<- raster("Canopy_Cover.tif")
names(CC) <- c("CC")
NAvalue(CC) <- -9999

##Aggregating Lidar CC
Canopy_Cover <- aggregate(CC,fact=4,fun=mean, filename=paste0(CC,"_100m.img"),overwrite=TRUE)

##Sampling points
sample.pnts <- rasterToPoints(CC, function(x)x > 0 & x <= 100,sp=T)

##Modelling RandomForest

##select the SAR images to include in the model
Orig.SAR.list <- list.files(path=getwd(),pattern =".tif$", full.names=TRUE)
SAR.bands.varlist <- select.list(Orig.SAR.list,multiple=TRUE,graphics=TRUE,title="SAR Images to use") #goes to a character list


n_bands <- length(SAR.bands.varlist)*3
for(p in 1:length(SAR.bands.varlist)){
  
  ## Subset the list of files
  img <- SAR.bands.varlist[p]  
  ## Read in the image 
  SAR.image <- stack(img)
  Jan2017 <- names(SAR.image)
  SAR_names <- substr(Jan2017,1,9)
  SAR_names <- c(paste0(SAR_names[1],"_VV"),paste0(SAR_names[2],"_VH"),paste0(SAR_names[3],"_NDI"))
  names(SAR.image) <- SAR_names
  SAR_band <- substr(SAR_names[2],1,9)
  All_band_names <- cbind(All_band_names,SAR_names) 
  ## To process only one band, specify it only

  ## Aggregate the 20m original SAR images to same resolution as the LiDAR data (i.e. 100m)
  SAR100 <- aggregate(SAR.image,fact=5,fun=mean, filename=paste0(SAR_band,"_100m.img"),overwrite=TRUE)
  names(SAR100) <- SAR_names
  
  ######RF MTRY (=Sqrt(number of variables)###################################
  var.num <- length(SAR_names)
  var.num <- (var.num)-1
  RF.mtry <- round(sqrt(var.num),0)
  
  ## Use the sample points to extract the SAR values 
  master.data <- raster::extract(SAR100,sample.pnts,sp=T)
  master.data.df <- data.frame(master.data@data)
  multi_temp_master <- cbind(multi_temp_master,master.data.df[,2])
  multi_temp_master <- cbind(multi_temp_master,master.data.df[,3])
  multi_temp_master <- cbind(multi_temp_master,master.data.df[,4])
  
  ################Number of K-fold cross validation runs
  a <- 2
  
  for(i in 1:a){
    
    ## Randomly sample 70% of the dataset
    n.sample <- sample(1:nrow(master.data.df), nrow(master.data.df)*0.7, replace=F) 
    ## Subset Training dataset
    training <- master.data.df[n.sample,]
    colnames(training) <- c("CC",SAR_names)
    ## Subset test dataset
    test <- master.data.df[-n.sample,]
    colnames(test) <- c("CC",SAR_names)
    
    ## Set aside the observed test variable
    test.obs <- data.frame(Obs=test$CC)
    
    
    ## Training the Random Forest model
    set.seed(1)
    CC.rf.tr <- randomForest(CC ~ .,
                             data=training, 
                             ntree=500,
                             mtry=RF.mtry,
                             na.action = na.omit,
                             importance=TRUE,
                             replace=FALSE)
    
    rf.Impor [i] <- list(round(importance(CC.rf.tr), 2))
    
    ##Predict for the models, using the test data 
    predCC.rf <- predict(CC.rf.tr,test)
    predCC.rf.df = data.frame(predCC.rf)
    
    ## Create the Predicted and Observed dataset for validation
    predVobs = data.frame(Predrf=predCC.rf.df[,1],Obs=test.obs)
    predVobs.rf = lm(Predrf~Obs,data=predVobs)
    
    n_obs <- nrow(test)
    
    ## Store the validation stats for each iteration
    R.sq_test.rf[i] <- summary(predVobs.rf)$r.squared
    RMSE_test.rf[i] <- sqrt(mean((predVobs$Obs-predVobs$Predrf)^2,na.rm=T))
    rRMSE_test.rf[i] <- sqrt(mean((predVobs$Obs-predVobs$Predrf)^2,na.rm=T))/(mean(predVobs$Obs,na.rm=T))*100
    
    
  }
  
  ## Aggregate the stored validation stats for final results 
  results <- data.frame(SAR_band,
                        R_Sq=mean(R.sq_test.rf),
                        RMSE=mean(RMSE_test.rf),
                        rRMSE=mean(rRMSE_test.rf))
  
  results_final <- rbind(results_final,results)
  View(results_final)
  
  band_names_melt <- melt(All_band_names)
  band_names_melt <- as.character(band_names_melt[,3])
  colnames(multi_temp_master) <- band_names_melt
  
  ##Apply the random forest model to the SAR image and map the vegetation variable
  predict(SAR100, model=CC.rf.tr, filename = paste0(SAR_band,"_RF_Model.tif"), fun=predict, progress="window",overwrite=TRUE)
  
  
}