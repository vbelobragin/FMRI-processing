#---this is a main file to work with
#---Neuroconductor in Windows 10
#---ALL functions of fslr package are in Ubuntu script file
#---named NeuroUb.R
#---so this file code is used to process all *.nii files
#---with functions that are accessible in W10, thus minimize work in Ubuntu
#-------------------------------------------------------------

#library(dslabs)
#library(plyr)
#detach(package:plyr)
#library (dplyr)
#library(SnowballC)
#library(ggmap)
#library(maps)
#library (mice)
#library(zoo)
library (tidyverse)
library(caTools)
library(ROCR)
library(rpart)
library(rpart.plot)
library(randomForest)
library(caret)
library(gam)
library(e1071)
library(tm)
library(devtools)
library(flexclust)
library(igraph)
library(rgl)
library(lubridate)
library(matrixStats)
library(oro.dicom)
library(oro.nifti)
library(fslr)
library(scales)
library(reshape2)

nhack.dir=getwd()
neuroData.dir="E:/Neuro_data/Neurocoductor_data/NeuroH_MOOC_files"
neuroTemplate.dir="E:/Neuro_data/Neurocoductor_data/Template"
#setwd(dir=neuroData.dir)
#setwd(dir=neuroTemplate.dir)
#setwd(dir=nhack.dir)

#---load as001-0005-00001-000001-01.nii---
nhack.dir=getwd()
neuroSkFunc.dir="E:/Neuro_data/test_data_Skolkovo/functional/P_[001]/"
setwd(dir=neuroSkFunc.dir)
fname1="as001-0005-00001-000001-01.nii"
fslnii1=readNIfTI(fname1, reorient=F)
fname2="as001-0005-00001-001001-01.nii"
fslnii2=readNIfTI(fname2, reorient=F)
double_ortho(fslnii1, fslnii2)

InFiles=c("as001-0005-00001-000001-01.nii", "as001-0005-00001-001001-01.nii")
fslmerge(infiles=InFiles, direction="t", outfile="probfunc", retimg=F)


#---load 113-01-MPRAGE.nii.gz---
pathname="E:/Neuro_data/Neurohacking_data/kirby21/visit_1/113/"
fname="113-01-MPRAGE.nii.gz"
niifile=paste0(pathname,fname)
fslnii=readNIfTI(niifile, reorient=F)

#---several simple checks---
print({fslnii})
dd12=dim(fslnii)
image(1:dd12[1], 1:dd12[2],fslnii[,,100],  xlab='',ylab='')
mean(fslnii)

#---obtain a file "fast_img" after function "fsl_biascorrect"---
niifile="tempnii1.nii.gz"
fast_img=readNIfTI(niifile, reorient=F)
print({fast_img}) #---simple check

#---l 3.4.---
sub.bias<-niftiarr(fslnii,fslnii-fast_img)
q=quantile(sub.bias[sub.bias!=0], probs=seq(0,1,by=0.1))
fcol=div_gradient_pal(low="blue", mid="yellow", high='red')
ortho2(fslnii, sub.bias, col.y = alpha(fcol(seq(0,1,length=10)),0.5),
      ybreaks=q, ycolorbar=T, text=paste0("Original Image Minus N4", "n\ Bias-Corrected Image"))
jpeg('rplot.jpg', width=1100, height=1100)
ortho2(fslnii, sub.bias, col.y = alpha(fcol(seq(0,1,length=10)),0.5),
       ybreaks=q, ycolorbar=T, text=paste0("Original Image Minus N4", "n\ Bias-Corrected Image"))
dev.off()
slices=c(2,6,10,14,18)
vals=lapply(slices, function(x){
  cbind(img=c(fslnii[,,x]), fast=c(fast_img[,,x]),
  slice=x)
})
vals=do.call("rbind",vals)
vals=data.frame(vals)
vals=vals[vals$img>0.8 & vals$fast>0,]
colnames(vals)[1:2]=c("Original Value", "Bias-Corrected Value")
vv=melt(vals, id.vars="slice")
gg=ggplot(aes(x=value, colour=factor(slice)), data=vv)+geom_line(stat="density")+facet_wrap(~variable)
gg=gg+scale_colour_discrete(name="Slice #")
gg

#---l 3.5.---
#---make a mask from a skull-stripped image file---
niifile="tempnii2.nii.gz"
bet_fast=readNIfTI(niifile, reorient=F)
print({bet_fast}) #---simple check
bet_fast_mask<-niftiarr(bet_fast, 1)
is_in_mask=bet_fast>0
bet_fast_mask[!is_in_mask]<-NA
orthographic(bet_fast)
orthographic(fast_img, bet_fast_mask)
#---

#---l 3.6.---
#---load template un-skulled MNI152_T1_1mm_brain.nii.gz---
pathname="Template/"
fname="MNI152_T1_1mm_brain.nii.gz"
niifile=paste0(pathname,fname)
nhack.dir=getwd()
neuroData.dir="E:/Neuro_data/Neurocoductor_data/"
setwd(dir=neuroData.dir)
template.nii=readNIfTI(niifile, reorient=F)
setwd(dir=nhack.dir)
#---several simple checks---
print({template.nii})
dd12=dim(template.nii)
image(1:dd12[1], 1:dd12[2],template.nii[,,100],  xlab='',ylab='')
mean(template.nii)
orthographic(template.nii)
#-------------------------
#---L 4.9.
atlas = "JHU_MNI_SS_WMPM_Type-I"
txtfile = paste0(atlas, "_SlicerLUT.txt")
### read look up table (LUT)
jhut1.df = read.table(txtfile, stringsAsFactors=FALSE)
jhut1.df = jhut1.df[, 1:2]
colnames(jhut1.df) = c("index", "Label")
jhut1.df$index = as.numeric(jhut1.df$index)
jhut1.df[1:4,]
jhut1.img = readNIfTI(paste0(atlas, ".nii.gz"))
uimg = sort(unique(c(jhut1.img)))
all.ind = jhut1.df$index
stopifnot(all(uimg %in% all.ind))
#--------
##Make a data frame with the index of the atlas
##and the value of the ROI at that voxel
roi.df = data.frame(index = jhut1.img[syn_roi> 0],
                    roi = syn_roi[ syn_roi > 0])
##Obtain the number (sum) of voxels that have an roi
##value >0.5 in the roi by the index of labels
label_sums = ddply(roi.df, .(index), summarize,
                   sum_roi = sum(roi), sum_roi_thresh = sum(roi > 0.5))
label_sums = merge(label_sums, jhut1.df, by="index")




#-------------------------------------
#---do NOT try this code - it will never run
# if (!require(devtools)){install.packages(devtools)}
# devtools::install_github("stnava/cmaker")
# devtools::install_github("stnava/ITKR")
# devtools::install_github("stnava/ANTsR")