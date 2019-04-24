#---this is a processing of 3 test files from Skolkovo:
# P [001] Series 006 [MR - 3D MP RAGE]
# P [002] Series 006 [MR - 3D MP RAGE]
# P [003] Series 006 [MR - 3D MP RAGE]
#---next procedures implemented:
#---1) ortoimages
#---2)"heat maps" of 3 structural MRI files from Skolkovo in orderto check if they must be bias-corrected
#---3) skull-stripping
#---
#---the Ubuntu pair is in the file "Skolkovo1_Ub.R"
#---all images are kept in "E:\Neuro_data\test_data_Skolkovo\pictures_test_Skolkovo"
#---all processed files are kept in @E:\Neuro_data\test_data_Skolkovo\processed_structural_Skolkovo"
#-------------------------------------------------------------------

library (tidyverse)
library(oro.nifti)
library(fslr)
library(scales)
library(reshape2)

#---load 113-01-MPRAGE.nii.gz---
nhack.dir=getwd()
neuroData1.dir="E:/Neuro_data/test_data_Skolkovo/structural"
setwd(dir=neuroData1.dir)

#---this is for test purposes and may be removed---
nii=readNIfTI("ROI.nii.gz", reorient=F)
dim(nii)
orthographic(nii, xyz=c(200, 144, 11))
#--------------------------------------------------

#---orthoimages---------------------------------------------
#---1st file---
fname="P [001] Series 006 [MR - 3D MP RAGE].nii"
nii=readNIfTI(fname, reorient=F)
print({nii})
mean(nii)
dd12=dim(nii)
image(1:dd12[1], 1:dd12[2],nii[,,211],  xlab='',ylab='')
orthographic(nii, xyz=c(200, 220, 211))
#---2nd file---
fname="P [002] Series 006 [MR - 3D MP RAGE].nii"
nii=readNIfTI(fname, reorient=F)
orthographic(nii, xyz=c(200, 220, 211))
#---3rd file---
orthographic(nii, xyz=c(200, 220, 211))
fname="P [003] Series 006 [MR - 3D MP RAGE].nii"
nii=readNIfTI(fname, reorient=F)
orthographic(nii, xyz=c(200, 220, 211))
#------------------------------------------------------------
rm(nii)
gc()

#---heatmaps of corrected images-----------------------------
#---1st image---
nii.corr=readNIfTI("001_structural_corrected.nii.gz", reorient=F)
nii=readNIfTI("P [001] Series 006 [MR - 3D MP RAGE].nii", reorient=F)
sub.bias<-niftiarr(nii.corr,nii.corr-nii)
q=quantile(sub.bias[sub.bias!=0], probs=seq(0,1,by=0.1))
fcol=div_gradient_pal(low="blue", mid="yellow", high='red')
ortho2(nii.corr, sub.bias, col.y = alpha(fcol(seq(0,1,length=10)),0.5),
       ybreaks=q, ycolorbar=T, text=paste0("Original Image Minus N4", "n\ Bias-Corrected Image"))
#---2nd image---
nii.corr=readNIfTI("002_structural_corrected.nii.gz", reorient=F)
nii=readNIfTI("P [002] Series 006 [MR - 3D MP RAGE].nii", reorient=F)
sub.bias<-niftiarr(nii.corr,nii.corr-nii)
q=quantile(sub.bias[sub.bias!=0], probs=seq(0,1,by=0.1))
fcol=div_gradient_pal(low="blue", mid="yellow", high='red')
ortho2(nii.corr, sub.bias, col.y = alpha(fcol(seq(0,1,length=10)),0.5),
       ybreaks=q, ycolorbar=T, text=paste0("Original Image Minus N4", "n\ Bias-Corrected Image"))
#---3rd image---
nii.corr=readNIfTI("003_structural_corrected.nii.gz", reorient=F)
nii=readNIfTI("P [003] Series 006 [MR - 3D MP RAGE].nii", reorient=F)
sub.bias<-niftiarr(nii.corr,nii.corr-nii)
q=quantile(sub.bias[sub.bias!=0], probs=seq(0,1,by=0.1))
fcol=div_gradient_pal(low="blue", mid="yellow", high='red')
ortho2(nii.corr, sub.bias, col.y = alpha(fcol(seq(0,1,length=10)),0.5),
       ybreaks=q, ycolorbar=T, text=paste0("Original Image Minus N4", "n\ Bias-Corrected Image"))
#--------------------------------------------------------------------
rm(nii, nii.corr)
gc()
setwd(dir=nhack.dir)

nii2=readNIfTI("002_structural_corrected.nii.gz", reorient=F)
nii3=readNIfTI("003_structural_corrected.nii.gz", reorient=F)








#---white matter of the brain---
is_btw_300_400<-((nii_T1>300)&(nii_T1<400))
nii_T1_mask<-nii_T1
nii_T1_mask[!is_btw_300_400]=NA
overlay(nii_T1, nii_T1_mask, z=11, plot.type="single")
overlay(nii_T1, nii_T1_mask)
orthographic(nii_T1, nii_T1_mask, xyz=c(200,220,11), text="Image overlaid with mask", text.cex=1.5)

#---data manipulations---
path='Neurohacking_data/kirby21/visit_1/113/'
fname='113-01-MPRAGE.nii.gz'
fullname=paste0(path, fname)
T1<-readNIfTI(fullname, reorient=F)
orthographic(T1)
fname2='113-01-MPRAGE_mask.nii.gz'
maskname=paste0(path, fname2)
mask=readNIfTI(maskname, reorient=F)
orthographic(mask)
#---highlight only masked areas---
masked.T1=T1*mask
orthographic(masked.T1)

#--just for illustration, no deep sence---
fname3="113-01-MPRAGE_processed.nii.gz"
fullname3=paste0(path, fname3)
T1.follow=readNIfTI(fullname3, reorient=F)
orthographic(T1.follow)
subtract.T1<-T1.follow-T1
min(subtract.T1)
max(subtract.T1)
orthographic(subtract.T1)

#---transformation and smoothing---
#---allows to transform intensity values---
#---affectr how light or dark regions of the brain appear,
#---depending on their intensity value and the transfer function used
im_hist<-hist(T1, plot=F)
par(mar=c(5,4,4,4)+0.3)
coll=rgb(0,0,1,0.5)
plot(im_hist$mids, im_hist$count, log="y", type="h", lwd=10, lend=2,
     col=coll, xlab="Intensity Values", ylab="Count (Log Scale)")
#---Log-scale Histogram with linear Transfer Function---
par(new=T)
curve(x*1, axes=F, xlab="", ylab=", col=2, lwd=3)")
axis(side=4, at=pretty(range(im_hist$mids))/max(T1), 
     labels=pretty(range(im_hist$mids)))
mtext("Original Intensity", side=4, line=2)
#---this defines a linear spline. Other definitions are possible
lin.sp<-function(x, knots, slope){
  knots<-c(min(x), knots, max(x))
  slopeS<-slope[1]
  for(j in 2:length(slope)) {slopeS<-c(slopeS, slope[j]-sum(slopeS))}
  rvals<-numeric(length(x))
  for (i in 2:length(knots)) {
    rvals<-ifelse(x>=knots[i-1], 
    slopeS[i-1]*(x-knots[i-1])+rvals, rvals)}
  return(rvals)
}
#TRANSFORMATIONS - Define a spline with two notes and three slopes
knot.vals<-c(.3,.6)
slp.vals<-c(1,.5,.25)
par(new=T)
curve(lin.sp(x,knot.vals, slp.vals), axes=F, xlab="",ylab="", col=2,lwd=3)
axis(side=4, at=pretty(range(im_hist$mids))/max(T1), 
     labels=pretty(range(im_hist$mids)))
mtext("Transformed Intensitiy", side=4, line=2)
#transformations
trans_T1<-lin.sp(T1, knot.vals*max(T1), slp.vals)
image(T1, z=150, plot.type='single', main="Original Image")
image(trans_T1, z=150, plot.type='single, main="Transformed Image')
#---defining knot and slopes, we obtain a great control over transformation---

orthographic(T1)

#---SMOOTHING---
library(AnalyzeFMRI)
smooth.T1<-GaussSmoothArray(T1, voxdim = c(1,1,1),
            ksize=11, sigma=diag(3,3), mask=NULL,
            var.norm=F)
orthographic(smooth.T1)
#--reasonable smoothing can reduce an amount of noise

#---MRI contrasts
path1='Neurohacking_data/BRAINIX/NIfTI/'
#fluid-attenuated inversion recovery (FLAIR)
fname4='FLAIR.nii.gz'
fullname=paste0(path1, fname4)
volume.f<-readNIfTI(fullname, reorient=F)
orthographic(volume.f)
volume.f<-cal_img(volume.f)
image(volume.f, z=12, plot.type='single')
#---T1 image
fname5='T1.nii.gz'
fullname=paste0(path1, fname5)
volume.t1<-readNIfTI(fullname, reorient=F)
volume.t1<-cal_img(volume.t1)
image(volume.t1, z=12, plot.type='single')
#---T2 image
fname6='T2.nii.gz'
fullname=paste0(path1, fname6)
volume.t2<-readNIfTI(fullname, reorient=F)
volume.t2<-cal_img(volume.t2)
image(volume.t2, z=12, plot.type='single')

#---PREPROCESSING---
library(fslr)
fname7='113-01-MPRAGE.nii.gz'
#path2='Neurohacking_data/kirby21/visit_2/113/'
fullname=paste0(path, fname7)
mim<-readNIfTI(fullname, reorient=F)
mean(mim)
fslstats(mim, opts= "-m")

devtools::install_github("stnava/ITKR")
#---transpose data using t(): faces "up", versus "right"---
dd1=dim(t(slice$img[[1]]))
image(1:dd1[1], 1:dd1[2], t(slice$img[[1]]), col=gray(0:64/64))
hist(slice$img[[1]][,],breaks=50, xlab='FLAIR', prob=T, col=rgb(0,0,1,1/4),main="")

#---header information---
hdr=slice$hdr[[1]]
class(hdr)
names(hdr)
fg=hdr$name
length(fg)
hdr[hdr$name=="PixelSpacing", "value"]
hdr[hdr$name=="FlipAngle",]

#---reading all slides of an image---
all_slices_T1=readDICOM("Neurohacking_data/BRAINIX/DICOM/T1/")
dim(all_slices_T1$img[[11]])
hdr1=all_slices_T1$hdr[[11]]
hdr1[hdr$name=="PixelSpacing", "value"]

#---converting DICOM to Nifti---
nii_T1=dicom2nifti(all_slices_T1)
dd12=dim(nii_T1)
class(nii_T1)
image(1:dd12[1], 1:dd12[2],nii_T1[,,11], col=gray(0:64/64), xlab='',ylab='')
writeNIfTI(nim=nii_T1, filename=fname)