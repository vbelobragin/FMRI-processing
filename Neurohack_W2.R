library (tidyverse)
library(plyr)
detach(package:plyr)
library (dplyr)

library("oro.dicom")
slice=readDICOM("IM-0001-0011.dcm")

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

#---working with nifti files---
library(oro.nifti)
fname="Neurohacking_data/BRAINIX/NIfTI/Output_3D_File"
writeNIfTI(nim=nii_T1, filename=fname)
nii_T2=readNIfTI("Neurohacking_data/BRAINIX/NIfTI/T2.nii.gz", reorient=F)
#---to not compress, obtain gzipped=F in the function writeNIfTI
#---r2agui.sourceforge.net

#---basic visualization---
print({nii_T1})
image(1:dd12[1], 1:dd12[2],nii_T1[,,11],  xlab='',ylab='')
graphics::image
heat.colors(12)
image(nii_T1, z=11, plot.type='single')
oro.nifti::image
image(nii_T1)
orthographic(nii_T1, xyz=c(200, 220, 11))
#---xyz defines the cross of red lines---

#---vizualisation---
par(mfrow=c(1,2))
o<-par(mar=c(4,4,0,0))
hist(nii_T1, breaks=75, prob=T, xlab="T1 intensities", col=rgb(0,0,1,0.5), main="")
hist(nii_T1[nii_T1>20], breaks=75, prob=T, xlab="T1 intensities>20", col=rgb(0,0,1,0.5), main="")

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

#---template image processing---
pathname="Neurohacking_data/Template/"
fname="MNI152_T1_1mm_brain.nii.gz"
niifile=paste0(pathname,fname)
templNii=readNIfTI(niifile, reorient=F)
print({templNii})
dd12=dim(templNii)
image(1:dd12[1], 1:dd12[2],templNii[,,100],  xlab='',ylab='')
graphics::image
heat.colors(12)
image(templNii, z=11, plot.type='single')
oro.nifti::image
image(templNii)
orthographic(nii_T1, xyz=c(200, 220, 11))