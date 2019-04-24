#---this is a copy-to-Ubuntu file for Rstudio under Ubuntu run
#---this way we process in the devish skewed Ubunty Rstudio
#---fslr functions, which is impossible under W10

#----------------------------------------------------------------
library(fslr)
library(oro.nifti)
library(matrixStats)
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
fname="as001-0005-00001-000001-01.nii"
fslnii=readNIfTI(fname, reorient=F)
orthographic(fslnii)

InFiles012=c("as001-0005-00001-000001-01.nii", "as001-0005-00001-001001-01.nii", "as001-0005-00001-002001-01.nii")
InFiles02=c("as001-0005-00001-000001-01.nii",  "as001-0005-00001-002001-01.nii")
InFiles12=c("as001-0005-00001-001001-01.nii", "as001-0005-00001-002001-01.nii")
InFiles123=c("as001-0005-00001-001001-01.nii", "as001-0005-00001-002001-01.nii",'as001-0005-00001-003001-01.nii')
InFiles.all="*.nii"

fslmerge(infiles=InFiles.all, direction="t", outfile="probfunc", retimg=F)
fname="probfunc.nii.gz"
fslnii=readNIfTI(fname, reorient=F)
orthographic(fslnii)










fname1="001_structural.nii"
fname2="002_structural.nii"
fname3="003_structural.nii"
fslnii1=readNIfTI(fname1, reorient=F)
fslnii2=readNIfTI(fname2, reorient=F)
fslnii3=readNIfTI(fname3, reorient=F)
fast_img1=readNIfTI("001_structural_corrected.nii.gz", reorient=F)
fast_img2=readNIfTI("002_structural_corrected.nii.gz", reorient=F)
fast_img3=readNIfTI("003_structural_corrected.nii.gz", reorient=F)

#---image intensity correction-----------------------------------
fast_img=fsl_biascorrect(fslnii1, retimg=T)
wrfilename='001_structural_corrected'
writeNIfTI(fast_img, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/001_structural_corrected.nii.gz /mnt/e/Neuro_data/test_data_Skolkovo/structural")

fast_img=fsl_biascorrect(fslnii2, retimg=T)
wrfilename='002_structural_corrected'
writeNIfTI(fast_img, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/002_structural_corrected.nii.gz /mnt/e/Neuro_data/test_data_Skolkovo/structural")

fast_img=fsl_biascorrect(fslnii3, retimg=T)
wrfilename='003_structural_corrected'
writeNIfTI(fast_img, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/003_structural_corrected.nii.gz /mnt/e/Neuro_data/test_data_Skolkovo/structural")
#---final steps--------------------------------------------------
rm(fslnii1, fslnii2, fslnii3, fast_img)
gc()

#---skull stripping----------------------------------------------
#---with the center of gravity-----------------------------------
bet_fast=fslbet(infile=fast_img1, retimg=TRUE)
cog1=cog(bet_fast, ceil=T)
cog1=paste("-c", paste(cog1, collapse=" "))
bet_fast.cog=fslbet(infile=fast_img1, retimg=TRUE,opts=cog1)
wrfilename='001_structural_corrected_skull.stripped'
writeNIfTI(bet_fast.cog, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/001_structural_corrected_skull.stripped.nii.gz /mnt/e/Neuro_data/test_data_Skolkovo/structural")

bet_fast=fslbet(infile=fast_img2, retimg=TRUE)
cog1=cog(bet_fast, ceil=T)
cog1=paste("-c", paste(cog1, collapse=" "))
bet_fast.cog=fslbet(infile=fast_img2, retimg=TRUE,opts=cog1)
wrfilename='002_structural_corrected_skull.stripped'
writeNIfTI(bet_fast.cog, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/002_structural_corrected_skull.stripped.nii.gz /mnt/e/Neuro_data/test_data_Skolkovo/structural")

bet_fast=fslbet(infile=fast_img3, retimg=TRUE)
cog1=cog(bet_fast, ceil=T)
cog1=paste("-c", paste(cog1, collapse=" "))
bet_fast.cog=fslbet(infile=fast_img3, retimg=TRUE,opts=cog1)
wrfilename='003_structural_corrected_skull.stripped'
writeNIfTI(bet_fast.cog, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/003_structural_corrected_skull.stripped.nii.gz /mnt/e/Neuro_data/test_data_Skolkovo/structural")
