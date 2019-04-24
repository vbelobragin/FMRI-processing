#---this is a copy-to-Ubuntu file for Rstudio under Ubuntu run
#---for the package ANTsR
#---be cautious: NO ANTsR package is installed under W10---------------------------
#---so, it is not possible to split work with ANTsR among W10 and devish Ubuntu---
#---ALL functions of ANTsR must be implemented in this file
#----------------------------------------------------------------
library(oro.nifti)
library(extrantsr)
library(fslr)
T1_file="113-01-MPRAGE.nii.gz"
T1=readNIfTI(T1_file, reorient=FALSE)
T2_file="113-01-T2w.nii.gz"
T2w=readNIfTI(T2_file)
flirt_reg_t2_img=flirt(infile=T2_file, reffile=T1, dof=6, verbose = FALSE)
double_ortho(T1, flirt_reg_t2_img)



#---bias field correction---
#---lecture 3.9.------------
#---extrantsr::bias_correct wraps n3BiasFieldCorrection and n4BiasFieldCorrection from ANTsR for bias field correction^
n3img=
fname="113-01-MPRAGE.nii.gz"
fslnii=readNIfTI(fname, reorient=F)
print({fslnii}) #---simple check

#---image intensity correction-----------------------------------
fast_img=fsl_biascorrect(fslnii, retimg=T)

#---skull stripping----------------------------------------------
bet_fast=fslbet(infile=fast_img, retimg=TRUE)
#---now with the center of gravity-------------------------------
cog1=cog(bet_fast, ceil=T)
cog1=paste("-c", paste(cog1, collapse=" "))
bet_fast2=fslbet(infile=fast_img, retimg=TRUE,opts=cog1)

#---non-linear image registration
niifile="MNI152_T1_1mm_brain.nii.gz"
template.nii=readNIfTI(niifile, reorient=F)
orthographic(template.nii)
fnirt_fast=fnirt_with_affine(infile=bet_fast2, reffile=template, outfile="FNIRT_to_Template", retimg=TRUE)
orthographic(template.nii, fnirt_fast)

#---final steps--------------------------------------------------
wrfilename='tempnii1'#--'1' after l 3.4
writeNIfTI(fast_img, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/tempnii1.nii.gz /mnt/e/Neuro_data/Neurocoductor_data/NeuroH_MOOC_files")

wrfilename='tempnii2'#--'2' after l 3.5
writeNIfTI(bet_fast, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/tempnii2.nii.gz /mnt/e/Neuro_data/Neurocoductor_data/NeuroH_MOOC_files")

wrfilename='tempnii21'#--'21' after cog l 3.5
writeNIfTI(bet_fast, filename = wrfilename, verbose=TRUE)
system("cp /home/bvv/tempnii21.nii.gz /mnt/e/Neuro_data/Neurocoductor_data/NeuroH_MOOC_files")
