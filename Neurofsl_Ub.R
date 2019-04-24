#---this is a copy-to-Ubuntu file for Rstudio under Ubuntu run
#---this way we process in the devish skewed Ubunty Rstudio
#---fslr functions, which is impossible under W10
#---BE CAUTIOUS:any function which CAN be processed in W10 is EXCLUDED from this devish Ubuntu file
#---and processed in W10 in the file "NeuroWin10.R"
#---------------------------------------------------------------
#---so the steps follow:
#---edit this file in normal Rstudio for W10 and save;
#---after that copy this file to Ubuntu home directory via Ubuntu terminal commands:
#---cd /mnt/c/Users/Belobragin\ Vladimir/Documents/R/working/Neurohacking/
#---cp NeuroUb.R /home/bvv
#---then load this file in Rstudio under Ubuntu (do your best) and run it 
#---after acheived result close Rstudio under Ubunty without save and without regret
#---and copy back to W10 processed *.nii file (or some other) with
#---cp /home/tempnii(*).nii.gz /mnt/c/Users/Belobragin\ Vladimir/Documents/R/working/Neurohacking/
#---etc.
#---and forget about Ubuntu till next time
#---so, once again: this file contains only commands which are absolutely nessesary to run under Ubuntu
#---and read - write commands
#----------------------------------------------------------------
library(fslr)
library(oro.nifti)
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
