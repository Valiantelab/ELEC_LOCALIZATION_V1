#!/bin/sh

# ct2mri.sh
# 
#The second argument is the nii.gz file of the CT scan that you need to create with something like Matlab or FSL's mri_convert function.
#This script uses FSL's flirt command with max flexibility (i.e., 12 degrees of freedom mapping) transform the CT scan so that it lines up with the preimplant MRI by maximizing the mutual information between the volumes.
#Images of the two volumes overlayed are automatically generated so that you can inspect the quality of the coregistration. 
#In the process the elec_recon subfolder in the patient's FreeSurfer folder is created along with the following files:
#    T1.nii.gz: The full head MRI
#    brainmask.nii.gz: The skull stripped MRI
#    ctINt1.nii.gz: The post-implant CT coregistered to the pre-implant MRI
#
# Created by David Groppe on 2/11/15.
# Questions? Email: david.m.groppe@gmail.com
# Copyright 2015 __MyCompanyName__. All rights reserved.

usage='\nUSAGE:\n  ct2mri12dfFast.sh freesurferSubjectName ctNiiFile\n\nEXAMPLE:\n ct2mri12df.sh TWH014 /Users/dgroppe/Desktop/TWH_14_DICOMS/ct.nii.gz\n'

if [[ "$#" = 0 ]]; then
 echo $usage
 exit 2
fi

sub=$1
fsDir=$SUBJECTS_DIR$sub
if [ ! -d $fsDir ]; then
  echo
  echo "... ${fsDir} is not a directory."
  echo "...you have the wrong FreeSurfer codename for the subject or you have not yet run recon-all on this subject's MRI " 
  echo
  exit 2
fi

if [ ! -f  $2 ]; then
 echo
 echo "...File ${2} not found. Exit."
 echo
 exit 2
fi

elecReconPath=$SUBJECTS_DIR$sub/elec_recon
echo 'Creating directory ' $elecReconPath
mkdir $elecReconPath

echo 'Creating T1.nii.gz in elec_recon folder for coregistration.'
mriPath=$SUBJECTS_DIR$sub/mri
mri_convert $mriPath/T1.mgz $elecReconPath/T1.nii.gz

echo 'Creating brainmask.nii.gz in elec_recon folder for use with BioImageSuite later.'
mri_convert $mriPath/brainmask.mgz $elecReconPath/brainmask.nii.gz

echo 'Copying CT nii.gz file to elec_recon folder.'
cp $2 $elecReconPath/.

echo 'Registering ' $2 ' to T1.nii.gz with a maximally flexible linear (12 degrees of freedom) transformation that maximizes mutual information between the volumes. This takes awhile....'
flirt -in $2  -ref $elecReconPath/T1.nii.gz -out $elecReconPath/ctINt1.nii.gz -omat $elecReconPath/ct2t1.mat -interp trilinear -cost mutualinfo -dof 12
slices $elecReconPath/ctINt1.nii.gz $elecReconPath/T1.nii.gz
slices $elecReconPath/T1.nii.gz  $elecReconPath/ctINt1.nii.gz  
echo 'Run the command below to interactively inspect the coregistration:'
echo "fslview ${elecReconPath}/T1.nii.gz ${elecReconPath}/ctINt1.nii.gz" 
