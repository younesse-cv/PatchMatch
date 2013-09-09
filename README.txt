=====================================
README:   Inpainting based PatchMatch    
=====================================

@Author: Younesse ANDAM

@Contact: younesse.andam@gmail.com 
   
Description: This project is a personal implementation of an algorithm called PATCHMATCH that restores missing areas in an image. 
The algorithm is presented in the following paper 
 PatchMatch  A Randomized Correspondence Algorithm          
               for Structural Image Editing                       
   by C.Barnes,E.Shechtman,A.Finkelstein and Dan B.Goldman        
   ACM Transactions on Graphics (Proc. SIGGRAPH), vol.28, aug-2009  
   
 For more information please refer to                
 http://www.cs.princeton.edu/gfx/pubs/Barnes_2009_PAR/index.php   

Copyright (c) 2010-2011   

Requirements 
============

To run the project you need to install Opencv library and link it to your project. 
Opencv can be download it here 
http://opencv.org/downloads.html

How to use 
===========

The project accepts two images 
1- The original image 
2- The pruned image you can delete a part of interest in the image. The algorithm will patch the remaining image to give a natural result. 

The project contains some example of images to try it. You may find them in image_files. 
After choosing the image file, enter the paths of those image files in main.c 

char fileNameInput[50] = YOUR_PATH_HERE_OF_ORIGINAL_IMAGE;
char fileNameMasked[50] = YOUR_PATH_HERE_OF_PRUNED_IMAGE;


Enjoy!!
