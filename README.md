# radiosity
Matlab code for solving the radiosity equation for lighting in computer graphics virtual spaces. 

The code is related to the examples I show in the YouTube video https://youtu.be/krIVZvzlxUQ Note that the video language is Finnish but that it has English subtitles. 

Note that Matlab's own betacdf function is in a toolbox that not everyone necessarily have. So I included the file Scaled_betaCDF.m by Greig A. Paterson, making this code package independent of that toolbox. Also, the license of Paterson's code is included. The role of the beta cdf function is to optimise the colors of virtually lit spaces for display. 

1. Simplest example: empty room. 

First run "radiosity_emptyroom_Fcomp.m" in Matlab to create the geometric form factor matrix F, which will be saved to a file in the subfolder ./data/. Then you can run "radiosity_emptyroom_BW.m" to create a black-and-white image of a lit room, and "radiosity_emptyroom_color.m" for a color image. If you want more details (smaller patches), you can make n larger in "radiosity_emptyroom_Fcomp.m". Note that too large n value will cause your computer to run out of memory, so it is a good idea to increase n gradually. 

2. Another example: room with a dividing wall. 

First run "radiosity_wall_Fcomp.m" in Matlab to create the geometric form factor matrix F, which will be saved to a file in the subfolder ./data/. Then you can run "radiosity_wall_color.m" to create a color image. If you want more details (smaller patches), you can make halfn larger in "radiosity_wall_Fcomp.m". Note that too large halfn value will cause your computer to run out of memory, so it is a good idea to increase halfn gradually. 

3. Third example: room with a (levitating) table. 

First run "radiosity_table_Fcomp.m" in Matlab to create the geometric form factor matrix F, which will be saved to a file in the subfolder ./data/. Then you can run "radiosity_table_color.m" to create a color image. If you want more details (smaller patches), you can make halfn larger in "radiosity_wall_Fcomp.m". Note that too large halfn value will cause your computer to run out of memory, so it is a good idea to increase halfn gradually. 
