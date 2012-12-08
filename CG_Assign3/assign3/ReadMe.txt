Assignment #3: Ray tracing



FULL NAME: Shengzhe Chen  USCID: 4341396265

MANDATORY FEATURES

------------------



<Under "Status" please indicate whether it has been implemented and is

functioning correctly.  If not, please explain the current status.>



Feature:                                 Status: finish? (yes/no)

-------------------------------------    -------------------------

1) Ray tracing triangles                  yes
2) Ray tracing sphere                     yes

3) Triangle Phong Shading                 yes


4) Sphere Phong Shading                   yes


5) Shadows rays                           yes

6) Still images                           yes
   

7) Extra Credit (up to 20 points)
   

Besides with the above madantory features, I have added soft shadow, recursive ray tracing and, anti-alias.

For soft shadow part, the volume of a point light is set to 0.05, you can simply modify the LIGHT_VOLUME to have better result if the light position has been changed.And a 
intersected point will emit 20 random vectors to check the visibility of light, modify the RANDOM_LIGHTS to match your requirement.

For recursive ray tracing part, depth of ray trace is hard coded to 5 times.

For anti-alias part,the number of samples of each pixel is set to 9, and samples are randomly picked up within each pixel cell, you can modify the SAMPLING_TIMES to have better result if you increase the SAMPLING_TIMES.