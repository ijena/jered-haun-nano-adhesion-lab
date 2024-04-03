#2/28/24 
Questions and errors

How exactly is the gaussian distribution supposed to work?
We talked about np.linsppace but should it be np.random.normal? If yes then is midpoint 500 and what is standard deviation

#2/28/24

Calculate constant sigma and use that to create a range of  - 3 sigma to 3 sigma.
Pick one sigma in that range and use it to pass as the sigma parameter for the gaussian function and call that gaussian function in position/velocity calculations
Fix how are you keeping your particle in bounds

Also use previous velocity and position for the mean in the gaussian

#3/6/24

What's wrong with velocity graphs?
Sigma was supposed to be between -3sigma and 3 sigma but I took absolute value because sigma cannot be negative
UROP thing
BME 199 course form
# TO DO - add normal distribution to the histogram

#3/6/24

Fix velocity by changing bounds and checking logic
Try playing around with the bins of historgram
Put an overlay of normal distribution on the histogram
Find papers of people who have done nanorods, nanospheres simulations
Make ppt of the plots
Work on the UROP proposal
Plot the average position at every time as a line graph on the position vs time graph
#Goals
#Working nanosphere simulator which is validated and have plots with it in 3 dimensions
Hooke's and Spring development of attachment dynamics on the bond surface
Same things for rods instead of spheres

#3/19/24
Questions 
Figure how to get the overlay in the distribution (maybe scale factor?)
Planning to change bins how i get the overlay done 
Was unable to fix velocity issues but can look into it again

#4/3/24

Work on Project Abstract (overview, your role and the impact and results with the detachment curve)
Try creating your own function for normal distribution histogram overlay instead of using norm.pdf
Create an absolute position of the particle vs time graph, Slope should be ndt = 1dt d is the diffusion coefficient
Continue working on the velocity graphs 
Stokes-einstein equation to get diffusion coefficient for a nanoparticle in water


