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

4/10/24
Which Stokes-Einstein equation to use
Implemented my own normal distribution but it isnt working correctly

Meeting notes
Take the position values given from the simulation and add up the absolute values of the position at each second and plot a linear graph of that and calculate the slope of that instant
Fix the slope to work on the analysis graph
UROP proposal 
Use stokes-einstein-sutherland equation from wikepedia
Plot the analysis graph for velocities too but keep changing the mean of the velocities from 10^1 to 10^10 <br>

4/16/24
Is the slope analysis for velocity also the same formula of Slope = nDt?
Some values cannot be plotted, is it because of infinite values?
What is the meeting location
Continue debugging

Meeting notes <br>
4/17/24 <br>
t in ndt is temperature <br>
Use scipy.optimize.curve_fit to calculate slope of the cumualtive position graph <br>
Change position vs time to a scatter <br>
Put graphs into slideshow <br>
Change total time and see if the histograms change <br>
Change gaussian distribution mean in position and velocity to be the previous position and velocity <br>

Questions <br>
4/23/24 <br>
Equation to find the ideal slope of the cumulative absolute velocity graphs <br>

Meeting Notes 4/24

Implement bivariate normal distribution, pick 2 values from normal distribution and put it into the bivariate normal distribution https://mathworld.wolfram.com/BivariateNormalDistribution.html
Look into Mersenne Twister algorithm to get random numbers initially
Cumulative position ideal slope (refer to picture) we are trying to solve for x and t is total time taken and T is temperature so fix ideal slope based on that
Continue debugging

Questions for 5/8 <br>
How exactly is the bivariate normal distribution working since I get 2 values and I am averaging them <br>
Overflow error for velocity <br>
How exactly is all the normal distribution stuff working? <br>
UROP Poster graphs and printing <br>

Meeting Notes for 5/9 <br>
Refer to Page 347 of Allen- Tildesley paper to pick the normal distribution number and use those equations to put in the velocity and position equations, let the mean be 0 initially <br>
Try mean of previous position and velocity too <br>

Questions for 5/15 <br>
Issues with normal distribution being previous position and velocity <br>
Confirm understanding of box and muller method <br>

Meeting notes for 5/23 <br>
Remove unnecessary comments and code <br>
Fix formatting <br>
Delta t in sigma position and sigma velocity should be current time interval <br>
Make a function that returns every value that is varying <br>
In the loop just call that function <br>
Remove random_sigma_position and random-sigma_velocity and just use sigma_position and sigma_velocity as the standard deviation as that gaussian distribution <br>
Normal distribution mean should be at 500 for position and 0 for velocity (check midpoint) <br>

Questions for 5/30 Meeting <br>
Particle mass should include a -12 I believe but it leads to math domain error <br>
Should c0, c1 and c2 also vary with time? <br>
c0,c1 and c2 are constant values because inverse viscous relaxation time and is way too small compared to time <br>
Sigma velocity is 0 because of the e term with inverse viscous relaxation time * time <br>
