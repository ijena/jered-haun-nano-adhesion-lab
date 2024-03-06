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
