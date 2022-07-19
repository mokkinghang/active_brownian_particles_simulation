# Active Brownian particles simulation

This is a simulation code for simulating a group of interacting active Brownian particles governed by the Langevin equations

<img src="https://latex.codecogs.com/svg.image?\dot&space;\vec&space;r&space;=&space;v_0&space;\hat{e}&plus;\sqrt{2D_T}\vec\xi&space;" title="\dot \vec r = v_0 \hat{e}+\sqrt{2D_T}\vec\xi " />
<img src="https://latex.codecogs.com/svg.image?\dot&space;\varphi&space;=&space;\sqrt{2D_R}\chi," title="\dot \varphi = \sqrt{2D_R}\chi," />

where

<img src="https://latex.codecogs.com/svg.image?\hat{e}&space;=&space;(\cos&space;\varphi,&space;\sin&space;\varphi)" title="\hat{e} = (\cos \varphi, \sin \varphi)" />. 

They are subjected to the WCA potential, which is a variant of the Lennard-Jones potential:

<img src="https://latex.codecogs.com/svg.image?&space;&space;&space;&space;V_{WCA}(r)=\left\{&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;\begin{array}{ll}&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;4\varepsilon((\frac{\sigma}{r})^{12}-(\frac{\sigma}{r})^6)&plus;\varepsilon,&space;&&space;r&space;\leq&space;2^{1/6}\sigma&space;\\&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;0,&space;&&space;r&space;>&space;2^{1/6}\sigma&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;\end{array}&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;\right." title=" V_{WCA}(r)=\left\{ \begin{array}{ll} 4\varepsilon((\frac{\sigma}{r})^{12}-(\frac{\sigma}{r})^6)+\varepsilon, & r \leq 2^{1/6}\sigma \\ 0, & r > 2^{1/6}\sigma \end{array} \right." />

This interaction potential should retrieve the collision dynamics of hard spheres. 

In the simulation, we need a nondimensionalized system of equations of motion, which can be expressed in terms of the Peclet number \textit{Pe}, defined as

<img src="https://latex.codecogs.com/svg.image?Pe&space;:=&space;\sigma&space;v_0/D_T." title="Pe := \sigma v_0/D_T." />

The nondimensionalized Langevin equations are

<img src="https://latex.codecogs.com/svg.image?\frac{d&space;\vec&space;r'}{d&space; t'}&space;=&space;Pe&space;\;&space;\hat{e}&space;&plus;&space;\sqrt{2}{\vec&space;\xi'}&space;" title="\frac{d \vec r'}{d \bar t} = Pe \; \hat{e} + \sqrt{2}{\vec \xi'} " />
<img src="https://latex.codecogs.com/svg.image?\frac{d&space;\varphi'}{dt'}&space;=&space;\sqrt{\frac{3}{2}}\frac{\sigma}{R}&space;\chi'&space;" title="\frac{d \varphi'}{dt'} = \sqrt{\frac{3}{2}}\frac{\sigma}{R} \chi' " />


## Task 1: Langevin Equation

In this 
To plot the trajectories of five particles for Pe=0 and Pe=20, run in command line `julia plot_trajectories.jl`.

To plot the MSD and MSAD plots, run in command line `julia plot_msd_msad.jl`.
