#TreeGen

Interactive tree growth model with environment sensitivity

Owner: Damith Jinasena

Email: damithsj@gmail.com


##Description

I designed this model based on the DBM model which simulates the propagation of a lightning discharge. In this model I used Honda model to generate the binary tree geometry and random growth of branches was modeled using Monte-Carlo method where it was possible to simulate the biasing of branches with the resource availability.
First version of this algorithm was presented in [SIGGRAPH Asia 2012 posters](http://dl.acm.org/citation.cfm?id=2407196)

For my Mphil research I combined binary tree growth algorithm with DBM and created a new growth model which simulates the growth of a tree with environment interactions which allows a tree to adapt to the distribution of light and proximity to solid obstacles it its growth cycle.
This program generates POV-Ray file for each step of the tree growth so by rendering the generated POV file sequence, tree growth can be modeled step by step

##How to generate trees

*  `conf\initialize-params.cfg`  -this file has all the parameters you can control with the program. You can define the grid size (environment), Light Source configurations, Obstacle configurations, basic branch properties etc…
* Run the `Release\laplace_3d_sth_3_rand.exe`
* file will generate in the `output` folder with the addition of each branch to the tree.

I know this program is too buggy and need so much developments. But as a prototype it works great for my research scope and gives considerably good results. I really appreciate your comments and feedback for further improvements of this project.

Have a nice time!

DSJ
