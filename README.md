# variable-stiffness
An Open-Source Framework for the FE Modelling and Optimal Design of Fiber-Steered Variable-Stiffness Composite Cylinders Using Water Strider Algorithm

* For running the algorithm, one should first create a file in the following directory:
D:\optLaminatedComp
Then, copy and paste the containing files in the algorithm folder into the above directory.

* The current algorithm ought to runs the verification case (L/R=2) of the article, having PSO algorithm. 

* Replace "pythonScript_Spline" with "pythonScript" in the following line of the file bucklingMoment,
in order to incorporate the ability to replace the spline fit of theta 
with linear: 
system(['abaqus cae nogui=pythonScript.py'])
