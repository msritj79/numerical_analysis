# numerical_analysis

Calculate the hydrodynamic pressure around dimples under full-film lubrication

1. Solving the Reynolds equation in cylindrical coordinates
2. Ignore surface roughness and asperity contact
3. The cavitation boundary condition follows mass-conservative JFO model


## Equation  
Reynolds equation coupled with a mass-conservative cavitation algorithm

<img src="image/formula/1.png">

　　For complete film zones:  
　　　<img src="image/formula/2.png">,　　<img src="image/formula/3.png">, 　<img src="image/formula/4.png">


　　For cavitated zones:  
　　　<img src="image/formula/5.png">,　　<img src="image/formula/6.png">,　　<img src="image/formula/7.png">


## DEMO


**Parameter**  

<img src="image/parameter.png" width=350px>


**Output**

|**pressure_3D_grapgh**|**pressure_2D_pressure**|**depth-pressure_graph**|
|---|---|---|
|<img src="image/pressure_3Dgraph_depth025.png">|<img src="image/pressure_2Dgraph_depth025.png">|<img src="image/depth-meanp_n750.png" width=500px>|
