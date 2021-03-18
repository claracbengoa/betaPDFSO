# β-PDFSO

Reliability Index based strategy for the **P**robability **D**amage approach in **F**ail-**S**afe **O**ptimization

The β-PDFSO is a new fail-safe optimization strategy to include available information on the probability of occurrence of different accidental scenarios as well as uncertainty in parameters affecting structural responses. The optimization approach avoids obtaining oversized designs, as the value of the objective function is reduced compared to the fail-safe RBDO. Due to the lack of knowledge of which damaged configuration will occur, a new random reliability index 
<img src="https://render.githubusercontent.com/render/math?math=\hat{\beta}">
is defined for each limit-state of the damaged structure. This new random variable is constructed through the 
<img src="https://render.githubusercontent.com/render/math?math=\beta">
in each damaged configuration. The method guarantees 
<img src="https://render.githubusercontent.com/render/math?math=\beta^T">
in the limit-states of the intact configuration and 
<img src="https://render.githubusercontent.com/render/math?math=P[\beta < \beta^T] \le p_f^T">
in limit-states affecting the damaged structure.

## Example

An optimization problem is defined using mathematical expressions to simulate the objective function 
<img src="https://render.githubusercontent.com/render/math?math=F">
and the structural responses 
<img src="https://render.githubusercontent.com/render/math?math=R">. 
These structural responses are defined by a set of equations, which depend on the design variables 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{d}=[d_1, d_2]">, 
the random variables 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{x}=[x_1, x_2]"> 
and a vector of coefficients 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{c}=[c_1, c_2, c_3, c_4, c_5, c_6]">, 
as shown in Eq.1, where a response 
<img src="https://render.githubusercontent.com/render/math?math=R_{h,k}^{M_i}"> 
is defined per load case 
<img src="https://render.githubusercontent.com/render/math?math=h"> 
and limit-state 
<img src="https://render.githubusercontent.com/render/math?math=k">. 
By varying the vector 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{c}"> 
the structural response 
<img src="https://render.githubusercontent.com/render/math?math=R_{h,k}^{M_i}"> 
is modified, being possible to obtain a different value for the intact and damaged configurations: the vector 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{c}"> 
is taken as the unit vector to represent the structural response in the intact model (
<img src="https://render.githubusercontent.com/render/math?math=M_0">
), while values between 0 and 1 are used to simulate the response in the damaged configurations (
<img src="https://render.githubusercontent.com/render/math?math=M_i, i=1,...,D">
). By adopting this approach, the structural responses defined in Eq. 1 increase for values of 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{c}"> 
lower than 1, being possible to simulate a loss of structural capacity in the damaged configurations. In this problem, 
<img src="https://render.githubusercontent.com/render/math?math=D=60"> 
was taken as the number of damaged configurations.

![](refs/eq1.png)

The formulation of the optimization problem is presented in Eq. 2. The objective is to minimize the objective function 
<img src="https://render.githubusercontent.com/render/math?math=F">, defined as the sum of the design variables 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{d}=(d_1, d_2)">.
The random variables are defined as normally distributed with mean 
<img src="https://render.githubusercontent.com/render/math?math=\bm{\mu}_{\textbf{x}}=(0.25, 1)">
and standard deviation 
<img src="https://render.githubusercontent.com/render/math?math=\bm{\sigma}_{\textbf{x}}=(0.025, 0.1)">. 
The probability of occurrence 
<img src="https://render.githubusercontent.com/render/math?math=P_{M_i}"> 
of each damaged configuration is shown in Table 1.

![](refs/eq2.png)
![](refs/table1.png)
![](refs/figure1.png)

The vector of coefficients defined in Eq. 1, 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{c}=[c_1, c_2, c_3, c_4, c_5, c_6]">, 
is set for each load case <img src="https://render.githubusercontent.com/render/math?math=h">, 
limit-state <img src="https://render.githubusercontent.com/render/math?math=k"> and each model 
<img src="https://render.githubusercontent.com/render/math?math=M_i">, 
expressed as <img src="https://render.githubusercontent.com/render/math?math=\textbf{c}_{h,k}^{M_{i}}"> . Only some values of these coefficients are summarized in the arrays presented in Eq. 3. The full length matrices are available in the matlab file main.m.

![](refs/eq3.png)

The optimization problem was solved for a target reliability index <img src="https://render.githubusercontent.com/render/math?math=\beta^T=3.7190"> and a target probability of failure of <img src="https://render.githubusercontent.com/render/math?math=p_{f}^{T}=0.10">, which correspond to a confidence level, 
<img src="https://render.githubusercontent.com/render/math?math=R^{T}=0.90">.




## How to use it

Edit the `defineInputParameters` function, where the following inputs are defined:

Problem Parameters       | Description
:-----------------       | :----------
d0                       | Initial value of design varaibles
nDV                      | Number of design variables
nLC                      | Number of load cases
nDcon                    | Number of limit-states
nDamages                 | Number of damaged configurations
pDamages                 | Probability of occurrence of each damaged configuration
pfPDFSO                  | Target probability of failure in the β-PDFSO 
mu_xi                    | mean of the random variable i
sig_xi                   | standard deviation of the random variable i
beta_min                 | target reliability index for RA
resp                     | Structural Response (cell) 
respMax                  | Maximum allowable response (cell)  
GFun                     | Limit state (cell)
Cfd0                     | coefficients for the intact model (cell)
CfDamages                | coefficients for the damaged configurations (cell)
dGFun_dresp              | derivative of GFun with respect to resp
dresp_dxi                | derivative of resp with respect to xi    

Optimization Parameters       | Description
:----------------------       | :----------
lb                            | Lower bounds of design variables
ub                            | Upper bounds of design variables
TolCon                        | Constraint tolerance
TolFun                        | Objective function tolerance
DiffMinChange                 | Minimum finite difference step
DiffMaxChange                 | Maximum finite difference step


## Results

The program generates a CFD plot of each probabilistic constraint and the file opt_results_it.txt.

![](refs/table2.png)

![](refs/figure2.png)
