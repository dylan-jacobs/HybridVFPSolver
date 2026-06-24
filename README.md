# Hybrid Vlasov-Fokker-Planck Solver
Solves the Vlasov Fokker-Planck equation in cylindrical coordinates using a hybrid model; the ion distributions are described using a kinetic model, while the electrons are described using a fluid model. We test our solver on three problems: the heat equation, the 0D2V Lenard-Bernstein-Fokker-Planck equation, and the 1D2V Vlasov-Fokker-Planck equation.

## Heat Equation
Everything for the heat equation test is located in the ./Heat-Equation/main_heat_equation.m MATLAB script. It can be set up to test either the solver's temporal accuracy or its mass conservation and rank-adaptivity. There are different parameters for each case, commented accordingly. To test the temporal accuracy, initialize the variable ```lambdavals``` to a vector (line 14) and ensure you are plotting the temporal accuracy graph in the "Figures" section of the script. 

### Example figures
<p float="left">
  <img width="30%" src="https://github.com/user-attachments/assets/53ae59a1-6904-43c0-89be-5d527f69ac84" />
  <img width="30%" src="https://github.com/user-attachments/assets/fb2e7e37-7aba-4b07-845e-25180a6954f9" />
</p>

## 0D2V Lenard-Bernstein-Fokker-Planck (or Dougherty-Fokker-Planck) Equation
The Dougherty-Fokker-Planck (DFP) equation tests are located in the ./Dougherty-Fokker-Planck/main_DFP.m MATLAB script. The script can be set up to test either the solver's temporal accuracy or its other tests, such as mass, momentum, and energy conservation, L1 drive to equilibrium, relative entropy decay, and rank-adaptivity. There are different parameters for each case, commented accordingly in the script. To test the temporal accuracy, initialize the variable ```lambdavals``` to a vector (line 3) and ensure you are plotting the temporal accuracy graph in the "Figures" section of the script. Because the temporal accuracy test also requires the reference solution (./Dougherty-Fokker-Planck/DFP-Reference-Solution/dfp_refsoln_dirk3.mat), it is necessary to ensure that the mesh parameters in the tested simulation match those of the reference solution. You can run and save your own reference solutions using line 138. The ```main_DFP.m``` script uses the quadrature-corrected-moments (QCM) procedure to slightly perturb the solution's numerical moments so that the moments generating the Maxwellian distribution exactly match that distribution's _numerical_ moments. Thus, it is important that the ```QCM.m``` script is located in the same directory as ```main_DFP.m```.

### Example figures
<p float="left">
  <img width="30%" src="https://github.com/user-attachments/assets/56fe9572-8e4f-40f8-af25-92eb7d0d8d82" />
  <img width="30%" src="https://github.com/user-attachments/assets/093b686e-a17e-45a9-a787-4a993e1f3a81" />
  <img width="30%" src="https://github.com/user-attachments/assets/30ee4fc6-0158-4b55-84f0-f675bd79a88e" />
</p>

## 1D2V Vlasov-Fokker-Planck Equation
### Standing Shock Test
The VFP Standing Shock tests are located in the ./Standing-Shock/main_VFP.m MATLAB script. The script can be set up to test either the solver's temporal accuracy or its other tests, such as mass, momentum, and energy conservation, and rank-adaptivity. There are different parameters for each case, commented accordingly in the script. To test the temporal accuracy, initialize the variable ```DTvals``` to a vector (line 7) and ensure you are plotting the temporal accuracy graph in the "Figures and Tables" section of the script. Because the temporal accuracy test also requires the reference solution (./Standing-Shock-Problem/VFP-Reference-Solution/vdfp_refsoln_IMEX222_QCM_Nx160_Nvperp200_Nvpar200_T1_tole-8_dt10-3.mat), it is necessary to ensure that the mesh parameters in the tested simulation match those of the reference solution. It is also necessary to uncomment lines 413-434 to calculate the temporal errors during the simulation. When running all other tests besides temporal accuracy, these lines (413-434) should be commented out because the reference solution mesh parameters may not match the simulated parameters and will cause errors. You can run and save your own reference solutions using line 138. The ```main_VFP.m``` script uses the quadrature-corrected-moments (QCM) procedure to slightly perturb the solution's numerical moments so that the moments generating the Maxwellian distribution exactly match that distribution's _numerical_ moments. Thus, it is important that the ```QCM.m``` script is located in the same directory as ```main_VFP.m```.
