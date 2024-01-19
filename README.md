# ECMO_Sim

Here is my MATLAB script simulating V-V ECMO circuit and its various effects on overall oxygen balance. I implemented the code based on researchers' prior works 
- Zanella et al Journal of Critical Care 36 (2016) 178–186
- Spinelli and Bartlett, ASAIO Journal 60:6 (2014) 688-693

The implementation is relatively simple and based on oxygen mass balance and Fick principle. Based on specified oxygen consumption rate, cardiac output, hemoglobin, ECMO blood flow, and other input parameters, the code outputs the steady-state true SaO2 and SvO2 (arterial and venous saturation). The code also accounts for recirculation effect that is commonly seen in V-V support.

Some assumptions made for simplicity:
- dissolved oxygen component is neglected, which ignores about ~10% of total oxygen
- Outlet saturation of membrane oxygenator is set at 100% 
