This folder is where the parameters of the simulation are mostly set

Some parameters change with the case that we are running.
These are held in the runParams structure.

getParameters is the main function and calls the rest of the functions.

In the function getHydraulic, the functions makeSwitchLossMap and makeEHALossMap are called
It is in this function that you decide whether to make a new SwitchMap.mat or use the previous one