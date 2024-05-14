**Title**
Optimal sizing and operation of hydrogen generation sites accounting for waste heat recovery

**Abstract**
This work proposes a mixed integer linear program designed to optimize the size and operation of the components of a hydrogen generation site for minimal cost. The model investigates the recovery of heat generated during electrolyser operation as a strategy to lower the levelized cost of hydrogen. The recovered heat can be exported to a medium-temperature (35 degC) and/or high-temperature (65 degC) district heating network. Furthermore, it incorporates accurate representations of the component cost and operation curves, with a particular focus on the electrolyser.

**Installation**
In order to run the optimization, MATLAB (including the optimization toolbox) and Gurobi should be installed. 

**Usage**
Begin by updating the path to Gurobi in the code to match its location on your laptop.
Navigate to the ‘INPUT DATA’ section and update the path to point to the location of your input data on your laptop.
If you wish to modify the specifics of the optimization problem, you can do so by adjusting the values of the parameters in the ‘INPUT PARAMETERS’ section.
Decide whether you want to export a full set of results or just a summary. You can make this selection in the ‘EXPORT SUMMARY’ and ‘FULL EXPORT’ sections. Don’t forget to update the path to the location where you want the exported data to be saved.

**Contributors**
Roxanne Vandenberghe, Gabriele Humbert
