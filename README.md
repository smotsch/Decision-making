Program *Decision-making*
=========================

**Table of Contents**
- [Installation](#installation)
- [Program execution](#program-execution)
- [Output](#output)
- [Parameters](#parameters)
- [The program](#the-program)

## Installation

Execute the `makefile`:
```bash
	> make
```
The compiler creates the executable file `./DM` and move it to the directory `bin`. By default, the makefile uses the *gfortran* compiler. If you want to use
*ifort*, open the `makefile` and use the instruction written at the bottom. 

## Program execution

To execute the program, run the binary file:
```bash
	> ./DM
```
The program is going to read the parameters of the model in the files:
* `PARAMETER_DM.txt`: for most of the parameters (number of particles, payoff matrix...)
* `PARAMETER_init.txt`: for the initial condition (uniform, Gaussian...)

The parameters are written in external files because we do not need to recompile when we
change one parameter by doing so. The problem with this method is that we have to write
the value of the parameters at a given line. We cannot add a comment line in the file
`PARAMETER_DM.txt`, otherwise the numbering is ruined.

During the execution, the program  `DM` is going to compute the positions
(`X`) and the strategies (`S`) of the particles at each time step. In the
terminal, the program gives some information about the parameters used for the
simulation. At the end, it displays the computation time.

## Output

The program `DM` save the trajectory and the strategy of each particle over
time. At each time step, it creates in the directory 'data' the files:
* `particleX_******`     : x coordinate of the particles
* `particleY_******`     : y coordinate of the particles
* `particleS_******`     : s strategy of the particles,

with `******` a counter of time step. The program writes only if the parameter
`isTrajectorySave` is `True` (line 33 in `PARAMETER_DM.txt`).


## Parameters

* `PARAMETER_MicroVic_flat.txt`
 * `N`          : number of particles
 * `R_coop`     : radius cooperating particle
 * `R_def`      : radius defecting particle
 * `choiceModel` : 5 different rules
 * `alpha`       : intensity jump process
 * `a00...`      : payoff matrix
 * `Threshold`   : use for threshold model
 * `Lx`         : size of the domain in x
 * `Ly`         : same thing in y
 * `dt`         : time step
 * `Time`       : total time simulation
 * `isGrid`     : use or not the trick of the grid
 * `isInitRand` : initialize or not the random
 * `nbSeed`     : seed number for the randomness
 * `isTrajectorySave` : save or note the trajectories
 * `jumpPrint`  : to save only at certain time step
* `PARAMETER_init.txt`
 * `initCondX`     : choice for the initial condition for `X`
 * `Lx,Ly`         : size domain for uniform distribution
 * `xMean`,`yMean`,`xStd`,`yStd` : parameters for the Gaussian
 * `ratioCoop`    : proportion of cooperating particles


## The program

The main program is the file *main_DM.f90*.
It uses different modules (defined in separated files):

| File                              | Description   |
| ----------------------------------|:-------------:|
| `toolkit`                         | contains the usual functions (RandNorm...)
| `input_output_DM`                 | for the input/output (Lecture, FilePrint...)
| `initial_condition_DM`            | to initialize with the proper initial conditions (InitCond)
| `grid_interaction`                | to use the super grid (CellNumber, ListAdd...)
| `interaction_DM`                  | to update `X` and `S`



The architecture of the program is the following:
```bash
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%    main_MicroVic_flat                                                      %%
  %%      -> declaration of variables                                           %%
  %%      -> lecture of parameters ("Lecture")                                  %%
  %%      -> initialization of variables ("InitCond")                           %%
  %%                                                                            %%
  %%      -> loop in time                                                       %%
  %%        A) update strategy S                                                %%
  %%        B) update position X                                                %%
  %%        C) update of the grid (if necessary)                                %%
  %%        D) write trajectories  ("FilePrint")                                %%
  %%      -<                                                                    %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```
There is no Jedi mind tricks in the program. It supposes to be easy understandable.
There might be two points that needs some explanations:
* A structure `PARAM_DM` is used in order to save all the parameters of
the model in only one variable (just called `P`). This structure is defined in
the file `input_output_DM.f90`. This avoid to write 36 arguments each time a
subroutine is called.
* The grid method (also called the **Verlet list** by Wikipedia) consists to
allocate a  number at each particle depending on its position on a **virtual
grid** (`PosGrid`). Therefore, when we look for the neighbors of a particle,
we only have to   check the particles within a nearby square. The time-saving
is spectacular.


## Legal notice information

 This (modest) program is distributed under the GNU GPL license version 2. More
information are available in the file COPYING.txt. For any information or bugs,
please contact me at: `smotsch[at]asu.edu`.
