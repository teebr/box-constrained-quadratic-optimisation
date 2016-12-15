# box-constrained-quadratic-optimisation
solve quadratic minimisation subject to inequality constraints

This is an abstract class to allow you to define a generic quadratic minimisation problem with upper and lower limits on the element values. When the `qpBox` class is inherited, you specify two functions (`hCalc` and `fCalc`) which calculate the H and f matrices based on whatever inputs (and additional properties) you choose.

Since it is a `matlab.System` class, it can be run in MATLAB or Simulink (it was originally written to solve a control problem in Simulink where the H and f matrices change over time): examples are provided for both. the built-in MATLAB solver `quadprog` cannot be implemented in Simulink, and also has a lot of overhead so will typically run slower.
