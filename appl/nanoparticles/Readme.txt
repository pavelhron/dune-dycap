in this folder you can find:

estimation.cc
estimation.conf
fitdataclass.hh
fitmodel.hh
model.hh
reaction_model.hh
tracer.conf
transportparameters.hh

estimation.cc
- main program
- read initial parameters
- create classes for advection-diffusion-reaction solver
- and parameter estimation object

estimation.conf
- parameter file with options
- important

fitdataclass.hh
- read data from data file (in folder data)
- be carefull about the format
- not well documented, but not so important

fitmodel.hh
- compares experimental data with computed values
- to understand everything you need to have a look
  at files in dune/dycap/estimation, but it is rather
  technical and you do not need it

model.hh
- class to solve advection-diffusion-reaction problem
- important is function setParamValue (if you will change parameters you need to estimate, add some,..)
- solution print at the end of estimation
- still only a boilerplate class especially for advection-diffusion

reaction_model.hh
- important for reaction (right side of the equation)
- two main classes:
  Reaction parameter
  - set and get parameters you need in your model
  AdhesionModel
  - in setup get the new parameters
  - in function f compute the right hand side of the equation
  - to add the reaciton to advection-diffusion system, it is multiplied in
    f_adapter by water content
  for solid phase be careful:
      - we compute only the concentration S_np (see report)
      - from this reason it is divide by bulk_density and multiplied by
      	water content
