[![Build Status](https://travis-ci.org/TechnionYP5779/SimuMole.svg?branch=master)](https://travis-ci.org/TechnionYP5779/SimuMole)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/TechnionYP5779/SimuMole.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/TechnionYP5779/SimuMole/alerts/)
# Molecular dynamics

### Abstract
Molecular dynamics is a computer simulation method for studying the physical movements of atoms and molecules. The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the dynamic evolution of the system. In the most common version, the trajectories of atoms and molecules are determined by numerically solving Newton's equations of motion for a system of interacting particles, where forces between the particles and their potential energies are often calculated using interatomic potentials or molecular mechanics force fields.

### Client
Rotem Ornai, PHD candidate from the faculty of biotechnology and food engineering under the supervision of Yoval Shoham.
as a part of her research she needs to simulate the movement of molecules in biological processes under certain parameters such as temperature and time using the molecular dynamics algorithm and produce an animation of the process

### Goals
- Create a web application that is able to create an animation of the process with the parameters received from the user
- The user should be able to store and retrieve his simulations and possibly share them
- The user should be able to download the animation
- The user should be able to display an existing simulation that was created by our app

### Technologies
OpenMM: A high performance toolkit for molecular simulation. We use it for running the simulation to get the non-visualize results.

PyMOL: A powerful and comprehensive molecular visualization product for rendering and animating 3D molecular structures. We use it for fetching two PDB files to one containing those two and for displaying the simulations.

Django: An open source web framework written in Python. We use it for developing our web.

