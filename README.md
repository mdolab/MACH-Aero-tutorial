# MACH-Aero Tutorials
[![Documentation Status](https://readthedocs.com/projects/mdolab-mach-aero-tutorial/badge/?version=latest)](https://mdolab-mach-aero-tutorial.readthedocs-hosted.com/en/latest/?badge=latest)

This repository contains step-by-step tutorials for MDO of aircraft configurations with high fidelity (MACH).
The MACH framework is developed by the [MDO Lab](http://mdolab.engin.umich.edu).
It facilitates the aerodynamic design, analysis, and optimization of aircraft.
The tutorial covers the basics needed to optimize the aerodynamic surface of a wing.
It also includes an airfoil optimization example.

## Tutorial documentation
You can either view the [online tutorial](https://mdolab-mach-aero-tutorial.readthedocs-hosted.com) or build the tutorial documentation locally.
To generate the tutorial locally, open a terminal and enter the following commands:

    $ cd mach_aero_tutorials/doc
    $ make html

This generates html files in _build/html/. You can then open _build/html/index.html in your preferred browser to view the tutorial documentation.

## Required Packages to build
- Sphinx (1.6.6)
