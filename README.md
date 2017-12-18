QMoM: A 0D Aerosol Size Distribution Moment Model
===============================================
Toy box model for tracking aerosol size distribution moments.

Model Description
-----------------
This is a 0 dimensional implementation of the quadrature method of moments model described by McGraw, 1997.
Most operational aerosol microphysics models track the aerosol size distribution directly, either by perscribing a size
distribution to one or more "modes" or by breaking up the radius space into discrete bins. This model, instead of tracking
the distribution itself, tracks the radial moments of a distribution. This technique has advantages:
- QMoM makes no assumptions about the underlying size distribution
- QMoM is less expensive than sectional microphysical models
- Observations of aerosol distributions in the atmosphere often measure moments of that distribution, resulting in simpler comparisons with QMoM than other models
- Radiative properties of aerosols can be calculated directly from moments, facilitating simple integration into coupled aerosol-radiation models

Usage
-----
First, you will need to create all the directories required by the model and scripts:

```
./init.sh
```

To use this model, create or select one of the experiment .ini files and run qmom.py as follows:

```
python qmom.py [experiment_name]
```

To visualize the model output, run visualize.py as follows:

```
python visualize.py [experiment_name]
```

Literature
----------
McGraw, Robert. "Description of aerosol dynamics by the quadrature method of moments." Aerosol Science and Technology 27.2 (1997): 255-265.
