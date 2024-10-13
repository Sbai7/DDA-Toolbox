# The Double Delineation Approach Toolbox
The double delineation approach (DDA) is a new concept in groundwater modelling that was developed in the following couple of papers by Dr. M. A. Sbai. The primary focus was the delineation of capture zones in groundwater flow systems, however the developed technique goes far beyond that and has many other applications as briefly described in the following:

![Alt text]figures/fig1.png?raw=true "")

## Paper 1
[Sbai, M.A. (2018). A Practical Grid-Based Alternative Method to Advective Particle Tracking. Groundwater, 56(6), 881-892. https://doi.org/10.1111/gwat.12646](https://www.researchgate.net/publication/323339565_A_Practical_Grid-Based_Alternative_Method_to_Advective_Particle_Tracking)

### Objective
Proposes a grid-based alternative to the conventional advective particle tracking method used in groundwater flow simulations.

### Motivation
Conventional particle tracking methods are not always robust or practical for general applications, particularly for visualizing capture and swept zones around wells.

### Methodology
- The proposed method uses industry-standard finite-difference grid-based technique as an alternative. 
- It models the forward and backward travel times and residence times on the same grid used for flow solutions. 
- This approach is computationally efficient and easily interpretable for practitioners.

### Advantages
- It accurately delineates capture zones and swept zones around wells, with minimal computational effort.
- The method is more comprehensive and easier to use compared to particle tracking, especially in complex multi-layered aquifers.
- The model introduces a dispersion indicator to assess the quality of numerical results, improving reliability.

### Applications
- The method is tested in two examples: a synthetic pump-and-treat system and a geothermal doublet system in the Paris Basin (Dogger aquifer).
- It demonstrates advantages in delineating capture zones and understanding flow dynamics in groundwater systems.

### Results 
- The grid-based method produces smoother, more reliable results compared to particle-tracking methods, especially in heterogeneous media.
- It also allows the visualization of hydraulic connections and residence times in a clearer, more quantitative manner.

### Conclusion
- The method is a promising alternative to particle tracking, particularly for complex groundwater modelling studies.
- It can be applied to other transport dynamics, such as seawater intrusion and aquifer-river interactions, making it a valuable screening tool before more complex simulations.


## Paper 2
[Sbai, M.A. (2020). A Dual Delineation Approach to Untangle Interferences in Well Doublets Systems. Groundwater, 58(6), 894-902. https://doi.org/10.1111/gwat.12978](https://www.researchgate.net/publication/338175524_A_Dual_Delineation_Approach_to_Untangle_Interferences_in_Well_Doublets_Systems)

### Objective 
- The paper extends the DDA, which was initially introduced to delineate groundwater capture zones around wells. This extended approach helps to untangle interferences in well doublet systems, such as those in geothermal reservoirs, by delineating groundwater contribution areas.
- The method is particularly useful in managing low-enthalpy geothermal resources, helping to assess cold-water breakthroughs and temperature declines at production wells.

### Methodology
- Combines the forward and backward travel times of groundwater, allowing the simultaneous delineation of steady-state capture zones for each well.
- Provides an alternative to the traditional particle tracking method, offering both groundwater contribution areas and capture zones.
- The DDA solves a series of stationary advection equations to monitor the evolution of interference zones between well doublets in a field, such as those in the Dogger aquifer.
- Produces time-related capture zones for pumping wells and injection wells.
- Produces source-sink connections that map out how groundwater volumes move between source (injection) and sink (pumping) wells.
- The approach is computationally efficient, solving for capture zones and interference areas between well pairs through simple post-processing of dual outputs.

### Applications
#### Geothermal Reservoir Management
- The method forecasts the onset and rate of heat breakthroughs (cold water reaching production wells) in low-enthalpy geothermal systems.
- It assesses the connectivity between injection and production wells, quantifying groundwater volumes and flow rates moving between these wells over time.
- The study applies the DD approach to a geothermal well field in the Dogger aquifer (Paris Basin, France), providing insights into the temperature decline at production wells and the time-dependent evolution of flow rate fractions between wells.
#### Tracer Test Interpretation
- The method helps interpret slug tracer tests, which are often used to assess flow rates and cooling rates in geothermal systems.
- By analyzing flow paths and volumetric partitions, it can forecast cooling rates at pumping wells and assist in optimizing well placements and flow rates.

### Results
- The DDA accurately predicted the areas where cold water would break through at production wells at specific times.
- It showed the interconnected flow paths between injection and production wells, and the corresponding flow rate partitions over time.
- The results indicated that most geothermal production wells would experience significant temperature declines by 30 years, which would impact the economic viability of the project.
- The method predicts the rate at which temperatures at production wells will decline, helping in the sustainable management of geothermal systems by identifying which wells are most likely to contribute to cooling.
- It estimates thermal breakthrough times, providing a useful surrogate to complex coupled heat-transport simulations.

### Conclusions
- The DD approach provides a more detailed and comprehensive understanding of groundwater flow dynamics in complex systems like geothermal doublets.
- It is a powerful tool for forecasting temperature declines, understanding flow interferences between wells, and optimizing geothermal resource management.
- The approach can be easily implemented in existing groundwater models and extended to broader applications, including shallow aquifer systems, recharge zones, and other water management challenges.


## Getting started 
- Copy the contents of the downloaded zip file to any folder on your drive or home directory.
- Add this folder and all its subdirectories to MATLAB’s path. To do this, click on the "Set Path" icon in the "HOME" toolbar at the top of MATLAB’s Command Window. This will open the "Set Path" dialog box. Next, click the "Add with Subfolders..." button and select the folder where you installed the DDA package. Finally, click the "Save" button, followed by the "Close" button. 
- Navigate to the "examples" subdirectory within the installation folder and review the provided examples. Each example is available in two formats: (1) a plain text .m MATLAB script file, and (2) a live script file, where not only the descriptions but also the text and graphical results are precomputed. You can rerun these live scripts to interactively solve the examples.


## Citing

Cite the following publications if you got some results using this toolbox:

```
@article{sbai-2018,
  author  = {M. A. Sbai},
  title   = {{A Practical Grid-Based Alternative Method to Advective Particle Tracking}},
  journal = {Groundwater},
  volume  = {56},
  number  = {6},
  pages   = {881-892},
  year    = {2018},
  doi     = {https://doi.org/10.1111/gwat.12646},
  url     = {https://ngwa.onlinelibrary.wiley.com/doi/10.1111/gwat.12646},
}
```

```
@article{sbai-2020,
  author  = {M. A. Sbai},
  title   = {{A Dual Delineation Approach to Untangle Interferences in Well Doublets Systems}},
  journal = {Groundwater},
  volume  = {58},
  number  = {5},
  pages   = 822-830},
  year    = {2020},
  doi     = {https://doi.org/10.1111/gwat.12978},
  url     = {https://ngwa.onlinelibrary.wiley.com/doi/10.1111/gwat.12978},
}
```

