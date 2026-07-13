This project includes modifications to the original Biome-BGCMuSo model (version 6.1). All modified parts of the code have been annotated with comments in the style of "ZL" to indicate the changes made by Lei Zheng. These annotations aim to provide clarity and traceability for the modifications.

<img width="1535" height="863" alt="放牧示意图" src="https://github.com/user-attachments/assets/8ebfa6e4-c44e-4cda-8d01-13b31f397a79" />

Key features

（1）Simulation of winter grazing: The model uses the Standing Dead Biomass (STDB) pool to simulate winter grazing. This feature provides a more realistic simulation of carbon and nutrient fluxes during the non-growing season on the Tibetan Plateau. Our model incorporates the Seligman defoliation equation into the grazing process. This coupling allows a more nuanced representation of the interaction between livestock and vegetation. We have also improved the nitrogen (N) cycle within the grazing process. The model now accounts for all N uptake through livestock feeding and assumes that all N in excreta is derived from ingested plant biomass. The model also distinguishes between urinary and faecal nitrogen. Please refer to our publication: Zheng, L., Pei, Y., Yang, F., Yuan, X., Peng, J., Zhu, B., & Wang, G. (2026). Modeling grazing effects on carbon dynamics in alpine grasslands of the Qinghai-Tibetan Plateau using the Biome-BGCMuSo model. Agriculture, Ecosystems & Environment, 397, 110056.

（2）Improvements to Transient Run Process: We have made improvements to the transient run process in our model. Previously, the transient run was limited to using a single round of meteorological data. Our updated model now allows multiple rounds of meteorological data to be input. An innovative feature of our transient run is the ability to simulate a gradual increase in grazing intensity. Users can now set the model to start with zero grazing intensity and gradually adjust it to the normal first year intensity defined for the simulation. This feature is particularly useful for studying the effects of increasing grazing pressure on ecosystem dynamics and for creating more realistic scenarios of grazing management over time.

（3）Simplified regional simulation configuration: To streamline the setup for regional scale simulations, we have implemented a centralised configuration strategy that distinguishes between spatially invariant and variant parameters. The INI_CONST file now contains parameters that are common to all pixels and do not vary with spatial location. This approach eliminates the need to adjust these common parameters individually in each pixel's INI file.

（4）New phenology model: This project introduces new phenology calculation methods. Please refer to our publication: Zheng, L., Qi, Y., Wang, Y., Peng, J., & Qin, Z. (2022). Calibration and validation of phenological models for Biome-BGCMuSo in the grasslands of Tibetan Plateau using remote sensing data. Agricultural and Forest Meteorology, 322, 109001.









