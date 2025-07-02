# MATLAB Simulation of Leaky Integrate-and-Fire (LIF) Neurons

## Overview

This project implements a biologically inspired simulation of Leaky Integrate-and-Fire (LIF) neurons using MATLAB. It serves as a lightweight and fast alternative to traditional Cadence circuit simulations, offering significant speedup with acceptable accuracy.

## Motivation

Conventional neuromorphic circuit simulations using tools like **Cadence** are highly accurate but computationally expensive, often requiring **up to 24 hours** to simulate a single neuron. This project leverages MATLAB to achieve the same simulation goals in **just a few minutes**.

## Key Features

- ✅ **Biologically Inspired**: Based on the LIF model, mimicking neuronal charge-discharge behavior.
- ⚡ **Efficient**: Reduces simulation time from ~1 day to minutes.
- 🎯 **Accurate**: Achieves a **9.2% RMSPE (Root Mean Square Percentage Error)** compared to Cadence.
- 📈 **Data-Driven**: Uses interpolated data from Cadence to model frequency and power behavior.
- 🔧 **Customizable**: Spike shape and response frequency are configurable using polynomial fitting.

## Project Structure

/lif-neuron-matlab
├── lif_neuron_simulation.m # Main script for simulation
├── interpolation_module.m # Handles data interpolation (Cadence to MATLAB)
├── spike_generator.m # Spike waveform generator (exponential)
├── data/
│ └── cadence_data.csv # Example input from Cadence
└── README.md
## How to Use

1. Clone or download this repository.
2. Open the project folder in MATLAB.
3. Run `lif_neuron_simulation.m`.
4. Observe generated spike trains and performance metrics.
5. Compare results to Cadence output using included visualization tools.

## Results

- Simulated spikes in MATLAB closely match those from Cadence.
- For polynomial order 15 interpolation:
  - **RMSPE**: **9.2%**
- Allows real-time analysis of neuron models that would otherwise take hours to simulate.

## Authors

- **Edouard David**  
- **Marc Zhan**

## Supervisors

- M. **Aziz BenLarbi-DELAI**  
- M. **Pietro Ferreira**

**Laboratoire de Génie Électrique de Paris (GEEPS)**  
Sorbonne Université / CentraleSupélec

## License

For academic and research use only. Contact the authors for redistribution or commercial use.
