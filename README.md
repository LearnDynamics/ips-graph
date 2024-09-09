## Learning interacting particle systems on networks
Joint inference of the weight matrices and interaction kernels of interacting particle systems on networks

Authors: Quanjun Lang, Xiong Wang, Fei Lu, Mauro Maggioni, 

Reference paper: [https://arxiv.org/abs/2402.08412](arxiv)


The code implements the numerical tests using the ALS and ORALS algorithms. The tests examine the dependence of their accuracy on three key parameters: sample size, level of observation noise, and the strength of the stochastic force. They also examine the robustness with respect to possible misspecification of the model (i.e., the hypothesis space does not contain the true interaction kernel). 


## Usage

### Run a demonstration
```
typical_example_demo.m    
```
 It shows a typical joint estimation of the weight matrix and the kernel, as well as the estimated trajectories.

### Convergence tests:
```
typical_exp_convergence_M
typical_exp_convergence_observation_noise
typical_exp_convergence_stochastic_force
```
They show the convergence of the estimator as sample size M increases, the level of the observation noise decreases, or the strength of the stochastic force decays. 

### Kuramoto model on a network
```
Kuramoto_multiple_batch.m
```
It shows an application of the algorithms to the Kuramoto model on a network. The test includes both the ground truth inside and outside (model misspecification) of the hypothesis space. 


### Additional testing
```
main_demo_testing.m
```
It includes additional testing and comparison between the ALS and ORALS algorithms. 

### Cite this paper:
```
@article{lang2024interacting,
  title={Interacting Particle Systems on Networks: joint inference of the network and the interaction kernel},
  author={Lang, Quanjun and Wang, Xiong and Lu, Fei and Maggioni, Mauro},
  journal={arXiv preprint arXiv:2402.08412},
  year={2024}
}
```