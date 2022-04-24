# What is this repository?
This repository provides the implementation code for the paper, A Deep Generative Approach to Conditional Sampling, published on JASA in 2022.

You can find our paper on [JASA.](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.2016424?journalCode=uasa20)




# What is included in this repository?
The main function to define and train the conditional sampler is in *fitgcde.R*. If you just want to adopt our methods and apply to your own problem, this is the file you should look at.

To replicate the simulation results in the paper, you can find the corresponding code in the folder *simulation* and companion code under folder *utility*, 

# What is the development environment?
The code is developed under R 3.6.1, Tensorflow 2.0.0-gpu. You may need to modify the code accordingly if you are using different environment.

# Citation
Please cite our work in your publications if it helps your research:

    @article{zhou2022deep,
      title={A deep generative approach to conditional sampling},
      author={Zhou, Xingyu and Jiao, Yuling and Liu, Jin and Huang, Jian},
      journal={Journal of the American Statistical Association},
      pages={1--12},
      year={2022},
      publisher={Taylor \& Francis}
    }
