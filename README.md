# lattice_matching_oxides

### Executive Summary
The goal of our work is to find lattice-matched oxide films of different crystal types (i.e., rock salt vs perovskite) with high melting points. These constraints provide suitable conditions for high-temperature thermal stability by limiting oxidation, interlayer diffusion, and phase change. To do this, we rely on computed chemical and physical properties from the [**Materials Project Database**](https://materialsproject.org). For more information, please see our peer-reviewed manuscript here:(insert published work). A detailed description of our materials selection process is present in the supplementary section.

### Tutorial
Download our code:
```
git clone https://github.com/sean-mcsherry/lattice_matching_oxides.git
```
Open the interactive notebook [**find_new_films.ipynb**](find_new_films.ipynb), and follow the steps in the tutorial. 


### Requirements
This code uses materials project api 0.24.5 https://github.com/materialsproject/api. The new materials project api is great, but still in the early stages of development. I will try to keep up with their frequent releases, but please be patient! Small changes in the api may break my code. 

There are several other small requirements for this code including: numpy, pandas, matplotlib, math, gcd, functools, chemparse



