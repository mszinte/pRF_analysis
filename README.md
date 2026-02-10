# About
---
This repository contains code for analyzing the different MRI projects in our team. The analysis_code folder contains the main scripts for preprocessing and pRF analysis. Each project also has its own folder with project-specific analysis code.

# Authors (alphabetic order): 
---
Marco BEDINI, Sina KLING, Uriel LASCOMBES, Martin SZINTE


# Task-specific analysis
---
After pRF analysis each project has its own analysis. The project's read me can be found at : 

- RetinoMaps [README.md](RetinoMaps/README.md)
- Amblyo_prf [README.md](amblyo_prf/README.md)
- Amblyo7T_prf [README.md](amblyo7T_prf/README.md)
- centbids [README.md](centbids/README.md)
- nCSF [README.md](nCSF/README.md)
- amsterdam24 [README.md](amsterdam24/README.md)

# Environment Set Up 
---
To install dependencies run the following:
``` pip install -r requirements.txt ```
Afterwards we recommend working with a conda environment. 

Install pycortex: 

``` CFLAGS="-std=c99" pip install git+https://github.com/gallantlab/pycortex.git ```

Install prfpy:  

``` git clone https://github.com/VU-Cog-Sci/prfpy.git ```

cd to folder 

``` python installer.py ```
everything should be fine, so to check that all version are correct: 
``` pip install -r requirements.txt ```
