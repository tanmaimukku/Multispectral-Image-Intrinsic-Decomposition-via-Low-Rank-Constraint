# Multispectral Image Intrinsic Decomposition via Low Rank Constraint 

## Instructions on how to execute the code
- Down load the dataset from https://drive.google.com/file/d/18x-V4Xf_beFr3vVgq1STYb89ezrdFOEJ/view?usp=sharing .
- Place the dataset folder according to the following directory structure as shown below:
  ```
  ├── code
        ├── lriid
        ├── siid
        └── dataset
            ├── data
            ├── GroundTruth
            └── Masks
  ├── Documants
  └── README.md
  ```
- `data` will contain the input data.
- `GroundTruth` will contain ground truth files for reflectance and shading.
- `Masks` will contain the mask files.
- `lriid` folder contains the "main.m" file which runs the final code. To run open the project in MATLAB and then run "main.m" file. It will ask you for the input file name. 'cup', 'train' and 'plane' are the input file names. Input anyone of them without quotes e.g.- plane. The code take around 1 minute to execute. It will display the ground truth and output reflectance and shading images. The LMSE for LRIID and SIID algorithms will be displayed in the console. You have run each file separately.

## Code Explanantion

### main.m
- First we read the input data.
- Then we import the reflectance basis and use PCA to reduce dimensions of basis (low rank constraint).
- Then we use the illumination matrix and reflectance basis to calculate the shading basis.
- Then we calculate the weight matrices using the method given in the paper.
- Then we calculate a initial low rank estimate of the shading image using generic constraint.
- Then we calculate a initial low rank estimate of the reflectance image using generic constraint and the data constraint defined by the previous shading estimate.
- Then we calculate shading image using previous reflectance estimate and reflectance image using previous shading estimate in a loop until reflectance and shading output converge.
- Then reconstruct the shading and reflectance output images.

### Normalized
Takes in 1 input image array and returns another image array which contains normalized cosine distance.
to signify the differences between spectra of pixels in one neighborhood. This distance can be formulated
as dp,q
dp,q tends to approximate to 0 when pixel p and q have same spectra, and depart from 0 when spectra of p and q
are different.

### Resize
Takes in size and input image array and resizes it to the desired size.

### Weightmap
A function to calculate the weight choices based on the normalized distances. It takes the normalized distances, 
α and β as input parameters and returns an array of the corresponding weight choices.
and β are parameters of sigmoid function. To set α and β, We sample 20 value of α within [1000, 10000] and 50
values of β within [10^−5, 10^−2] and choose which perform best.

### LRIID_Specular
LRIID denotes the low rank multispectral image intrinsic decomposition from our algorithm. We compare the computation
results obtained with and without low rank constraint.
An important step to formulate a low rank constraint is to derive the low rank basis. With the help of multispectral
imaging systems as PMIS and CASSI , we can successfully get grouth-truth illumination spectra. Also, there
are plenty of work referring to how to extract illumination
from the image. For example, can be applied in multispectral domain and performs well in implementation. In
order not to complicate our method, We simplify the process of getting normalized illumination spectra Bs by using
the grouth-truth illumination data.

#### Note
All the specular functions of the same consider the specular, highlight part.
