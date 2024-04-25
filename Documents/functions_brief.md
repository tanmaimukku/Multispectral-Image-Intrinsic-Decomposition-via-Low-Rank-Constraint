## Normalized
Takes in 1 input image array and returns another image array which contains normalized cosine distance.
to signify the differences between spectra of pixels in one neighborhood. This distance can be formulated
as dp,q
dp,q tends to approximate to 0 when pixel p and q have same spectra, and depart from 0 when spectra of p and q
are different.

## Resize
Takes in size and input image array and resizes it to the desired size.

## Weightmap
A function to calculate the weight choices based on the normalized distances. It takes the normalized distances, 
α and β as input parameters and returns an array of the corresponding weight choices.
and β are parameters of sigmoid function. To set α and β, We sample 20 value of α within [1000, 10000] and 50
values of β within [10^−5, 10^−2] and choose which perform best.

## LRIID_Specular
LRIID denotes the low rank multispectral image intrinsic decomposition from our algorithm. We compare the computation
results obtained with and without low rank constraint.
An important step to formulate a low rank constraint is to derive the low rank basis. With the help of multispectral
imaging systems as PMIS and CASSI , we can successfully get grouth-truth illumination spectra. Also, there
are plenty of work referring to how to extract illumination
from the image. For example, can be applied in multispectral domain and performs well in implementation. In
order not to complicate our method, We simplify the process of getting normalized illumination spectra Bs by using
the grouth-truth illumination data.

### Note
All the specular functions of the same consider the specular, highlight part.