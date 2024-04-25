

# Project Proposal - Multispectral Image Intrinsic Decomposition via Low Rank Constraint

### Project ID : 21

### Repository link : https://github.com/Digital-Image-Processing-IIITH/project-random

### Team Members
- Amitesh Singh, 20171131 
- Tanmai Mukku, 20171145 
- Mehul Gupta, 20171156 
- Ankitha Eravelli, 2019900009

## Project Goals

- To define better and address the Intrinsic Image Decomposition problem of a whole multispectral image captured under general spectral illumination by propposing a low rank constraint.
- To decompose the shading and reflectance from a single multispectral image using Low Rank Multispectral Image Intrinsic Decomposition model (LRIID).
- To extend the Retinex model for multispectral domain for LRIID, assuming that the basis of Retinex theory would continue to take effect on multispectral domain.

## Problem Definition
Due to the complex illumination and the geometry structure of natural scenes, the spectra curves of the same surface can look very different and thus hinder the effectiveness of the clues that are obtainable from multispectral images. In this project, a Low Rank Multispectral Image Intrinsic Decomposition Model (LRIID) is presented to decompose the shading and reflectance from a single multispectral image. We extend the Retinex model, which is proposed for RGB image intrinsic decomposition, for multispectral domain. Based on this, a low rank constraint is proposed to reduce the ill-posedness of the problem and make the algorithm solvable. We will evaluate the method implemented using the dataset of 12 images with ground truths of shading and reflectance.

We assume the object surface as Lambertian and hence has diffuse reflection. In most prior work on intrinsic image decomposition, the captured luminance spectrum at every point lp is modelled as the product of Lambertian reflectance spectrum rp and shading spectrum sp, where sp   is used to characterize the combined effect of object geometry, illumination, occlusion and shadowing. Mathematically, this model can be expressed as: <br />
lp=sp.∗rp <br />
Where lp, rp and sp are all vectors with dimensions equal to the number of spectral bands of the captured image, .∗ denotes element-wise multiplication. The problem is to derive sp and rp from observed multispectral luminance vector lp. In this project, we will focus on recovering the reflectance spectrum using this model. Once rp is determined, the shading image can be derived by pointwise division. Different from the conventional approaches which operate in the logarithmic domain, we directly formulate the problem in the image domain.

## Project Results (Expected)
- For better visualization the result will be shown in pseudo-rgb and the image will be linearly normalized to the range [0, 1]. 
- We first show the performance of our proposed algorithm. 
- This is followed by our method on on-line dataset and visual results. 
- Finally, we test on our dataset with ground truth and compare our method with the results of an already tested method in [12] from the paper.

## Project Milestones and Expected Timeline
- Understanding the methodology and constraints given in the paper - 21 October
- Acquiring the dataset - 23 October
- Implementing steps 1, 2 and 3 of the algorithm - 30 October
- Implementing steps 4, 5 and 6 of the algorithm - 7 November
- Comparing result of our approach with other approaches - 10 November
- Complete documentation - 13 November
- Complete presentation - 17 November

## Dataset
- To test the algorithm we need a multispectral image dataset.
- We will try to acquire the 12 image dataset created by the authors of the paper. This is not available online.
- If we are not able to acquire the above dataset we will use Nayar multispectral image database. This is easily avalable on the net.
- We will also use some images available on the cornell website.
