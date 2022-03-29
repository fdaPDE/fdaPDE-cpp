# Welcome to fdaPDE
functional data analysis with partial differential equation regularization.
## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.


| Method | Description | 
| ------------ | ------------- | 
| Content Cell | Content Cell  | 
| Content Cell | Content Cell  | 

!!! note
	use a derivede class to iterate on the elements

!!! important
	important advice


We try to minimize the following functional
$$J_\lambda(\beta,f) = \sum_{i=0}^n |z_i - w_i^T \beta - f(p_i)|^2 - \lambda \int_\Omega(Lf - u)^2$$

