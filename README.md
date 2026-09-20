# HMMVLE

HMMVLE (pronounced as "humble" /ˈhʌmbəl/) is a toolkit for the creation, validation and enhancement of Hidden Markov Models (HMMs).



## Installation

Future releases of HMMVLE will include its installation as a Python package.
However, you can install it now from the source code as following:

```bash
mamba install -c bioconda python pandas pyfamsa pyhmmer biopython -y
mamba install matplotlib scikit-learn plotly -y
pip install --user pytrimal nbformat
```



## Getting started

HMMVLE improves existing HMM models by using the HMM-ModE methodology (Srivastava et al., 2007; Sinha & Lynn 2014).

An example Jupyter notebook can be found in the [examples folder](examples/HMMVLE-ModE.ipynb).



## Citation

The HMMVLE package can be cited as:

Climent Gargallo, G. (2026). HMMVLE: Hidden markov model VaLidation and Enhancement (Version 0.0.1) [Computer software]. https://github.com/guillecg/hmmvle



## References

Srivastava, P. K., Desai, D. K., Nandi, S., & Lynn, A. M. (2007). HMM-ModE - Improved classification using profile hidden Markov models by optimising the discrimination threshold and modifying emission probabilities with negative training sequences. BMC Bioinformatics, 8(1), 104. https://doi.org/10.1186/1471-2105-8-104

Sinha, S., & Lynn, A. (2014). HMM-ModE: Implementation, benchmarking and validation with HMMER3. BMC Research Notes, 7(1), 483. https://doi.org/10.1186/1756-0500-7-483



## Funding

This project has received funding from the European Union’s Horizon Europe Research and Innovation programme under the Marie Skłodowska-Curie Grant Agreement No. 101073271 - [SHINE project](https://www.shine-edn.eu/).
