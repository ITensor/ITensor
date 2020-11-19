Homepage: http://itensor.org/

An efficient and flexible C++ library for performing tensor network calculations.

The foundation of the library is the Intelligent Tensor or ITensor.
Contracting ITensors is no harder than multiplying scalars: matching indices
automatically find each other and contract. This makes it easy to transcribe
tensor network diagrams into correct, efficient code.

Installation instructions can be found in the [INSTALL](INSTALL.md) file.

## Citation

If you use ITensors.jl in your work, for now please cite the [arXiv preprint](https://arxiv.org/abs/2007.14822):

```
@misc{fishman2020itensor,
    title={The \mbox{ITensor} Software Library for Tensor Network Calculations},
    author={Matthew Fishman and Steven R. White and E. Miles Stoudenmire},
    year={2020},
    eprint={2007.14822},
    archivePrefix={arXiv},
    primaryClass={cs.MS}
}
```
