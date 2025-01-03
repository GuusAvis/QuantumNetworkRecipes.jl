# QuantumNetworkRecipes

QuantumNetworkRecipes can be used to specify the characteristics of a quantum network in the form of a "recipe".
Each recipe needs to be self contained and must fully specify both the physical characteristics of the network and any protocols executed in the network, leaving no room for ambiguity.
Recipes then provide a common language for computer programs to analyze and share information about quantum networks.
Moreover, QuantumNetworkRecipes provides an interface to access network properties and run performance estimates in an implementation-agnostic way.
For example, one could use the same function and recipe to determine the fidelity of quantum states created by a quantum-repeater chain using either analytical methods, simple numerical methods or discrete-event simulation by implementing different methods of that function.

Apart from just providing abstract types and ideas for defining recipes and functions of recipes, this package also contains concrete types and function implementations.
The goal is to provide tools for defining network recipes in a modular way by combining small, simple recipes, such as recipes for network nodes and edges, and to make it easy to reuse code that defines recipes or functions on recipes.
Towards this end, parametric types and functions that dispatch on type parameters are a common design pattern in this package.
QuantumNetworkRecipes contains a number of recipes for nodes (e.g., nodes with a depolarizing quantum memory and a photon source) and edges (e.g., based on optical fiber and a midpoint station where entanglement is heralded).
Moreover, it contains and implements a number of performance-estimate functions (e.g., the secret-key rate achievable on a repeater chain).

Currently, the package is primarily focussed on quantum-repeater chains, which are represented by the `ChainRecipe` type.
This type requires defining recipes for the edges and the nodes of the chain, as well as the protocol executed by the chain.
Right now, the swap-ASAP protocol (with cutoffs) is the primary protocol that is supported.
Its fidelity, entangling rate and secret-key rate can be determined numerically or, in some cases, analytically by the provided functions.

The future plan for this package is to make it easier to define more different types of networks (e.g., networks that are not chains), and provide a much larger library of functions to evaluate the performance of these networks.
In particular, we would like to collect analytical results from the domain of quantum-network performance analysis and implement them here, thereby providing a something of a repository for known analytical expressions.

We note that particular care has been taken here to ensure compatibility with [StochasticAD.jl](https://github.com/gaurav-arya/StochasticAD.jl).
This means that StochasticAD can be used to automatically extract derivatives from the numerical functions included in this package, e.g., to determine the derivative of the secret-key rate of a repeater chain with respect to the locations of its repeater nodes.


## Installation

Currently, QuantumNetworkRecipes has not been added to the central Julia registry, as it is not yet clear to what extent this package will be used and developed further.
However, it can still be installed easily using the github url either by running
```julia
using Pkg; Pkg.add(url="https://github.com/GuusAvis/QuantumNetworkRecipes.jl")
```
in Julia, or by using `]` in the REPL to go to Pkg mode and then typing
```julia
add https://github.com/GuusAvis/QuantumNetworkRecipes.jl
```


## Citation

This package was developed in the process of writing the paper
Optimization of Quantum-Repeater Networks using Stochastic Automatic Differentiation
by Guus Avis and Stefan Krastanov.
Please refer to that paper when discussing or using this software package.

