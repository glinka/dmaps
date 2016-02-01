Basic implementation of DMAPS algorithm in C++ and Python

===============================
Using the Python module
===============================

The module's methods are very basic, and can be used as follows

```
>>> from test_dmaps import gen_swissroll
>>> swissroll_data = gen_swissroll()
>>> k = 15; epsilon = 2.5
>>> eigvals, eigvects = dmaps.embed_data(data, k, epsilon)
>>> from plot_dmaps import plot_embeddings
>>> plot_embeddings(eigvects, eigvals, k=3)
```

where k is the number of eigenvector/eigenvalue pairs to compute, or equivalently, the dimension of the embedding. If each data vector *x* is *n*-dimensional, the preceeding code would embed into *k=15* dimensions, and the new coordinates for the first data point would be *v1(1), v2(1), v3(1),...,v15(1)* where *v i(j)* represents the *jth* component of the *ith* eigenvector (glossing over the optional dependence on the eigenvalues).

**Note:** It is highly recommended to begin any analysis by running

```
>>> epsilon_plot(epsilons, data)
```

to guide the selection of a proper epsilon value (chosen in the linearly increasing region of resulting figure).
