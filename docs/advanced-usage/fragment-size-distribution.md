# Fragment size distribution

Fragment size distribution $f_{i,j}$ is an integral part of the model. It defines how polymer fragmented from one size class is distributed to smaller size classes. When compared to the fragmentation rate $k_\text{frag}$:

$k_{\text{frag},i}$
: The fragmentation rate $k_\text{frag},i$ is the rate at which the mass of polymer in size class $i$ is lost from that size class due to fragmentation. For a given period of time $\Delta t$, a mass $m = \Delta t k_{\text{frag},i}$ of polymer fragments from size class $i$. However, this doesn't tell us how much of this mass goes to each daughter (smaller) size class. This information is given by $f_{i,j}$.

$f_{i,j}$
: The fragment size distribution is the fraction of this fragmented mass that is transfered from size class $i$ to size class $k$.


## General fragment size distribution matrix

$f_{i,j}$ is a 2D matrix, but as fragmentation only happens from larger to smaller size classes, the matrix is zero everywhere except its lower left triangle. Taking the example of four size classes, with the smallest indexed by 1 and the largest by 4:

![Image showing the 4x4 f_i_j matrix that is all zeros except the lower left triangle](../img/fragment-size-distribution.png)


## Example - only fragmentation from largest size class

Taking the example of a situation where only the largest size class fragments, the matrix can be simplified to be mostly zeros:

![Image showing the 4x4 f_i_j matrix that is all zeroes except the final row](../img/fragment-size-distribution_largest-only.png)

:::{note}
As shown here, the model uses a $f_{i,j}$ matrix with zeros on the diagonal (fragmentation from and to the same size class). For use as a proper mass balance matrix, the diagonal should instead be the fraction of polymer left in that size class after fragmentation.
:::


## Changing the fragment size distribution in the model

The model presumes an even split between all daughter size classes, which is calculated internally by:

```python
# Start with a zero-filled array of shape (N,N)
fsd = np.zeros((n_size_classes, n_size_classes))
# Fill with an equal split between daughter size classes
for i in np.arange(n_size_classes):
    fsd[i, :] = 1 / i if i != 0 else 0
# Get the lower triangle of this matrix, which effectively sets fsd
# to zero for size classes larger or equal to the current one
fsd = np.tril(fsd, k=-1)
```

:::{note}
This default of an even split between size classes is likely to change as our project develops and experimental data gives us a better idea of a sensible default (which might be polymer-specific).
:::

To change this, you must generate your own $f_{i,j}$ matrix and directly set the `fsd` variable. For example:

```python
fmnp = FragmentMNP(config, data)
fmnp.fsd = my_fsd_matrix
```

:::{caution}
Currently, no validation is performed on this matrix, and so the model will produce erroneous results if a matrix in which rows don't sum to unity is provided.
:::