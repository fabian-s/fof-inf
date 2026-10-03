# fastFGEE comparison: results

## A. fastFGEE's simulation design (Loewinger et al., arXiv 2506.20437v3, §4)

Replicates per cell: 199–200. Pointwise 95% coverage pooled over beta0, beta1, beta2 (their Table 3 statistic).

| design | model | fgee_base | fgee_base_wild | fgee_base_fixed | refund_cr_freq_centred | refund_cr | refund_cl2 | reml_cl2 | fgee_independent_independent | fgee_independent_independent_wild | fgee_ar1_exchangeable | fgee_ar1_exchangeable_wild |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| binomial N=25 n_i=1 | 0.740 | 0.686 | 0.959 | 0.868 | 0.868 | 0.895 | 0.915 | 0.916 | 0.773 | 0.904 | 0.770 | 0.907 |
| binomial N=100 n_i=1 | 0.783 | 0.680 | 0.998 | 0.930 | 0.930 | 0.944 | 0.949 | 0.948 | 0.883 | 0.941 | 0.877 | 0.940 |
| binomial N=500 n_i=1 | 0.795 | 0.620 | 1.000 | 0.945 | 0.945 | 0.953 | 0.954 | 0.954 | 0.918 | 0.951 | 0.915 | 0.949 |
| binomial N=25 n_i=5 | 0.541 | 0.764 | 0.963 | 0.897 | 0.897 | 0.904 | 0.922 | 0.922 | 0.790 | 0.959 | 0.800 | 0.961 |
| binomial N=100 n_i=5 | 0.558 | 0.751 | 0.999 | 0.941 | 0.941 | 0.944 | 0.948 | 0.948 | 0.901 | 0.958 | 0.900 | 0.958 |
| binomial N=500 n_i=5 | 0.575 | 0.715 | 1.000 | 0.949 | 0.949 | 0.950 | 0.951 | 0.951 | 0.932 | 0.962 | 0.927 | 0.956 |
| gaussian N=25 n_i=1 | 0.655 | 0.624 | 0.950 | 0.882 | 0.889 | 0.903 | 0.925 | 0.925 | 0.750 | 0.904 | 0.715 | 0.883 |
| gaussian N=100 n_i=1 | 0.694 | 0.631 | 0.994 | 0.934 | 0.935 | 0.944 | 0.949 | 0.949 | 0.873 | 0.931 | 0.855 | 0.922 |
| gaussian N=500 n_i=1 | 0.720 | 0.608 | 1.000 | 0.944 | 0.944 | 0.950 | 0.951 | 0.950 | 0.921 | 0.945 | 0.920 | 0.941 |
| gaussian N=25 n_i=5 | 0.434 | 0.644 | 0.959 | 0.906 | 0.909 | 0.912 | 0.931 | 0.930 | 0.740 | 0.952 | 0.747 | 0.951 |
| gaussian N=100 n_i=5 | 0.453 | 0.632 | 0.996 | 0.939 | 0.940 | 0.942 | 0.946 | 0.946 | 0.879 | 0.944 | 0.857 | 0.944 |
| gaussian N=500 n_i=5 | 0.463 | 0.597 | 1.000 | 0.948 | 0.948 | 0.949 | 0.950 | 0.948 | 0.930 | 0.955 | 0.919 | 0.951 |
| poisson N=25 n_i=1 | 0.735 | 0.772 | 0.979 | 0.891 | 0.891 | 0.904 | 0.932 | 0.931 | 0.877 | 0.919 | 0.873 | 0.919 |
| poisson N=100 n_i=1 | 0.742 | 0.788 | 1.000 | 0.926 | 0.926 | 0.932 | 0.943 | 0.930 | 0.914 | 0.930 | 0.914 | 0.926 |
| poisson N=500 n_i=1 | 0.764 | 0.802 | 1.000 | 0.944 | 0.944 | 0.946 | 0.949 | 0.938 | 0.938 | 0.944 | 0.939 | 0.942 |
| poisson N=25 n_i=5 | 0.493 | 0.844 | 0.937 | 0.908 | 0.908 | 0.911 | 0.934 | 0.934 | 0.891 | 0.939 | 0.893 | 0.951 |
| poisson N=100 n_i=5 | 0.491 | 0.858 | 0.973 | 0.931 | 0.931 | 0.932 | 0.940 | 0.939 | 0.916 | 0.933 | 0.922 | 0.944 |
| poisson N=500 n_i=5 | 0.496 | 0.870 | 0.997 | 0.943 | 0.943 | 0.943 | 0.947 | 0.942 | 0.938 | 0.942 | 0.939 | 0.946 |

Largest Monte Carlo SE of a pooled coverage: 0.010.

Ranges of pooled coverage over cells (n_i = 5 | n_i = 1):

- model: 0.43–0.57 | 0.65–0.79
- fgee_base: 0.60–0.87 | 0.61–0.80
- fgee_base_wild: 0.94–1.00 | 0.95–1.00
- fgee_base_fixed: 0.90–0.95 | 0.87–0.94
- refund_cr_freq_centred: 0.90–0.95 | 0.87–0.94
- refund_cr: 0.90–0.95 | 0.89–0.95
- refund_cl2: 0.92–0.95 | 0.92–0.95
- reml_cl2: 0.92–0.95 | 0.92–0.95
- fgee_independent_independent: 0.74–0.94 | 0.75–0.94
- fgee_independent_independent_wild: 0.93–0.96 | 0.90–0.95
- fgee_ar1_exchangeable: 0.75–0.94 | 0.71–0.94
- fgee_ar1_exchangeable_wild: 0.94–0.96 | 0.88–0.95

Per family, n_i = 5 (their design): fgee_base Gaussian 0.60–0.64, Poisson 0.84–0.87, binomial 0.72–0.76.

### SE calibration mean(SE)/sd(estimate), beta1

| design | model | fgee_base | fgee_base_wild | fgee_base_fixed | refund_cr_freq_centred | refund_cr | refund_cl2 | reml_cl2 | fgee_independent_independent | fgee_independent_independent_wild | fgee_ar1_exchangeable | fgee_ar1_exchangeable_wild |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| binomial N=25 n_i=1 | 0.59 | 0.53 | 0.53 | 0.81 | 0.81 | 0.87 | 0.95 | 0.95 | 0.73 | 0.73 | 0.73 | 0.73 |
| binomial N=100 n_i=1 | 0.64 | 0.51 | 0.51 | 0.95 | 0.95 | 1.00 | 1.02 | 1.02 | 0.91 | 0.91 | 0.91 | 0.91 |
| binomial N=500 n_i=1 | 0.65 | 0.45 | 0.45 | 0.98 | 0.98 | 1.01 | 1.02 | 1.02 | 0.96 | 0.96 | 0.96 | 0.96 |
| binomial N=25 n_i=5 | 0.34 | 0.62 | 0.62 | 0.86 | 0.86 | 0.88 | 0.95 | 0.95 | 0.77 | 0.77 | 0.79 | 0.79 |
| binomial N=100 n_i=5 | 0.37 | 0.60 | 0.60 | 0.98 | 0.98 | 0.99 | 1.01 | 1.01 | 0.94 | 0.94 | 0.95 | 0.95 |
| binomial N=500 n_i=5 | 0.36 | 0.53 | 0.53 | 0.97 | 0.97 | 0.98 | 0.98 | 0.98 | 0.95 | 0.95 | 0.96 | 0.96 |
| gaussian N=25 n_i=1 | 0.49 | 0.45 | 0.45 | 0.84 | 0.86 | 0.90 | 0.99 | 0.99 | 0.74 | 0.74 | 0.66 | 0.66 |
| gaussian N=100 n_i=1 | 0.52 | 0.46 | 0.46 | 0.94 | 0.95 | 0.98 | 1.00 | 1.00 | 0.86 | 0.86 | 0.82 | 0.82 |
| gaussian N=500 n_i=1 | 0.55 | 0.44 | 0.44 | 0.98 | 0.98 | 1.00 | 1.01 | 1.01 | 0.96 | 0.96 | 0.95 | 0.95 |
| gaussian N=25 n_i=5 | 0.27 | 0.48 | 0.48 | 0.88 | 0.90 | 0.90 | 0.99 | 0.99 | 0.68 | 0.68 | 0.65 | 0.65 |
| gaussian N=100 n_i=5 | 0.28 | 0.46 | 0.46 | 0.96 | 0.97 | 0.97 | 0.99 | 0.99 | 0.90 | 0.90 | 0.87 | 0.87 |
| gaussian N=500 n_i=5 | 0.29 | 0.44 | 0.44 | 1.00 | 1.00 | 1.00 | 1.01 | 1.01 | 0.97 | 0.97 | 0.97 | 0.97 |
| poisson N=25 n_i=1 | 0.54 | 0.59 | 0.59 | 0.81 | 0.81 | 0.84 | 0.95 | 0.95 | 0.78 | 0.78 | 0.78 | 0.78 |
| poisson N=100 n_i=1 | 0.58 | 0.64 | 0.64 | 0.93 | 0.93 | 0.94 | 1.00 | 1.00 | 0.91 | 0.91 | 0.91 | 0.91 |
| poisson N=500 n_i=1 | 0.60 | 0.65 | 0.65 | 0.97 | 0.97 | 0.98 | 1.01 | 1.01 | 0.97 | 0.97 | 0.97 | 0.97 |
| poisson N=25 n_i=5 | 0.30 | 0.69 | 0.69 | 0.84 | 0.84 | 0.84 | 0.94 | 0.94 | 0.81 | 0.81 | 0.80 | 0.80 |
| poisson N=100 n_i=5 | 0.31 | 0.76 | 0.76 | 0.95 | 0.95 | 0.95 | 1.00 | 1.00 | 0.92 | 0.92 | 0.92 | 0.92 |
| poisson N=500 n_i=5 | 0.31 | 0.76 | 0.76 | 0.95 | 0.95 | 0.95 | 0.97 | 0.97 | 0.94 | 0.94 | 0.95 | 0.95 |

### SE calibration mean(SE)/sd(estimate), beta2

| design | model | fgee_base | fgee_base_wild | fgee_base_fixed | refund_cr_freq_centred | refund_cr | refund_cl2 | reml_cl2 | fgee_independent_independent | fgee_independent_independent_wild | fgee_ar1_exchangeable | fgee_ar1_exchangeable_wild |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| binomial N=25 n_i=1 | 0.58 | 0.53 | 0.53 | 0.81 | 0.81 | 0.87 | 0.96 | 0.96 | 0.73 | 0.73 | 0.74 | 0.74 |
| binomial N=100 n_i=1 | 0.64 | 0.52 | 0.52 | 0.95 | 0.95 | 1.00 | 1.02 | 1.02 | 0.92 | 0.92 | 0.92 | 0.92 |
| binomial N=500 n_i=1 | 0.67 | 0.47 | 0.47 | 1.00 | 1.00 | 1.03 | 1.04 | 1.04 | 0.99 | 0.99 | 0.99 | 0.99 |
| binomial N=25 n_i=5 | 0.43 | 0.60 | 0.60 | 0.86 | 0.86 | 0.88 | 0.95 | 0.95 | 0.78 | 0.78 | 0.81 | 0.81 |
| binomial N=100 n_i=5 | 0.45 | 0.60 | 0.60 | 0.97 | 0.97 | 0.99 | 1.01 | 1.01 | 0.95 | 0.95 | 0.95 | 0.95 |
| binomial N=500 n_i=5 | 0.48 | 0.59 | 0.59 | 1.05 | 1.05 | 1.06 | 1.06 | 1.06 | 1.04 | 1.04 | 1.00 | 1.00 |
| gaussian N=25 n_i=1 | 0.48 | 0.44 | 0.44 | 0.82 | 0.85 | 0.88 | 0.97 | 0.97 | 0.72 | 0.72 | 0.64 | 0.64 |
| gaussian N=100 n_i=1 | 0.53 | 0.46 | 0.46 | 0.96 | 0.96 | 0.99 | 1.02 | 1.02 | 0.90 | 0.90 | 0.88 | 0.88 |
| gaussian N=500 n_i=1 | 0.56 | 0.44 | 0.44 | 0.99 | 0.99 | 1.01 | 1.02 | 1.02 | 0.96 | 0.96 | 0.96 | 0.96 |
| gaussian N=25 n_i=5 | 0.35 | 0.47 | 0.47 | 0.89 | 0.90 | 0.91 | 0.99 | 0.99 | 0.73 | 0.73 | 0.80 | 0.80 |
| gaussian N=100 n_i=5 | 0.36 | 0.46 | 0.46 | 0.97 | 0.97 | 0.98 | 1.00 | 1.00 | 0.92 | 0.92 | 0.87 | 0.87 |
| gaussian N=500 n_i=5 | 0.37 | 0.43 | 0.43 | 0.98 | 0.98 | 0.99 | 0.99 | 0.99 | 0.97 | 0.97 | 0.94 | 0.94 |
| poisson N=25 n_i=1 | 0.55 | 0.59 | 0.59 | 0.80 | 0.80 | 0.83 | 0.95 | 0.95 | 0.78 | 0.78 | 0.78 | 0.78 |
| poisson N=100 n_i=1 | 0.58 | 0.66 | 0.66 | 0.92 | 0.92 | 0.93 | 0.99 | 0.99 | 0.90 | 0.90 | 0.90 | 0.90 |
| poisson N=500 n_i=1 | 0.61 | 0.70 | 0.70 | 0.99 | 0.99 | 0.99 | 1.01 | 1.01 | 0.98 | 0.98 | 0.98 | 0.98 |
| poisson N=25 n_i=5 | 0.40 | 0.74 | 0.74 | 0.87 | 0.87 | 0.88 | 0.97 | 0.97 | 0.85 | 0.85 | 0.87 | 0.87 |
| poisson N=100 n_i=5 | 0.39 | 0.77 | 0.77 | 0.92 | 0.92 | 0.92 | 0.96 | 0.96 | 0.90 | 0.90 | 0.93 | 0.93 |
| poisson N=500 n_i=5 | 0.40 | 0.83 | 0.83 | 1.00 | 1.00 | 1.01 | 1.02 | 1.02 | 0.99 | 0.99 | 0.99 | 0.99 |

### RMSE relative to the pffr fit (their Table 1 statistic), beta1

| design | fgee_ar1_exchangeable | fgee_independent_independent | reml_cl2 |
|---|---|---|---|
| binomial N=100 n_i=1 | 1.00 | 1.00 | 1.00 |
| binomial N=100 n_i=5 | 0.91 | 0.90 | 1.00 |
| binomial N=25 n_i=1 | 0.97 | 0.98 | 1.01 |
| binomial N=25 n_i=5 | 0.83 | 0.85 | 1.00 |
| binomial N=500 n_i=1 | 1.00 | 0.98 | 1.00 |
| binomial N=500 n_i=5 | 0.92 | 0.92 | 1.00 |
| gaussian N=100 n_i=1 | 1.04 | 1.01 | 1.00 |
| gaussian N=100 n_i=5 | 0.93 | 0.91 | 1.00 |
| gaussian N=25 n_i=1 | 0.96 | 0.94 | 1.00 |
| gaussian N=25 n_i=5 | 0.86 | 0.87 | 1.00 |
| gaussian N=500 n_i=1 | 0.99 | 0.96 | 1.00 |
| gaussian N=500 n_i=5 | 0.95 | 0.91 | 1.00 |
| poisson N=100 n_i=1 | 0.97 | 0.93 | 1.00 |
| poisson N=100 n_i=5 | 0.84 | 0.92 | 1.00 |
| poisson N=25 n_i=1 | 0.95 | 0.94 | 1.00 |
| poisson N=25 n_i=5 | 0.81 | 0.87 | 1.00 |
| poisson N=500 n_i=1 | 0.98 | 0.97 | 1.00 |
| poisson N=500 n_i=5 | 0.84 | 0.96 | 1.00 |

## B. Restart study data sets (M1: ff + z, one curve per subject, core cells at mid signal)

Replicates per cell: 200–200.

### Coverage, beta

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| binomial fpc G=100 | 0.779 | 0.952 | 0.845 | 0.844 | 0.813 |
| binomial fpc G=40 | 0.749 | 0.926 | 0.732 | 0.724 | 0.672 |
| binomial iid G=100 | 0.994 | 0.990 | 0.872 | 0.872 | 0.869 |
| binomial iid G=40 | 0.984 | 0.973 | 0.828 | 0.828 | 0.824 |
| binomial ou G=100 | 0.819 | 0.951 | 0.851 | 0.849 | 0.831 |
| binomial ou G=40 | 0.800 | 0.937 | 0.747 | 0.742 | 0.710 |
| gaussian fpc G=100 | 0.657 | 0.944 | 0.890 | 0.875 | 0.856 |
| gaussian fpc G=40 | 0.615 | 0.926 | 0.833 | 0.793 | 0.707 |
| gaussian iid G=100 | 0.989 | 0.986 | 0.896 | 0.896 | 0.892 |
| gaussian iid G=40 | 0.992 | 0.989 | 0.859 | 0.859 | 0.853 |
| gaussian ou G=100 | 0.672 | 0.940 | 0.885 | 0.876 | 0.846 |
| gaussian ou G=40 | 0.638 | 0.927 | 0.837 | 0.803 | 0.718 |
| poisson fpc G=100 | 0.682 | 0.945 | 0.888 | 0.889 | 0.847 |
| poisson fpc G=40 | 0.658 | 0.928 | 0.841 | 0.836 | 0.701 |
| poisson iid G=100 | 0.989 | 0.986 | 0.891 | 0.891 | 0.890 |
| poisson iid G=40 | 0.993 | 0.989 | 0.855 | 0.855 | 0.851 |
| poisson ou G=100 | 0.698 | 0.942 | 0.888 | 0.887 | 0.857 |
| poisson ou G=40 | 0.684 | 0.927 | 0.837 | 0.833 | 0.743 |

### Coverage, mean

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| binomial fpc G=100 | 0.722 | 0.941 | 0.862 | 0.861 | 0.825 |
| binomial fpc G=40 | 0.704 | 0.921 | 0.819 | 0.814 | 0.714 |
| binomial iid G=100 | 0.963 | 0.949 | 0.872 | 0.872 | 0.870 |
| binomial iid G=40 | 0.952 | 0.925 | 0.828 | 0.828 | 0.827 |
| binomial ou G=100 | 0.752 | 0.940 | 0.855 | 0.854 | 0.836 |
| binomial ou G=40 | 0.747 | 0.927 | 0.808 | 0.805 | 0.761 |
| gaussian fpc G=100 | 0.634 | 0.944 | 0.907 | 0.891 | 0.868 |
| gaussian fpc G=40 | 0.606 | 0.931 | 0.866 | 0.828 | 0.757 |
| gaussian iid G=100 | 0.967 | 0.961 | 0.904 | 0.904 | 0.903 |
| gaussian iid G=40 | 0.968 | 0.958 | 0.868 | 0.868 | 0.865 |
| gaussian ou G=100 | 0.646 | 0.942 | 0.906 | 0.897 | 0.870 |
| gaussian ou G=40 | 0.624 | 0.931 | 0.866 | 0.842 | 0.765 |
| poisson fpc G=100 | 0.653 | 0.945 | 0.904 | 0.902 | 0.864 |
| poisson fpc G=40 | 0.642 | 0.934 | 0.867 | 0.862 | 0.745 |
| poisson iid G=100 | 0.966 | 0.960 | 0.901 | 0.901 | 0.900 |
| poisson iid G=40 | 0.968 | 0.958 | 0.868 | 0.868 | 0.865 |
| poisson ou G=100 | 0.664 | 0.943 | 0.905 | 0.904 | 0.878 |
| poisson ou G=40 | 0.657 | 0.930 | 0.862 | 0.859 | 0.784 |

### Coverage, alpha

| design | reml_model | reml_cl2 | fgee_independent | fgee_independent_wild | fgee_ar1 | fgee_ar1_wild | fgee_fpca | fgee_fpca_wild |
|---|---|---|---|---|---|---|---|---|
| binomial fpc G=100 | 0.702 | 0.920 | 0.846 | 0.939 | 0.850 | 0.940 | 0.814 | 0.910 |
| binomial fpc G=40 | 0.668 | 0.885 | 0.829 | 0.945 | 0.823 | 0.944 | 0.655 | 0.822 |
| binomial iid G=100 | 0.943 | 0.924 | 0.856 | 0.894 | 0.856 | 0.894 | 0.847 | 0.895 |
| binomial iid G=40 | 0.903 | 0.859 | 0.788 | 0.830 | 0.788 | 0.830 | 0.787 | 0.838 |
| binomial ou G=100 | 0.749 | 0.926 | 0.833 | 0.932 | 0.831 | 0.932 | 0.828 | 0.917 |
| binomial ou G=40 | 0.712 | 0.904 | 0.795 | 0.934 | 0.790 | 0.937 | 0.737 | 0.899 |
| gaussian fpc G=100 | 0.641 | 0.950 | 0.921 | 0.940 | 0.913 | 0.935 | 0.891 | 0.913 |
| gaussian fpc G=40 | 0.614 | 0.937 | 0.881 | 0.909 | 0.853 | 0.898 | 0.786 | 0.856 |
| gaussian iid G=100 | 0.956 | 0.956 | 0.911 | 0.931 | 0.912 | 0.932 | 0.910 | 0.931 |
| gaussian iid G=40 | 0.956 | 0.947 | 0.878 | 0.912 | 0.878 | 0.912 | 0.877 | 0.911 |
| gaussian ou G=100 | 0.659 | 0.951 | 0.921 | 0.932 | 0.910 | 0.927 | 0.889 | 0.907 |
| gaussian ou G=40 | 0.632 | 0.940 | 0.885 | 0.913 | 0.864 | 0.901 | 0.791 | 0.861 |
| poisson fpc G=100 | 0.665 | 0.948 | 0.914 | 0.948 | 0.903 | 0.952 | 0.877 | 0.943 |
| poisson fpc G=40 | 0.641 | 0.933 | 0.883 | 0.935 | 0.873 | 0.940 | 0.777 | 0.893 |
| poisson iid G=100 | 0.954 | 0.951 | 0.905 | 0.936 | 0.905 | 0.936 | 0.897 | 0.937 |
| poisson iid G=40 | 0.960 | 0.953 | 0.887 | 0.928 | 0.887 | 0.928 | 0.881 | 0.930 |
| poisson ou G=100 | 0.668 | 0.950 | 0.920 | 0.950 | 0.918 | 0.950 | 0.887 | 0.940 |
| poisson ou G=40 | 0.655 | 0.935 | 0.881 | 0.935 | 0.875 | 0.938 | 0.804 | 0.907 |

### Coverage, gamma

| design | reml_model | reml_cl2 | fgee_independent | fgee_independent_wild | fgee_ar1 | fgee_ar1_wild | fgee_fpca | fgee_fpca_wild |
|---|---|---|---|---|---|---|---|---|
| binomial fpc G=100 | 0.730 | 0.939 | 0.832 | 0.919 | 0.831 | 0.921 | 0.808 | 0.898 |
| binomial fpc G=40 | 0.691 | 0.922 | 0.801 | 0.909 | 0.799 | 0.912 | 0.698 | 0.825 |
| binomial iid G=100 | 0.945 | 0.923 | 0.841 | 0.885 | 0.841 | 0.885 | 0.845 | 0.889 |
| binomial iid G=40 | 0.931 | 0.890 | 0.779 | 0.832 | 0.780 | 0.833 | 0.782 | 0.840 |
| binomial ou G=100 | 0.740 | 0.924 | 0.810 | 0.895 | 0.811 | 0.894 | 0.786 | 0.880 |
| binomial ou G=40 | 0.724 | 0.920 | 0.767 | 0.890 | 0.761 | 0.890 | 0.715 | 0.859 |
| gaussian fpc G=100 | 0.655 | 0.948 | 0.910 | 0.931 | 0.891 | 0.915 | 0.855 | 0.900 |
| gaussian fpc G=40 | 0.621 | 0.946 | 0.877 | 0.916 | 0.839 | 0.891 | 0.761 | 0.856 |
| gaussian iid G=100 | 0.960 | 0.950 | 0.906 | 0.927 | 0.905 | 0.927 | 0.903 | 0.925 |
| gaussian iid G=40 | 0.959 | 0.951 | 0.865 | 0.901 | 0.864 | 0.901 | 0.860 | 0.900 |
| gaussian ou G=100 | 0.663 | 0.938 | 0.900 | 0.917 | 0.891 | 0.912 | 0.865 | 0.895 |
| gaussian ou G=40 | 0.653 | 0.938 | 0.866 | 0.909 | 0.843 | 0.901 | 0.770 | 0.850 |
| poisson fpc G=100 | 0.662 | 0.945 | 0.905 | 0.931 | 0.899 | 0.931 | 0.847 | 0.896 |
| poisson fpc G=40 | 0.647 | 0.943 | 0.865 | 0.917 | 0.855 | 0.912 | 0.733 | 0.835 |
| poisson iid G=100 | 0.963 | 0.955 | 0.900 | 0.926 | 0.900 | 0.926 | 0.901 | 0.926 |
| poisson iid G=40 | 0.955 | 0.948 | 0.862 | 0.904 | 0.862 | 0.904 | 0.857 | 0.899 |
| poisson ou G=100 | 0.677 | 0.937 | 0.902 | 0.920 | 0.898 | 0.921 | 0.872 | 0.906 |
| poisson ou G=40 | 0.667 | 0.937 | 0.856 | 0.905 | 0.850 | 0.905 | 0.770 | 0.855 |

### beta: RMSE ratio to REML

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| binomial fpc G=100 | 0.55 | 0.55 | 0.63 |
| binomial fpc G=40 | 0.69 | 0.69 | 0.88 |
| binomial iid G=100 | 1.16 | 1.16 | 1.16 |
| binomial iid G=40 | 1.49 | 1.49 | 1.45 |
| binomial ou G=100 | 0.62 | 0.61 | 0.66 |
| binomial ou G=40 | 0.73 | 0.69 | 0.79 |
| gaussian fpc G=100 | 0.48 | 0.52 | 0.57 |
| gaussian fpc G=40 | 0.60 | 0.67 | 0.80 |
| gaussian iid G=100 | 0.84 | 0.84 | 0.84 |
| gaussian iid G=40 | 1.03 | 1.03 | 1.02 |
| gaussian ou G=100 | 0.53 | 0.54 | 0.59 |
| gaussian ou G=40 | 0.57 | 0.59 | 0.71 |
| poisson fpc G=100 | 0.47 | 0.47 | 0.57 |
| poisson fpc G=40 | 0.58 | 0.58 | 0.78 |
| poisson iid G=100 | 0.84 | 0.85 | 0.84 |
| poisson iid G=40 | 1.17 | 1.17 | 1.18 |
| poisson ou G=100 | 0.50 | 0.49 | 0.54 |
| poisson ou G=40 | 0.62 | 0.62 | 0.74 |

### beta: interval width ratio to REML + CL2

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| binomial fpc G=100 | 0.34 | 0.35 | 0.37 |
| binomial fpc G=40 | 0.35 | 0.35 | 0.42 |
| binomial iid G=100 | 0.61 | 0.61 | 0.60 |
| binomial iid G=40 | 0.71 | 0.71 | 0.69 |
| binomial ou G=100 | 0.39 | 0.39 | 0.41 |
| binomial ou G=40 | 0.38 | 0.37 | 0.40 |
| gaussian fpc G=100 | 0.38 | 0.40 | 0.40 |
| gaussian fpc G=40 | 0.38 | 0.40 | 0.42 |
| gaussian iid G=100 | 0.52 | 0.52 | 0.51 |
| gaussian iid G=40 | 0.52 | 0.52 | 0.52 |
| gaussian ou G=100 | 0.40 | 0.40 | 0.41 |
| gaussian ou G=40 | 0.37 | 0.36 | 0.38 |
| poisson fpc G=100 | 0.36 | 0.37 | 0.38 |
| poisson fpc G=40 | 0.36 | 0.36 | 0.40 |
| poisson iid G=100 | 0.51 | 0.51 | 0.50 |
| poisson iid G=40 | 0.53 | 0.53 | 0.53 |
| poisson ou G=100 | 0.37 | 0.37 | 0.38 |
| poisson ou G=40 | 0.38 | 0.37 | 0.39 |

### beta: mean(SE)/sd(estimate)

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| binomial fpc G=100 | 0.57 | 0.93 | 0.65 | 0.65 | 0.58 |
| binomial fpc G=40 | 0.50 | 0.83 | 0.46 | 0.45 | 0.41 |
| binomial iid G=100 | 1.50 | 1.43 | 0.76 | 0.76 | 0.76 |
| binomial iid G=40 | 1.50 | 1.38 | 0.62 | 0.62 | 0.62 |
| binomial ou G=100 | 0.62 | 0.94 | 0.67 | 0.67 | 0.64 |
| binomial ou G=40 | 0.57 | 0.86 | 0.49 | 0.49 | 0.47 |
| gaussian fpc G=100 | 0.45 | 0.93 | 0.76 | 0.73 | 0.68 |
| gaussian fpc G=40 | 0.41 | 0.88 | 0.59 | 0.54 | 0.48 |
| gaussian iid G=100 | 1.32 | 1.30 | 0.88 | 0.88 | 0.88 |
| gaussian iid G=40 | 1.38 | 1.36 | 0.74 | 0.74 | 0.74 |
| gaussian ou G=100 | 0.46 | 0.93 | 0.74 | 0.72 | 0.67 |
| gaussian ou G=40 | 0.43 | 0.89 | 0.61 | 0.57 | 0.49 |
| poisson fpc G=100 | 0.47 | 0.93 | 0.75 | 0.76 | 0.65 |
| poisson fpc G=40 | 0.44 | 0.88 | 0.58 | 0.58 | 0.46 |
| poisson iid G=100 | 1.34 | 1.31 | 0.87 | 0.87 | 0.87 |
| poisson iid G=40 | 1.41 | 1.38 | 0.69 | 0.69 | 0.68 |
| poisson ou G=100 | 0.48 | 0.93 | 0.73 | 0.74 | 0.68 |
| poisson ou G=40 | 0.47 | 0.88 | 0.58 | 0.58 | 0.49 |

### mean: RMSE ratio to REML

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| binomial fpc G=100 | 0.83 | 0.83 | 0.85 |
| binomial fpc G=40 | 0.78 | 0.78 | 0.90 |
| binomial iid G=100 | 1.05 | 1.05 | 1.04 |
| binomial iid G=40 | 1.08 | 1.08 | 1.07 |
| binomial ou G=100 | 0.87 | 0.86 | 0.88 |
| binomial ou G=40 | 0.82 | 0.81 | 0.85 |
| gaussian fpc G=100 | 0.85 | 0.89 | 0.89 |
| gaussian fpc G=40 | 0.83 | 0.90 | 0.92 |
| gaussian iid G=100 | 1.00 | 1.00 | 1.00 |
| gaussian iid G=40 | 1.03 | 1.03 | 1.03 |
| gaussian ou G=100 | 0.86 | 0.86 | 0.91 |
| gaussian ou G=40 | 0.83 | 0.84 | 0.92 |
| poisson fpc G=100 | 0.85 | 0.85 | 0.86 |
| poisson fpc G=40 | 0.82 | 0.83 | 0.89 |
| poisson iid G=100 | 1.00 | 1.00 | 1.00 |
| poisson iid G=40 | 1.03 | 1.03 | 1.03 |
| poisson ou G=100 | 0.85 | 0.84 | 0.88 |
| poisson ou G=40 | 0.84 | 0.83 | 0.89 |

### mean: interval width ratio to REML + CL2

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| binomial fpc G=100 | 0.66 | 0.66 | 0.61 |
| binomial fpc G=40 | 0.58 | 0.57 | 0.53 |
| binomial iid G=100 | 0.82 | 0.82 | 0.81 |
| binomial iid G=40 | 0.81 | 0.82 | 0.80 |
| binomial ou G=100 | 0.68 | 0.67 | 0.66 |
| binomial ou G=40 | 0.60 | 0.59 | 0.56 |
| gaussian fpc G=100 | 0.73 | 0.74 | 0.69 |
| gaussian fpc G=40 | 0.66 | 0.66 | 0.58 |
| gaussian iid G=100 | 0.79 | 0.79 | 0.79 |
| gaussian iid G=40 | 0.74 | 0.74 | 0.74 |
| gaussian ou G=100 | 0.74 | 0.73 | 0.71 |
| gaussian ou G=40 | 0.66 | 0.63 | 0.58 |
| poisson fpc G=100 | 0.73 | 0.73 | 0.66 |
| poisson fpc G=40 | 0.66 | 0.65 | 0.54 |
| poisson iid G=100 | 0.79 | 0.79 | 0.79 |
| poisson iid G=40 | 0.75 | 0.75 | 0.75 |
| poisson ou G=100 | 0.73 | 0.73 | 0.70 |
| poisson ou G=40 | 0.67 | 0.66 | 0.59 |

### mean: mean(SE)/sd(estimate)

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| binomial fpc G=100 | 0.56 | 0.98 | 0.83 | 0.83 | 0.74 |
| binomial fpc G=40 | 0.53 | 0.92 | 0.70 | 0.70 | 0.55 |
| binomial iid G=100 | 1.16 | 1.11 | 0.86 | 0.86 | 0.85 |
| binomial iid G=40 | 1.16 | 1.07 | 0.76 | 0.76 | 0.77 |
| binomial ou G=100 | 0.59 | 0.98 | 0.83 | 0.83 | 0.79 |
| binomial ou G=40 | 0.57 | 0.93 | 0.71 | 0.70 | 0.63 |
| gaussian fpc G=100 | 0.47 | 0.99 | 0.91 | 0.87 | 0.81 |
| gaussian fpc G=40 | 0.44 | 0.97 | 0.81 | 0.73 | 0.62 |
| gaussian iid G=100 | 1.10 | 1.09 | 0.91 | 0.91 | 0.91 |
| gaussian iid G=40 | 1.12 | 1.11 | 0.83 | 0.83 | 0.82 |
| gaussian ou G=100 | 0.48 | 0.99 | 0.91 | 0.89 | 0.82 |
| gaussian ou G=40 | 0.45 | 0.97 | 0.81 | 0.76 | 0.63 |
| poisson fpc G=100 | 0.49 | 1.00 | 0.91 | 0.91 | 0.80 |
| poisson fpc G=40 | 0.47 | 0.98 | 0.81 | 0.80 | 0.60 |
| poisson iid G=100 | 1.11 | 1.09 | 0.91 | 0.91 | 0.91 |
| poisson iid G=40 | 1.13 | 1.11 | 0.83 | 0.83 | 0.83 |
| poisson ou G=100 | 0.50 | 1.00 | 0.91 | 0.91 | 0.84 |
| poisson ou G=40 | 0.48 | 0.97 | 0.80 | 0.80 | 0.66 |

Largest Monte Carlo SE of a coverage: 0.018.

## P. Plasmode (public datasets, MID truth, REML residual curves attached)

Replicates per cell: 200–200.

### Coverage, beta

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| dti G=40 | 0.631 | 0.922 | 0.764 | 0.759 | 0.715 |
| dti G=92 | 0.643 | 0.929 | 0.848 | 0.837 | 0.824 |
| ecg_roahd G=100 | 0.644 | 0.944 | 0.815 | 0.831 | 0.829 |
| ecg_roahd G=40 | 0.641 | 0.937 | 0.773 | 0.789 | 0.750 |
| gait_vc G=138 | 0.517 | 0.940 | 0.829 | 0.866 | 0.817 |
| gait_vc G=40 | 0.470 | 0.936 | 0.737 | 0.798 | 0.759 |
| running G=40 | 0.467 | 0.912 | 0.688 | 0.727 | 0.715 |
| running G=90 | 0.515 | 0.922 | 0.746 | 0.755 | 0.700 |

### Coverage, mean

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| dti G=40 | 0.633 | 0.919 | 0.829 | 0.818 | 0.766 |
| dti G=92 | 0.651 | 0.932 | 0.887 | 0.877 | 0.862 |
| ecg_roahd G=100 | 0.757 | 0.953 | 0.901 | 0.899 | 0.906 |
| ecg_roahd G=40 | 0.734 | 0.944 | 0.861 | 0.849 | 0.837 |
| gait_vc G=138 | 0.539 | 0.943 | 0.907 | 0.926 | 0.902 |
| gait_vc G=40 | 0.504 | 0.936 | 0.852 | 0.875 | 0.838 |
| running G=40 | 0.531 | 0.923 | 0.823 | 0.836 | 0.855 |
| running G=90 | 0.575 | 0.934 | 0.878 | 0.878 | 0.898 |

### Coverage, alpha

| design | reml_model | reml_cl2 | fgee_independent | fgee_independent_wild | fgee_ar1 | fgee_ar1_wild | fgee_fpca | fgee_fpca_wild |
|---|---|---|---|---|---|---|---|---|
| dti G=40 | 0.673 | 0.933 | 0.872 | 0.915 | 0.849 | 0.913 | 0.777 | 0.845 |
| dti G=92 | 0.703 | 0.937 | 0.916 | 0.934 | 0.905 | 0.933 | 0.887 | 0.908 |
| ecg_roahd G=100 | 0.718 | 0.955 | 0.897 | 0.967 | 0.883 | 0.956 | 0.873 | 0.948 |
| ecg_roahd G=40 | 0.693 | 0.955 | 0.863 | 1.000 | 0.845 | 1.000 | 0.828 | 1.000 |
| gait_vc G=138 | 0.415 | 0.930 | 0.910 | 0.909 | 0.926 | 0.927 | 0.883 | 0.913 |
| gait_vc G=40 | 0.419 | 0.893 | 0.818 | 0.861 | 0.839 | 0.887 | 0.757 | 0.835 |
| running G=40 | 0.483 | 0.931 | 0.727 | 0.819 | 0.728 | 0.843 | 0.735 | 0.891 |
| running G=90 | 0.533 | 0.939 | 0.845 | 0.887 | 0.792 | 0.886 | 0.816 | 0.925 |

### Coverage, gamma

| design | reml_model | reml_cl2 | fgee_independent | fgee_independent_wild | fgee_ar1 | fgee_ar1_wild | fgee_fpca | fgee_fpca_wild |
|---|---|---|---|---|---|---|---|---|
| dti G=40 | 0.670 | 0.934 | 0.859 | 0.927 | 0.837 | 0.929 | 0.755 | 0.851 |
| dti G=92 | 0.682 | 0.934 | 0.861 | 0.926 | 0.828 | 0.927 | 0.840 | 0.908 |
| ecg_roahd G=100 | 0.790 | 0.961 | 0.908 | 0.946 | 0.908 | 0.951 | 0.915 | 0.953 |
| ecg_roahd G=40 | 0.766 | 0.952 | 0.885 | 0.964 | 0.875 | 0.968 | 0.841 | 0.945 |
| gait_vc G=138 | 0.534 | 0.933 | 0.778 | 0.891 | 0.875 | 0.932 | 0.740 | 0.916 |
| gait_vc G=40 | 0.501 | 0.946 | 0.649 | 0.911 | 0.760 | 0.930 | 0.635 | 0.902 |
| running G=40 | 0.499 | 0.922 | 0.689 | 0.801 | 0.688 | 0.818 | 0.737 | 0.873 |
| running G=90 | 0.550 | 0.934 | 0.790 | 0.857 | 0.777 | 0.855 | 0.767 | 0.917 |

### beta: RMSE ratio to REML

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| dti G=40 | 0.82 | 0.86 | 0.95 |
| dti G=92 | 0.81 | 0.84 | 0.90 |
| ecg_roahd G=100 | 0.67 | 0.54 | 0.50 |
| ecg_roahd G=40 | 0.78 | 0.64 | 0.72 |
| gait_vc G=138 | 0.78 | 0.81 | 0.72 |
| gait_vc G=40 | 0.73 | 0.76 | 0.67 |
| running G=40 | 0.77 | 0.80 | 0.90 |
| running G=90 | 0.74 | 0.76 | 0.93 |

### beta: interval width ratio to REML + CL2

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| dti G=40 | 0.54 | 0.56 | 0.57 |
| dti G=92 | 0.65 | 0.66 | 0.68 |
| ecg_roahd G=100 | 0.42 | 0.38 | 0.38 |
| ecg_roahd G=40 | 0.43 | 0.40 | 0.43 |
| gait_vc G=138 | 0.58 | 0.74 | 0.68 |
| gait_vc G=40 | 0.43 | 0.56 | 0.50 |
| running G=40 | 0.46 | 0.55 | 0.70 |
| running G=90 | 0.53 | 0.61 | 0.79 |

### beta: mean(SE)/sd(estimate)

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| dti G=40 | 0.48 | 0.95 | 0.68 | 0.66 | 0.59 |
| dti G=92 | 0.50 | 0.97 | 0.82 | 0.80 | 0.76 |
| ecg_roahd G=100 | 0.80 | 1.02 | 0.71 | 0.73 | 0.73 |
| ecg_roahd G=40 | 0.72 | 0.98 | 0.61 | 0.62 | 0.58 |
| gait_vc G=138 | 0.42 | 0.98 | 0.80 | 0.87 | 0.83 |
| gait_vc G=40 | 0.35 | 0.96 | 0.65 | 0.76 | 0.68 |
| running G=40 | 0.32 | 0.88 | 0.56 | 0.65 | 0.66 |
| running G=90 | 0.40 | 0.93 | 0.69 | 0.75 | 0.79 |

### mean: RMSE ratio to REML

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| dti G=40 | 0.92 | 0.93 | 1.02 |
| dti G=92 | 0.94 | 0.95 | 1.01 |
| ecg_roahd G=100 | 0.84 | 0.85 | 0.84 |
| ecg_roahd G=40 | 0.83 | 0.84 | 0.85 |
| gait_vc G=138 | 0.98 | 1.02 | 1.33 |
| gait_vc G=40 | 0.93 | 0.98 | 1.12 |
| running G=40 | 0.97 | 1.00 | 1.38 |
| running G=90 | 1.00 | 1.02 | 1.70 |

### mean: interval width ratio to REML + CL2

| design | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|
| dti G=40 | 0.70 | 0.69 | 0.67 |
| dti G=92 | 0.80 | 0.80 | 0.81 |
| ecg_roahd G=100 | 0.70 | 0.70 | 0.70 |
| ecg_roahd G=40 | 0.63 | 0.63 | 0.62 |
| gait_vc G=138 | 0.85 | 0.98 | 1.17 |
| gait_vc G=40 | 0.71 | 0.81 | 0.86 |
| running G=40 | 0.71 | 0.75 | 1.18 |
| running G=90 | 0.83 | 0.84 | 1.65 |

### mean: mean(SE)/sd(estimate)

| design | reml_model | reml_cl2 | fgee_independent | fgee_ar1 | fgee_fpca |
|---|---|---|---|---|---|
| dti G=40 | 0.49 | 0.97 | 0.78 | 0.75 | 0.65 |
| dti G=92 | 0.52 | 0.99 | 0.89 | 0.87 | 0.83 |
| ecg_roahd G=100 | 1.21 | 1.09 | 0.89 | 0.89 | 0.93 |
| ecg_roahd G=40 | 1.07 | 1.07 | 0.78 | 0.76 | 0.74 |
| gait_vc G=138 | 0.45 | 1.01 | 0.91 | 0.97 | 0.92 |
| gait_vc G=40 | 0.39 | 0.99 | 0.78 | 0.84 | 0.76 |
| running G=40 | 0.43 | 0.97 | 0.75 | 0.77 | 0.83 |
| running G=90 | 0.51 | 1.00 | 0.85 | 0.86 | 0.99 |

Largest Monte Carlo SE of a coverage: 0.020.

