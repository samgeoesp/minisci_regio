# Additional Notes

## Hyperparameter Search Space.

For this work we used a relatively narrow search space for the majority of our model training:

```
ridge :
    alpha :
      - 0.001
      - 0.01
      - 0.1
      - 1
supportvector :
    kernel :
      - rbf
    gamma :
      - auto
    epsilon :
      - 0.01
      - 0.1
      - 1
    C :
      - 0.1
      - 1
      - 10
      - 100
    degree :
      - 1
      - 2
      - 3
randomforestreg :
    n_estimators :
      - 100
```

## Feature Selection

All features selected for individual models can be found within [selected_features.json](selected_features.json) which contains each individual model and the features selected. Each model has the following name structure:

```<model>_<split>_<random_seed>```