<img src="assets/astex_logo_0.png" alt="drawing" width="180" style="float: right;"/>

# Predicting Regioselectivity of Minisci Reactions with QM/ML Approaches 

![work_explain](assets/github.png)

---

### [Hybrid QM/ML approaches for the prediction of regioselectivity in Minisci C-H Functionalisation Reactions](linkhere)

The models presented within this work utilise both the heteroarene and radical structures to make predictions. Each structure is fed into a bespoke toolkit for generation of QM descriptors. 

## Installation and Set-up


### Create Environment
```bash
conda env create -f env.yml --prefix ./venv
```
*or*
```bash
mamba create -f env.yml -p ./venv -y -c conda-forge
```

### Install qm_descriptors.
To make predictions, QM descriptors are required. This can be done using the [qm_descriptors](https://github.com/samgeoesp/qm_descriptors) package. Set up can be founnd within that package and changes to make this pip installable will be happening in the future.

### Run Example

```bash
python example/run_example.py
```

## HTE Set

The data created for the HTE set can be found in [data/experimental_reactions.json](data/experimental_reactions.json). 

