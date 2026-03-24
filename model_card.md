---
# Doc / guide: https://huggingface.co/docs/hub/model-cards
---

# Model Card for CKB Cirrhosis Score (CCS)

The CKB Cirrhosis Score (CCS) is a  risk prediction model for estimating the probability of cirrhosis using eight easily accessible predictors. It was developed in the China Kadoorie Biobank (CKB) and externally validated in Mass General Hospital (MHG) and NHANES. A benchmark ensemble SuperLearner (eSL) model was also developed for comparison; CCS achieved similar discrimination with greater simplicity and interpretability.

## Model Details

### Model Description

CCS is a penalized multivariable logistic regression model for prediction of cirrhosis, defined primarily by liver stiffness measurement (LSM) ≥13.5 kPa. The final model includes age, sex, body mass index (BMI), diabetes, total cholesterol (TC), high-density lipoprotein cholesterol (HDL-C), aspartate transaminase (AST), γ-glutamyltransferase (GGT), and an interaction between sex and BMI.

The model was designed for practical use in primary care and population-based risk stratification using variables that are widely available in routine clinical settings. In parallel, an ensemble SuperLearner model was built using multiple machine-learning base learners to assess whether substantially better performance could be achieved with a more complex modeling strategy. The simpler CCS performed comparably to the eSL model.

## Uses

### Direct Use

CCS is intended for direct estimation of an individual’s probability of prevalent cirrhosis using routinely measured clinical variables:

- age  
- sex  
- BMI  
- diabetes status  
- TC  
- HDL-C  
- AST  
- GGT  

Potential direct uses include:

- population-level screening and risk stratification
- prioritizing follow-up liver assessment in primary care or large cohorts
- supporting epidemiologic studies of cirrhosis risk

### Downstream Use

Potential downstream uses include:

- embedding CCS into web tools or electronic clinical decision support systems
- use as a baseline comparator in future liver disease prediction studies
- enrichment of higher-risk participants for imaging, elastography, or hepatology referral
- use as an upstream triage component in broader liver health screening workflows

### Out-of-Scope Use

CCS is not intended for:

- replacing clinical diagnosis of cirrhosis or specialist evaluation
- use as the sole basis for treatment decisions
- use in pediatric populations
- use in populations with characteristics very different from the development/validation cohorts without additional validation
- predicting outcomes unrelated to liver disease
- use when required predictor measurements are missing, measured inconsistently, or obtained under incompatible laboratory standards

## Bias, Risks, and Limitations

First, the model was developed for cirrhosis defined by liver stiffness thresholds, primarily LSM ≥13.5 kPa, rather than biopsy-confirmed cirrhosis. Performance may vary when alternative clinical definitions are used, although sensitivity analyses using other thresholds showed broadly similar results.

Second, predictor effects were estimated from adult cohorts with different sex distributions, age structures, and metabolic profiles. Model behavior may therefore differ in underrepresented subgroups or in clinical settings with substantially different disease etiologies.

### Recommendations

Users should:

- validate and, where needed, recalibrate CCS before deployment in new populations
- use CCS as a triage or risk-stratification aid rather than a stand-alone diagnostic tool
- consider local cirrhosis prevalence when interpreting PPV and NPV
- ensure consistent units and laboratory measurement conventions
- avoid extrapolating use beyond adult populations and settings similar to the study cohorts without further evidence

## How to Get Started with the Model

Use the linear predictor below and convert it to a probability:

```python
import math

def ccs_probability(age, female, diabetes, bmi, tc, hdl, ast, ggt):
    def pos_cube(x):
        return max(x, 0) ** 3

    lp = (
        -5.239626
        + 0.046680734 * age
        - 5.6458036 * female
        + 0.46411484 * diabetes
        - 0.16098654 * bmi
        + 0.0026687438 * pos_cube(bmi - 20.200001)
        - 0.0049893916 * pos_cube(bmi - 24.200001)
        + 0.0023206478 * pos_cube(bmi - 28.799999)
        + female * (
            0.26032348 * bmi
            - 0.0031969233 * pos_cube(bmi - 20.200001)
            + 0.0059768578 * pos_cube(bmi - 24.200001)
            - 0.0027799345 * pos_cube(bmi - 28.799999)
        )
        - 0.1056094 * tc
        - 0.40231523 * hdl
        + 0.013602785 * ast
        + 0.10879714 * ggt
        - 0.0001017148 * pos_cube(ggt - 13.4)
        + 0.00013169036 * pos_cube(ggt - 22.3)
        - 2.9975555e-5 * pos_cube(ggt - 52.5)
    )
    return 1 / (1 + math.exp(-lp))
```

## Training Details

### Training Data

The main model was developed in the CKB 3rd resurvey, including 23,187 participants. Candidate predictors were chosen a priori based on known or probable cirrhosis risk factors, clinical rationale, primary care usability, and data availabilityy.

### Training Procedure

#### Training Hyperparameters

- **Training regime:** Penalized maximum likelihood estimation for logistic regression; spline terms with up to 2 degrees of freedom for selected continuous predictors; tuning parameters selected by grid search using AIC

#### Speeds, Sizes, Times

CCS is a lightweight tabular risk model with minimal computational burden.

## Evaluation

### Testing Data, Factors & Metrics

#### Testing Data

External validation was conducted in:
- MHG validation cohort: 167,726 participants
- NHANES validation cohort: 3,132 participants
  
#### Factors

Performance was examined overall and across subgroups.
- age groups
- sex
- alcohol consumption
- diabetes status
- fatty liver status

#### Metrics

- AUROC with 95% confidence intervals
- negative predictive value (NPV)
- calibration intercept and slope
- integrated calibration index (ICI)
- E50 and E90
- continuous net reclassification improvement (NRI)
- decision-curve net benefit

### Results

Internal validation in CKB derivation cohort
- CCS 10-fold cross-validated AUROC: 0.762 (95% CI 0.744, 0.780)
- eSL 10-fold cross-validated AUROC: 0.760 (95% CI 0.742, 0.777)
- CCS NPV: 98.9%
- CCS subgroup AUROC range: 0.732–0.777

External validation in MHG
- CCS AUROC: 0.816 (95% CI 0.799, 0.834)
- eSL AUROC: 0.815 (95% CI 0.796, 0.833)
- LiverRisk AUROC: 0.757 (95% CI 0.735, 0.779)
- APRI AUROC: 0.757 (95% CI 0.735, 0.779)
- Forn’s score AUROC: 0.710 (95% CI 0.688, 0.731)
- FIB-4 AUROC: 0.684 (95% CI 0.661, 0.708)
- CCS NPV: 99.8%

External validation in NHANES
- CCS AUROC: 0.846 (95% CI 0.784, 0.905)
- eSL AUROC: 0.825 (95% CI 0.753, 0.897)
- NFS AUROC: 0.798 (95% CI 0.725, 0.861)
- refitted LiverRisk AUROC: 0.804 (95% CI 0.736, 0.866)
- CCS NPV: 99.4%

#### Summary

CCS demonstrated consistent and strong discrimination for prevalent cirrhosis across internal and external validation cohorts, with AUROC values ranging from 0.762 in internal validation to 0.846 in NHANES. It outperformed or matched conventional non-invasive fibrosis scores and performed similarly to a more complex ensemble SuperLearner benchmark.
