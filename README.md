# NIAGARA Trial: Mixture Cure Model Analysis

Mixture cure model analysis of event-free survival from the NIAGARA trial
(Powles et al., *NEJM* 2024;391:1773–1786), using reconstructed individual
patient data (rIPD).

**Reference:** Guyot P et al. *BMC Med Res Methodol* 2012;12:9.

---

## Repository structure

```
digitized_coordinates_control.csv    # Digitized KM coordinates, Control arm
digitized_coordinates_durvalumab.csv # Digitized KM coordinates, Durvalumab arm
efs_ripd_excel.csv                   # Pre-computed rIPD (reference)
00_guyot_reconstruction.R            # Step 0: generate efs_ripd_excel.csv
01_main_analysis.R                   # Main analysis (KM, cure model, regression)
02_sensitivity_analyses.R            # Sensitivity analyses (Supp Figs S1-S5)
README.md
```

Digitized coordinates were extracted from the published EFS figure using
WebPlotDigitizer (no header, two columns: time in months, survival in %).
Numbers at risk (3-month intervals, 0–57 months) were read from the published figure.

---

## How to run

Set the working directory to the repository root, then run scripts in order:

```r
source("00_guyot_reconstruction.R")   # generates data/efs_ripd_excel.csv
source("01_main_analysis.R")
source("02_sensitivity_analyses.R")
```

`data/efs_ripd_excel.csv` is provided as a pre-computed reference.
Running `00_guyot_reconstruction.R` will overwrite it with a new reconstruction;
minor differences in individual patient times are expected due to rounding.

Outputs are written to `outputs/` (created automatically).

---

## R packages required

```r
install.packages(c("survival", "ggplot2", "gridExtra", "scales"))
```

Developed under R 4.4.1.

---

## Citation

If you use this code, please cite the primary trial publication:

> Powles T, et al. Perioperative durvalumab with neoadjuvant chemotherapy
> in operable bladder cancer. *N Engl J Med* 2024;391:1773–1786.

and the reconstruction method:

> Guyot P, et al. Enhanced secondary analysis of survival data: reconstructing
> the data from published Kaplan-Meier survival curves.
> *BMC Med Res Methodol* 2012;12:9.
