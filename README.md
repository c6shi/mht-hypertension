---
output:
  html_document: default
  pdf_document: default
---
# Finding the Effect of MHT Initiation on Hypertension Incidence in Menopausal Women Using L-TMLE
PBHLTH 240B Fall 2026 Research Project

We further explore different types of interventions (initiation rules).

# Data:
We use the SWAN data publicly available on ICPSR: https://www.icpsr.umich.edu/web/ICPSR/series/253/studies. 
- Screener: ICPSR 04368
- Baseline: ICPSR 28762
- Visit 01: ICPSR 29221
- Visit 02: ICPSR 29401
- Visit 03: ICPSR 29701
- Visit 04: ICPSR 30142
- Visit 05: ICPSR 30501
- Visit 06: ICPSR 31181
- Visit 07: ICPSR 31901
- Visit 08: ICPSR 32122
- Visit 09: ICPSR 32721
- Visit 10: ICPSR 32961

Data cleaning is performed in `data-cleaning.R`. The final clean dataset can be found as a csv file: `data/mht_final.csv`.

## Data Dictionary:
`STATUS`: Menopausal status
0. Postmenopausal by bilateral salpingo oophorectomy
1. Natural postmenopause (no bleeding in 12 months not due to hysterectomy)
2. Late perimenopause (no bleeding in 3-11 months)
3. Early perimenopause (change in length of bleed or interbleed interval)
4. Pre-menopausal (no change in bleeding patterns)
5. Pregnant/breastfeeding
6. Unknown due to HT use
7. Unknown due to hysterectomy

Bilateral salpingo-oophorectomy (BSO): remove both ovaries & fallopian tubes

Hysterectomy: remove uterus & cervix


`HOTFLAS`: Thinking back over the past two weeks, please circle the number corresponding to how often you experienced any of the following: How often have you had hot flashes or flushes?
0. Not at all
1. 1-5 days
2. 6-8 days
3. 9-13 days
4. Every day