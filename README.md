# pteroparamyxo-BRT

Code for BRT (boosted regression trees) model of paramyxovirus hosts in Pteropodidae. This model generates host suitability predictions for a family of ~200 bat species, as part of a larger study (<a href="https://doi.org/10.1101/2025.09.11.675601" target="_blank">Juman et al. 2025</a>).

## Model Details

### Model Description

We trained a BRT model via the gbm package (Greenwell et al. 2022) to classify positives from unsampled or negative species. This is a pattern recognition, hypothesis-generating approach that identifies predictors associated with a binary response variable (i.e., paramyxovirus PCR positivity) and extrapolates those patterns to negative or unsampled species to predict their likelihood of viral positivity. We train our BRTs on 80% of the complete dataset, using 25 splits, 3000 trees with a Bernoulli error distribution, and five-fold cross-validation to prevent overfitting (Elith et al. 2008). We calculate the corrected area under the curve (AUC) as a performance metric on both our training and test dataset using target shuffling methods (Han et al. 2016). Model parameters are selected based on comparisons of the test AUC, sensitivity, and specificity of models across a grid search of factorial combinations of hyperparameters. Each variable is assigned a score indicating its relative contribution to classification, with higher values indicating stronger influence on the model and the sum of all relative influence scores adding up to 100. The contribution of each variable can be visualized using a partial dependence plot, which indicates the effect of the predictor on the response when all other predictors are averaged. We calculated the final relative influence and partial dependence of each variable by averaging these outputs across all 25 model splits. Finally, we predict the probability of paramyxovirus hosting across unsampled or negative bat species, including predictions where both citation count measures were held at their mean. We likewise average these predictions for each species across the 25 model splits and use a 90% sensitivity threshold with the PresenceAbsence package to stratify results into binary predictions of suspect and non-suspect paramyxovirus hosts (Freeman and Moisen 2008). This package can be applied to the raw probability values from our model to examine other potential cutoffs.

- **Developed by:** Maya M. Juman, Daniel J. Becker, Barbara A. Han
- **Model type:** Boosted Regression Trees
- **Language(s) (NLP):** R
- **License:** CC-BY 4.0

## Uses

This model can be adapted to host prediction studies in other systems or used to replicate results from Juman et al. 2025. This is a hypothesis-generating rather than hypothesis-testing approach, and it therefore does not reveal underlying mechanisms of host suitability.

## Bias, Risks, and Limitations

To account for the possibility of the model recapitulating sampling biases, we ran separate models that held citation count measures (a proxy for study effort) at their mean.

## How to Get Started with the Model

Three R scripts are provided here. The first (01_virusdata.R) collects and distills viral data from VIRION (Carlson et al. 2022) and Juman et al. 2025. The second (02_traits.R) harmonizes trait data from public trait databases to the virus data. The third (03_BRTs.R) contains code for hyperparameter tuning and model development.

## Training Details

### Training Data

Virus data were distilled from a publicly available database (<a href="https://doi.org/10.1101/2025.03.10.642350" target="_blank">Juman et al. 2025</a>). Trait data was obtained from public databases (full references available in <a href="https://doi.org/10.1101/2025.09.11.675601" target="_blank">Juman et al. 2025</a>). All training data is located in subfolder "training data" in this repository.

### Training Procedure and Hyperparameters

Hyperparameter tuning and model training code is included in 03_BRTs.R. We train our BRTs on 80% of the complete dataset, using 25 splits, 3000 trees with a Bernoulli error distribution, and five-fold cross-validation to prevent overfitting.

## Evaluation

### Testing Data, Factors & Metrics

#### Testing Data

The model was tested on 20% of the dataset reserved for testing purposes (03_BRTs.R).

#### Metrics

We calculate the corrected area under the curve (AUC) as a performance metric on both our training and test dataset using target shuffling methods. 

### Results

Our BRT analysis identified paramyxovirus PCR-positive pteropodid bat species with up to 92% accuracy (AUC = 0.92, SE = 0.01; corrected AUC = 0.79, SE = 0.02). The top six variables identified by the model as influential predictors of host status were virus-related PubMed citations, total range area, all PubMed citations, adult body length, ear length, and evolutionary distance (equal splits); all other predictors had <3% importance. Partial dependence plots revealed that both citation counts, as well as range area and ear length, were positively associated with paramyxovirus host status, while evolutionary distance was negatively associated with host status. Adult body length was more non-linearly associated with host status, suggesting that pteropodid species of small-to-medium sizes are more likely to be hosts. To account for the fact that species with higher citations are more likely to be predicted as hosts, we ran our predictions when fixing total and virus-related citations at their mean. Here, the three “novel” predicted species were Epomops franqueti, Epomophorus wahlbergi, and Macroglossus minimus, all of which had positivity probabilities above 15%.

## Citation

Juman, M.M., McDonough, M.M., Ferguson, A.W., et al. Museum collections and machine learning guide discovery of novel coronaviruses and paramyxoviruses. 2025. Preprint (bioRxiv): 10.1101/2025.09.11.675601.

## Model Card Contact

mmj38@cam.ac.uk
