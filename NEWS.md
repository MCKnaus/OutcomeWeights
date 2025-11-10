# OutcomeWeights 0.2.0

* Average treatment effect on the treated and on the untreated are now supported by dml_with_smoother

* get_outcome_weights() now compatible with all average_treatment_effect() options for causal_forest and instrumental_forest


# OutcomeWeights 0.1.1

* Fixed bug throwing an error if dml_with_smoother is called without Z, regardless whether it is a required input or not.

* Test added to rule out similar errors in the future


# OutcomeWeights 0.1.0

* Initial CRAN submission.
