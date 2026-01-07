# Methodology

This project adopts a **Bayesian time-series framework** to study **time-varying risk transmission and tail dependence** in multivariate financial systems.  
Rather than relying on a single static specification, the methodology is designed to progressively relax modeling assumptions, moving from average dynamics to time-varying and state-dependent risk propagation.

The overall strategy bridges rigorous econometric modeling with quantitative risk analysis, ensuring both interpretability and practical relevance.

---

## 1. Baseline VAR Framework

The analysis begins with a standard Vector Autoregression (VAR) model to characterize linear interdependencies among multiple financial variables.

Under the baseline specification:
- Model parameters are assumed to be time-invariant.
- Innovations follow a multivariate Gaussian distribution.
- Dynamic interactions are summarized using impulse response functions and forecast error variance decompositions.

This framework provides a benchmark for understanding **average shock transmission mechanisms** before introducing additional layers of complexity.

---

## 2. Bayesian Estimation Strategy

All VAR-based models are estimated within a **Bayesian framework**, which offers several advantages:
- Coherent probabilistic inference.
- Explicit incorporation of prior information.
- Full posterior distributions for parameters and latent states.

Conjugate prior structures are employed for regression coefficients and covariance matrices to ensure computational tractability.  
Posterior inference is conducted using **Markov Chain Monte Carlo (MCMC)** sampling, with appropriate burn-in periods discarded to mitigate initialization effects.

---

## 3. Time-Varying Parameter VAR (TVP-VAR)

To capture structural change and evolving market dynamics, the model is extended to a **Time-Varying Parameter VAR (TVP-VAR)**.

In this framework:
- Regression coefficients evolve according to stochastic state equations.
- Parameter dynamics follow random walk processes.
- The model accommodates gradual regime shifts and evolving dependence structures.

The TVP-VAR allows the analysis of **how shock propagation mechanisms change over time**, reflecting shifts in macroeconomic conditions, market sentiment, and risk perception.

---

## 4. Conditional VAR and Tail Risk Analysis

Beyond time variation, the methodology emphasizes **tail risk and state dependence** through conditional VAR analysis.

Conditional specifications focus on periods associated with:
- Elevated volatility.
- Extreme negative returns.
- Stressed market environments.

This approach isolates **risk transmission under adverse conditions**, enabling comparison between normal and tail states and highlighting potential amplification effects that are not captured by unconditional models.

---

## 5. Dynamic Risk Transmission and Interpretation

The constant-parameter VAR and TVP-VAR frameworks are compared along multiple dimensions:
- Stability of estimated relationships.
- Sensitivity to regime changes.
- Evolution of impulse response functions over time.

Rather than emphasizing model fit alone, the analysis focuses on the **economic interpretation of evolving shock transmission patterns**, particularly under tail events.

---

## 6. Relevance for Quantitative Risk Modeling

By integrating Bayesian VAR, TVP-VAR, and conditional analysis, the methodology provides a unified framework to study:
- Time-varying dependence structures.
- State-dependent risk propagation.
- Tail-driven amplification of shocks.

The framework is directly relevant for quantitative applications such as portfolio risk management, stress testing, and forecasting under uncertainty.
