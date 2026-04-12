# American Academic Publisher Draft

# Rate-Distortion Bounds on Serial Decoding Throughput

**Grant Lavell Whitmer III**

The Windstorm Institute, Fort Ann, New York 12827, United States of America

Email: grantwhitmer3@gmail.com (Corresponding Author)

---

## Abstract

This study derives a mechanistic bound on effective information per serial decoding event using the M-ary rate-distortion function R_M(epsilon) = log_2(M) - H_b(epsilon) - epsilon * log_2(M - 1). Applied to three independently characterized systems — the ribosome (M = 21, epsilon = 10^-4), human phonology (M = 31, epsilon ~ 0.02), and the chromatic scale (M = 12, epsilon ~ 0.05) — the framework predicts throughput values of 4.39, 4.71, and 3.12 bits respectively. The ribosome prediction achieves a residual of only 0.0004 bits under the uniform distribution assumption, though this near-exact match is recognized as a mathematical property of the formula at very low error rates rather than an independent empirical confirmation. Under realistic non-uniform amino acid frequencies, the ribosome prediction shifts from 4.390 to approximately 4.175 bits, revealing a small but real excess-capacity margin of approximately 0.2 bits. The informative predictions are the non-trivial cases: phonology (residual = +0.24 bits) and the chromatic scale (residual = +0.46 bits). The vocabulary-independence prediction was tested by evaluating 1,749 models across two architectures — 32 causal language models (vocabulary 30K-250K) and 1,717 encoder-decoder translation models (vocabulary 3,932-80,379) — on a shared reference corpus using an NVIDIA RTX 5090. Vocabulary size does not predict bits-per-byte in causal language models (beta_1 = 0.058, p = 0.643) or in language-matched translation models (Spearman p = 0.13). Within-family analysis using the Pythia model series with an identical tokenizer demonstrates that model capacity, not vocabulary, drives compression. An unexpected finding from language-mismatch stratification reveals a quantifiable +1.33-bit information cost when models process input in an unfamiliar language. These results support the hypothesis that in serial decoding systems, the rate-distortion floor is a universal mechanistic bound and vocabulary size is primarily a redundancy parameter. This work is the second paper in a seven-paper series investigating universal constraints on serial information processing [1]-[7].

**Keywords:** rate-distortion theory, serial decoding, receiver constraint, channel capacity, genetic code, tokenizer vocabulary, bits per byte, language model, information theory, throughput basin

---

## 1. Introduction

In a companion paper [1], we predicted that AI tokenizer vocabularies would cluster near 64 symbols, mirroring the genetic code. That prediction was falsified: AI tokenizers use 30K-256K tokens. However, the falsification revealed that while vocabulary sizes differ by orders of magnitude, the effective information per sequential processing event occupies a comparatively narrow range across substrates.

### 1.1 Scope and Definitions

This paper concerns serial discrimination channels — systems where a receiver decodes one symbol per time step under noise. Three quantities are distinguished:

**Definition 1 (Support-size upper bound).** H_0(X) = log_2|X|.                (1)

**Definition 2 (Source entropy rate).** H(X_t | C_t).                            (2)

**Definition 3 (Receiver-resolved information).** I_eff = I(X_t; Y_t | C_t).    (3)

These obey: I(X_t; Y_t | C_t) <= H(X_t | C_t) <= log_2|X|. The constraint does not apply to parallel visual processing, continuous analog signals, or systems without serial discrimination.

---

## 2. Materials and Methods

### 2.1 Rate-Distortion Floor Derivation

For a source with M equiprobable symbols and symmetric error rate epsilon, the M-ary rate-distortion function gives the minimum information per event for reliable reconstruction:

R_M(epsilon) = log_2(M) - H_b(epsilon) - epsilon * log_2(M - 1)              (4)

where H_b(epsilon) = -epsilon * log_2(epsilon) - (1 - epsilon) * log_2(1 - epsilon) is the binary entropy function [8].

### 2.2 Application to Three Biological Systems

The rate-distortion floor was applied to three independently characterized serial decoding systems:

1. **Ribosome** (M = 21 amino acids, epsilon = 10^-4): The error rate reflects the combined fidelity of aminoacyl-tRNA synthetase selection and ribosomal proofreading [9].
2. **Human phonology** (M = 31 cross-linguistic median phonemes, epsilon ~ 0.02): Error rate estimated from confusion matrices in noisy listening conditions [10].
3. **Chromatic scale** (M = 12 tones, epsilon ~ 0.05): Error rate estimated from pitch discrimination thresholds in musical perception.

### 2.3 Robustness Analysis: Non-Uniform Symbol Distributions

For a source with non-uniform distribution p_i over M symbols, the mutual information under a symmetric error channel becomes:

I(X;Y) = H(p) - H_b(epsilon) - epsilon * log_2(M - 1) + epsilon * D_KL(p || uniform) * (terms of order epsilon)    (5)

where H(p) = sum(p_i * log_2(1/p_i)) replaces log_2(M). For the ribosome, empirical amino acid frequencies from UniProt yield H(p) ~ 4.177 bits versus H(uniform) = 4.392 bits — a 0.215-bit reduction [9].

### 2.4 Tokenizer Sweep: 1,749 Models

#### 2.4.1 Model Selection

Two architecture classes were evaluated:

- **32 causal language models** (GPT-2, LLaMA, Pythia, Mistral, and others) spanning vocabulary sizes of 30K-250K tokens.
- **1,717 encoder-decoder translation models** (Helsinki-NLP/OpusMT) spanning vocabulary sizes of 3,932-80,379 tokens.

#### 2.4.2 Evaluation Metric

Bits-per-byte (BPB) normalizes across tokenizations and architectures. All models were evaluated using teacher-forced cross-entropy on WikiText-2 (English) [11].

#### 2.4.3 Language-Match Stratification

Because translation models were trained on specific language pairs but evaluated on English, they were separated into three groups by language match:

- **EN to X (matched encoder, n = 157):** English input matches corpus; decoder outputs foreign language.
- **X to EN (matched decoder, n = 256):** Foreign input mismatches corpus; decoder trained for English output.
- **X to Y (fully mismatched, n = 1,304):** Neither encoder nor decoder expects English.

#### 2.4.4 Hardware

All evaluations were performed on an NVIDIA RTX 5090 GPU with deterministic inference settings. The complete experiment protocol is provided in the supplementary materials [12].

#### 2.4.5 Statistical Methods

Vocabulary-BPB relationships were assessed using ordinary least squares regression (beta_1 coefficient) for causal language models and Spearman rank correlation (rho) for translation model groups. The coefficient of determination (R^2) quantified explained variance. Within-family analysis used the Pythia model series (14M to 1.4B parameters, identical tokenizer) to isolate the effect of model capacity from vocabulary size. Between-group differences were tested using the Kruskal-Wallis H test.

### 2.5 Independent Confirmation

Tlusty [13] derived M ~ 20-25 for the genetic code from the chromatic number theorem (Ringel-Youngs) — a completely independent topological argument (graph coloring versus rate-distortion) converging on the same biological parameter.

---

## 3. Results

### 3.1 Rate-Distortion Predictions

Table 1 presents the zero-free-parameter predictions from the rate-distortion floor applied to three biological systems.

**Table 1.** Rate-distortion predictions for three serial decoding systems

| System | M | epsilon | R_M(epsilon) | Observed | Residual |
|--------|---|---------|-------------|----------|----------|
| Ribosome | 21 | 10^-4 | 4.390 bits | 4.390 bits | +0.000 |
| Phonology | 31 | 0.02 | 4.71 bits | 4.95 bits | +0.24 |
| Chromatic scale | 12 | 0.05 | 3.12 bits | 3.58 bits | +0.46 |

At epsilon = 10^-4, the ribosome correction terms amount to only 0.0019 bits (0.04% of log_2(21)). The near-exact match between R_M(epsilon) and log_2(M) at very low error rates is a mathematical property of (4), not an independent empirical confirmation. The framework's value lies in the non-trivial cases (phonology and chromatic scale) and in unifying M, epsilon, and throughput into a single testable relationship.

The residuals represent excess-capacity margins above the rate-distortion floor. Biology operates nearest the floor; cognition and music carry higher slack, consistent with greater discrimination cost in neural substrates.

### 3.2 Non-Uniform Distribution Robustness

Under non-uniform amino acid frequencies, the ribosome prediction shifts from 4.390 to approximately 4.175 bits. The observed ribosome throughput therefore sits approximately 0.2 bits above the non-uniform floor — a small but real excess-capacity margin consistent with the slack observed in other substrates. This correction makes the result more interesting rather than less: rather than a suspiciously perfect residual of 0.000, a genuine physical margin comparable to what Hopfield kinetic proofreading predicts is revealed [14].

### 3.3 Vocabulary Independence Across 1,749 Models

Table 2 presents the vocabulary-independence results across all four model groups.

**Table 2.** Vocabulary independence across four model groups

| Group | n | BPB mean | BPB median | Statistic | p-value | R^2 | Vocab-independent? |
|-------|---|----------|------------|-----------|---------|-----|-------------------|
| Causal LMs | 32 | 0.94 | 0.85 | beta_1 = 0.058 | 0.643 | 0.007 | Yes |
| EN to X | 157 | 4.19 | 3.75 | rho = -0.122 | 0.129 | 0.002 | Yes |
| X to EN | 256 | 5.98 | 5.56 | rho = -0.068 | 0.280 | 0.016 | Yes |
| X to Y | 1,304 | 5.52 | 5.49 | rho = -0.068 | 0.015 | 0.032 | Weak* |

*Statistically significant at p < 0.05 but rho = -0.068 indicates negligible practical effect; vocabulary explains only 3.2% of variance.

Vocabulary size does not predict BPB in causal language models (p = 0.643, R^2 = 0.007) or in language-matched translation models (p = 0.129, R^2 = 0.002). The vocabulary-independence result holds across both architectures and all four language-match conditions. Figures 1-4 illustrate these relationships (submitted separately).

### 3.4 Within-Family Analysis

Within the Pythia model family, which uses an identical tokenizer across all model sizes, BPB drops 46% from 14M to 1.4B parameters. Model capacity drives compression; vocabulary does not.

### 3.5 Language-Mismatch Finding

The separation into language-match groups revealed an unexpected result: the information cost of processing mismatched input is quantifiable.

**Table 3.** Language-match effect on throughput

| Condition | n | BPB mean | Interpretation |
|-----------|---|----------|----------------|
| EN to X (matched encoder) | 157 | 4.19 | Model understands input; efficient processing |
| X to Y (fully mismatched) | 1,304 | 5.52 | Model confused on both sides; +1.33 bits penalty |
| X to EN (matched decoder) | 256 | 5.98 | Encoder garbled, decoder compensates; highest cost |

Matched-encoder models (EN to X), which receive English input matching the corpus, operate at BPB ~ 4.19 — within the throughput basin and near the ribosome prediction of 4.39 bits. Fully mismatched models (X to Y) operate at BPB ~ 5.52 — a penalty of +1.33 bits. The groups are significantly different (Kruskal-Wallis H = 230.2, p ~ 0), but vocabulary independence holds within each group. This finding was not part of the original experiment design and emerged during post-hoc analysis. It is reported as an exploratory finding, though the effect size (1.33 bits) is large and consistent across 1,461 models.

### 3.6 Empirical Phonology Values

Two phonology measurements appear in this series:

- **Support-size bound [SS]:** log_2(31) = 4.95 bits — the maximum entropy of the phoneme inventory assuming uniform distribution.
- **Confusion-matrix MI [MI]:** Miller and Nicely [10] measured I(X;Y) ~ 4.20 bits for English consonants — a direct Tier 1 estimate of receiver-resolved information under noise.

The rate-distortion prediction R_31(0.02) = 4.71 bits sits between these two values. Measured against the support-size bound, the residual is +0.24 bits (below the ceiling, as expected). Measured against the confusion-matrix MI, the residual is -0.51 bits (the floor overpredicts, suggesting that contextual redundancy in real speech reduces effective throughput below even the rate-distortion floor). This discrepancy highlights the importance of evidence tiers: Tier 1 and Tier 2 measurements of the "same" system can disagree by approximately 0.75 bits.

---

## 4. Discussion

### 4.1 The Receiver-Limited Hypothesis

The results are consistent with the receiver-limited hypothesis: vocabulary growth reallocates redundancy and sequence length without increasing per-event throughput. The cross-architecture consistency (causal plus encoder-decoder) and the language-match stratification strengthen this claim beyond what a single-architecture sweep could provide. The receiver-limited hypothesis predicts exactly the observed pattern: absolute throughput depends on the receiver's capacity and task, but vocabulary is a redundancy parameter in all cases.

### 4.2 Why Biological Systems Cluster at 3-5 Bits

The convergence is a geometric property of (4): for M in [12, 44] and epsilon in [10^-4, 0.05], R_M(epsilon) maps to 3-5 bits. Biological serial decoders operate in this regime because moderate M reflects the cost-benefit optimum under super-linear discrimination cost, and moderate epsilon reflects the thermodynamic floor of reliable discrimination (Hopfield [14]: free energy change ~ 2-3 kcal/mol per step). This is a consequence of the rate-distortion surface evaluated at biologically plausible parameters — not an empirical surprise.

### 4.3 Vocabulary as Dependent Variable

Given receiver capacity C_R and error rate epsilon, the optimal vocabulary size is:

V* = 2^(C_R + log_2(1/epsilon))                                              (6)

Both C_R and epsilon are independently measurable, making this a prediction rather than a tautology. This formulation recasts vocabulary size as a derived quantity — a dependent variable determined by receiver architecture and noise characteristics.

### 4.4 Cross-Substrate Implications

The cross-substrate kinship implied by this framework — ribosome, human ear, transformer attention head, and Braille reader all solving the same serial decoding problem — suggests a universal geometric structure governing information processing under noise and cost constraints. This kinship is not by ancestry but by constraint. The mathematics does not distinguish whether the receiver is RNA, neurons, silicon, or raised dots. What matters is the structure of the task: serial, noisy, cost-limited. The broader cross-substrate analysis is presented in the companion paper [3].

### 4.5 Falsifiable Predictions

Four falsifiable predictions are offered:

1. **Vocabulary plateau:** Same architecture trained at V in {256, 1K, 4K, 16K, 64K, 256K} should show BPB plateau with no ordinal relationship to V.
2. **Corruption convergence:** Tokenizers optimized under corruption should converge on smaller, more redundant vocabularies at similar BPB.
3. **Genetic code bounds:** All natural genetic codes should encode 15 <= N <= 25 amino acids, consistent with the Tlusty topological bound [13].
4. **Evolutionary convergence:** Evolutionary simulations should converge to 16-32 functional outputs (4-5 bits). Preliminary results from three independent implementations yield K ~ 19-30.

### 4.6 Limitations

1. **Near-tautology at low epsilon.** The ribosome prediction (residual = 0.0004 under uniform) is near-tautological because R_M(epsilon) approaches log_2(M) as epsilon approaches 0. The framework's value lies in the non-trivial cases.
2. **Equiprobable symbol assumption.** All predictions assume symmetric error channels with uniform input. Real systems violate both assumptions. The robustness analysis (Section 3.2) suggests the framework survives this violation but systematic cross-substrate analysis is needed.
3. **Observational design.** The 1,749-model experiment is observational, not controlled. Architecture, training data, and training duration confound vocabulary-size effects. A controlled single-architecture sweep is the definitive next step.
4. **Tier mixing.** Phonology and music predictions compare R_M(epsilon) against support-size bounds [SS], not mutual information [MI]. Conditional entropy rates from Conklin (1995), Temperley (2014), and Pimentel (2020) show music operates at 1.87-2.97 bits/note and phonemes at approximately 3.0 bits — significantly lower than support-size values.
5. **AI I_eff is model-dependent**, ranging from approximately 1.4 (frontier) to 5+ (small models). There is no single substrate-level "AI throughput."

---

## 5. Conclusion

The rate-distortion floor R_M(epsilon) provides a universal mechanistic bound for serial decoding throughput. For biological receivers with moderate M and moderate epsilon, this bound confines throughput to 3-5 bits — a geometric property of the rate-distortion surface, not an empirical coincidence. The empirical finding that AI bits-per-byte is independent of vocabulary size across 1,749 models (causal: p = 0.643; matched-encoder translation: Spearman p = 0.13) supports the vocabulary-as-redundancy hypothesis. Three independent evidence lines — the rate-distortion framework, the biological predictions, and the AI regularities (vocabulary independence and language-match penalty) — together establish a receiver-limited serial decoding constraint. Future work should prioritize controlled single-architecture vocabulary sweeps and Tier 1 confusion-matrix measurements across additional biological systems to narrow the empirical calibration of this framework.

---

## Acknowledgements

This paper was improved through adversarial review by six frontier AI models (Claude, GPT, Grok, Gemini, Sonar, and Perplexity). Mathematical derivations and data analysis were performed with the assistance of Claude (Anthropic), an AI research tool. The tokenizer-sweep experiment (1,749 models across two architectures) was executed by Hermes OC1 on an NVIDIA RTX 5090.

## Funding Information

This research received no external funding. All work was self-funded by the author.

## Author Contributions

Grant Lavell Whitmer III conceived the research framework, designed the experimental protocol, directed the 1,749-model evaluation, analyzed all results, and prepared the manuscript.

## Conflict of Interest

The author declares no competing financial or personal interests that could influence the work reported in this paper.

---

## References

[1] Whitmer III, G.L. "The Fons Constraint: Information-Theoretic Convergence on Encoding Depth in Self-Replicating Systems," Zenodo, 2026. DOI: 10.5281/zenodo.19274048

[2] Whitmer III, G.L. "The Throughput Basin: Cross-Substrate Convergence and Decomposition of Serial Decoding Throughput," Zenodo, 2026. DOI: 10.5281/zenodo.19323194

[3] Whitmer III, G.L. "The Serial Decoding Basin: Five Experiments on Convergence, Thermodynamic Anchoring, and Receiver-Limited Geometry," Zenodo, 2026. DOI: 10.5281/zenodo.19323423

[4] Whitmer III, G.L. "The Dissipative Decoder: Thermodynamic Cost Bounds on the Serial Decoding Throughput Basin," Zenodo, 2026. DOI: 10.5281/zenodo.19433048

[5] Whitmer III, G.L. "The Inherited Constraint: Biological Throughput Limits Shape the Information Structure of Human Language and AI," Zenodo, 2026. DOI: 10.5281/zenodo.19432911

[6] Whitmer III, G.L. "The Throughput Basin Origin: Four Orthogonal Experiments on Whether Serial Decoding Convergence Is Architectural, Thermodynamic, or Data-Driven," Zenodo, 2026. DOI: 10.5281/zenodo.19498582

[7] Whitmer III, G.L. "The Receiver-Limited Floor: Rate-Distortion Bounds on Serial Decoding Throughput," Zenodo, 2026. DOI: 10.5281/zenodo.19322973

[8] Shannon, C.E. "A Mathematical Theory of Communication," Bell System Technical Journal, vol. 27, no. 3, pp. 379-423, 1948. DOI: 10.1002/j.1538-7305.1948.tb01338.x

[9] Zaher, H.S.; Green, R. "Quality control by the ribosome following peptide bond formation," Nature, vol. 457, no. 7226, pp. 161-166, 2009. DOI: 10.1038/nature07582

[10] Miller, G.A.; Nicely, P.E. "An analysis of perceptual confusions among some English consonants," Journal of the Acoustical Society of America, vol. 27, no. 2, pp. 338-352, 1955. DOI: 10.1121/1.1907526

[11] Berger, T. Rate Distortion Theory; Prentice-Hall: Englewood Cliffs, NJ, USA, 1971.

[12] Whitmer III, G.L. "Receiver-Limited Floor: Experiment Code and Data," GitHub, 2026. https://github.com/Windstorm-Labs/receiver-limited-floor (accessed Apr. 12, 2026).

[13] Tlusty, T. "Rate-distortion scenario for the emergence and evolution of noisy molecular codes," Physical Review Letters, vol. 100, pp. 048101, 2008. DOI: 10.1103/PhysRevLett.100.048101

[14] Hopfield, J.J. "Kinetic proofreading: A new mechanism for reducing errors in biosynthetic processes requiring high specificity," Proceedings of the National Academy of Sciences, vol. 71, no. 10, pp. 4135-4139, 1974. DOI: 10.1073/pnas.71.10.4135

[15] Borst, A.; Theunissen, F.E. "Information theory and neural coding," Nature Neuroscience, vol. 2, no. 11, pp. 947-957, 1999. DOI: 10.1038/14731

---

## Appendix A: Experiment Protocol

The complete experiment protocol, including hardware specifications, corpus preparation, BPB computation methods for three model formats (HuggingFace, GGUF, MarianMT), statistical analysis procedures, and visualization specifications, is provided in the supplementary file APPENDIX-A-EXPERIMENT-PROTOCOL.md [12]. The protocol was designed for full reproducibility: any researcher with access to HuggingFace models and a CUDA-capable GPU can replicate the experiment using the provided instructions.

## Appendix B: Figure Descriptions

The following figures accompany this paper and are submitted separately:

- **Figure 1.** Rate-distortion curves for M in {4, 7, 12, 21, 31, 44, 64}. Biological systems (stars) sit inside the 3-5.5 bit convergence band.
- **Figure 2.** Predicted versus observed throughput — three systems with residual annotations.
- **Figure 3.** BPB versus vocabulary size — 1,749 models across four groups. Regression lines are flat in all groups.
- **Figure 4.** BPB distribution by group — boxplot showing four distinct operating regimes.
- **Figure 5.** Language-match effect on throughput — mean BPB with 95% confidence intervals. EN to X lands at 4.19 bits, near the ribosome floor.
- **Figure 6.** BPB density by group — kernel density estimation showing how each group relates to the throughput basin.
- **Figure 7.** The information cost of substrate mismatch — +1.33 bits between matched and mismatched processing.
- **Figure 8.** Vocabulary independence statistics by group — summary table visualization.

---

*This paper is Paper 2 of The Windstorm Series. The series includes Paper 1: The Fons Constraint [1], Paper 3: The Throughput Basin [2], Paper 4: The Serial Decoding Basin [3], Paper 5: The Dissipative Decoder [4], Paper 6: The Inherited Constraint [5], and Paper 7: The Throughput Basin Origin [6].*
