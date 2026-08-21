# Lista poprawek do pracy magisterskiej

**Dokument źródłowy:** `master_thesis (5).pdf`  
**Zakres przeglądu:** cała praca z wyjątkiem abstraktu  
**Sposób użycia:** realizuj zadania od sekcji **P0** do **P3**. Każde zadanie zawiera lokalizację, opis problemu, proponowaną zmianę i kryterium ukończenia.

---

## Legenda priorytetów

- **P0 — krytyczne:** błędy, które mogą podważać poprawność lub kompletność pracy.
- **P1 — wysokie:** niespójności metodologiczne, statystyczne lub interpretacyjne.
- **P2 — średnie:** poprawki merytoryczne i językowe zwiększające wiarygodność.
- **P3 — redakcyjne:** skład, terminologia, bibliografia i estetyka.

---

# P0 — poprawki krytyczne

## [ ] P0-01. Uzupełnić brakujące cytowanie oznaczone jako `[?]`

**Lokalizacja:** Chapter 1 → `Genetic Architecture`, fragment dotyczący `IL23R` i wariantu R381Q.

**Problem:**  
W tekście znajduje się nierozwiązane odwołanie bibliograficzne:

```text
reduced secretion of pro-inflammatory cytokines [?]
```

**Co zrobić:**

1. Zidentyfikować publikację źródłową dla funkcjonalnego efektu wariantu R381Q.
2. Dodać wpis do pliku `.bib`.
3. Zastąpić `[?]` prawidłowym `\cite{...}`.
4. Ponownie skompilować dokument i sprawdzić log LaTeX-a pod kątem `undefined citations`.

**Kryterium ukończenia:**  
W całym PDF nie występuje ciąg `[?]`, a bibliografia zawiera pełną pozycję.

---

## [ ] P0-02. Uzupełnić dokładne wersje oprogramowania

**Lokalizacja:** Methods → `Computational Environment and Software`, tabela oprogramowania.

**Obecnie:**

```text
GATK — 4
SAMtools — Not specified
BCFtools — Not specified
g:Profiler — [VERSION]
```

**Problem:**  
Deklarowana reprodukowalność nie jest zgodna z niepełnymi wersjami narzędzi.

**Co zrobić:**

- odzyskać wersje z logów, `conda env export`, `pip freeze`, historii terminala, plików środowiska lub metadanych workflow;
- wpisać pełne numery wersji;
- dodać cytowanie g:Profiler.

**Jeśli wersji nie da się odzyskać:**

Zamiast:

```text
Not specified
```

użyć:

```text
Exact patch version not retained
```

**Kryterium ukończenia:**  
W tabeli nie występują placeholdery `[VERSION]` ani `Not specified`.

---

## [ ] P0-03. Naprawić niewidoczny diagram metodologii

**Lokalizacja:** Methods, początek rozdziału, Figure 4.1 / okolice strony 31.

**Problem:**  
W PDF widoczny jest podpis lub odwołanie do ryciny, ale sam diagram nie renderuje się poprawnie.

**Najbezpieczniejsze rozwiązanie:**

1. Wyeksportować diagram jako wektorowy PDF.
2. Zastąpić:

```latex
\includesvg[width=\textwidth]{...}
```

przez:

```latex
\includegraphics[width=\textwidth]{images/methods_workflow.pdf}
```

3. Sprawdzić, czy plik jest rzeczywiście dodany do projektu.
4. Skompilować dokument od zera.

**Alternatywa dla SVG:**  
Kompilować z `--shell-escape` i upewnić się, że pakiet `svg` działa poprawnie.

**Kryterium ukończenia:**  
Rycina jest widoczna, czytelna przy szerokości A4 i ma prawidłowy numer oraz podpis.

---

## [ ] P0-04. Ponownie wygenerować macierze pomyłek z `Predicted_label`

**Lokalizacja:** Results → `Comparative Model Performance`, Figure 5.3.

**Problem:**  
Tabela metryk korzysta z:

```python
pipeline.predict()
```

a obecna macierz pomyłek została zbudowana z:

```python
Probability_UC >= 0.5
```

Dla SVM są to różne reguły klasyfikacji.

**Co zrobić:**

1. Użyć zapisanej kolumny `Predicted_label`.
2. Połączyć predykcje out-of-fold z wszystkich 50 zewnętrznych foldów.
3. Zbudować row-normalized confusion matrix.
4. Upewnić się, że SVM daje wartości zbliżone do:
   - sensitivity ≈ `0.612`,
   - specificity ≈ `0.407`.
5. Zaktualizować plik ryciny i opis.

**Poprawny podpis:**

```latex
\caption{Row-normalized confusion matrices calculated from pooled
out-of-fold class predictions across all 50 outer test folds. Each
individual contributed one prediction per repetition and therefore
appeared 10 times. Class labels were obtained from the native
\texttt{predict()} output of each classifier.}
```

**Kryterium ukończenia:**  
Wartości w macierzy pomyłek są zgodne z tabelą metryk.

---

## [ ] P0-05. Ujednolicić wartości na wykresie ROC/PR z metodą agregacji

**Lokalizacja:** Results → Figure 5.2 i tekst pod ryciną.

**Problem:**  
Krzywe są rysowane z połączonych 1440 predykcji out-of-fold, ale wartości w legendzie odpowiadają średnim z 10 repetition-level estimates.

**Rekomendowane rozwiązanie:**

- pozostawić pooled curves,
- policzyć AUROC i AP bezpośrednio z dokładnie tych samych 1440 prawdopodobieństw,
- takie wartości umieścić w legendzie,
- w tabeli zachować średnia ± SD z 10 powtórzeń.

**Dodać wyjaśnienie:**

```latex
The values reported in Table~\ref{tab:model_performance} summarize
repetition-level performance, whereas Figure~\ref{fig:model_roc_pr}
shows curves calculated from pooled out-of-fold probabilities across
all 50 outer test folds.
```

**Kryterium ukończenia:**  
Legenda opisuje wartości policzone z danych użytych do narysowania krzywej.

---

## [ ] P0-06. Przebudować stronę skrótów

**Lokalizacja:** strona `Shortcuts`.

**Błędy:**

- `Artifical Intelligence` → literówka,
- brak nawiasu w `GWAS`,
- niespójne liczby pojedyncze i mnogie,
- brak kluczowych skrótów,
- obecność skrótów nieużywanych.

**Zmienić tytuł:**

```text
Shortcuts
```

na:

```text
List of Abbreviations
```

**Usunąć, jeśli nieużywane:** LLM, NLP, API, EHR.

**Dodać co najmniej:** UC, IBD, WES, NGS, BAM, VCF, GATK, VQSR, VEP, CADD, AUROC, TPE, SHAP, SVM, k-NN, SNV, CNV, indel.

**Poprawić:**

```text
AI (Artifical Intelligence)
```

na:

```text
AI (Artificial Intelligence)
```

**Kryterium ukończenia:**  
Lista zawiera wszystkie kluczowe skróty użyte w pracy i nie zawiera literówek.

---

# P1 — spójność metodologiczna i statystyczna

## [ ] P1-01. Usunąć niepotwierdzone określenie `unrelated individuals`

**Lokalizacja:** Results → `Study Cohort and Final Modelling Matrix`.

**Obecnie:**

```text
The final analysis included 144 unrelated individuals
```

**Zamienić na:**

```text
The final analysis included 144 individuals
```

**Kryterium ukończenia:**  
Nigdzie w pracy nie jest deklarowany brak pokrewieństwa bez analizy kinship.

---

## [ ] P1-02. Usunąć `accuracy` z listy raportowanych metryk

**Lokalizacja:** Results → początek `Comparative Model Performance`.

**Obecnie:**

```text
balanced accuracy, sensitivity, specificity, precision, F1-score,
and accuracy
```

**Zamienić na:**

```text
balanced accuracy, sensitivity, specificity, precision, and F1-score
```

---

## [ ] P1-03. Poprawić podpis tabeli wydajności modeli

**Lokalizacja:** Results → Table 5.2.

**Zamienić na:**

```latex
Values are presented as mean $\pm$ standard deviation across
10 repetitions of nested five-fold cross-validation. Within each
repetition, metric values were first averaged across the five outer
test folds.
```

---

## [ ] P1-04. Poprawić podpis boxplotu AUROC

**Lokalizacja:** Results → Figure 5.1.

```latex
\caption{Distribution of AUROC values across 10 repetitions of nested
five-fold cross-validation. Each point represents one repetition-level
estimate obtained by averaging AUROC across the five corresponding
outer test folds. Boxes show the median and interquartile range. The
dashed line indicates chance-level discrimination.}
```

---

## [ ] P1-05. Zmienić `planned comparisons` na `pairwise` lub `post-hoc`

**Lokalizacja:** Results → sekcja statystyczna i podpis Table 5.3.

**Zamienić:**

```text
five planned comparisons
```

na:

```text
five post-hoc pairwise comparisons
```

lub neutralnie:

```text
five pairwise comparisons
```

---

## [ ] P1-06. Poprawić opis identycznych wartości Wilcoxona

**Zamienić:**

```text
the same Wilcoxon rank configuration
```

na:

```text
the same sign pattern of paired differences
```

**Zalecany fragment:**

```latex
In every comparison, Random Forest achieved a higher AUROC in each of
the 10 matched repetitions; all paired differences were positive and
none were equal to zero. This identical sign pattern produced the same
minimum exact two-sided Wilcoxon \(p\)-value for every comparison.
```

---

## [ ] P1-07. Usunąć albo zaraportować `Feature-selection Stability`

**Lokalizacja:** Methods → `Model Interpretation` → `Feature-selection Stability`.

**Rekomendacja:**  
Usunąć podsekcję z Methods i wszystkie zdania o selection frequency, jeśli analiza nie jest raportowana.

**Alternatywa:**  
Dodać do Results tabelę z median feature-selection frequency, top-10 frequency i top-20 frequency.

---

## [ ] P1-08. Zmienić hierarchię sekcji Pathway Analysis i Summary

**Obecny układ:**

```latex
\section{Summary of Main Findings}
\subsection{Functional Pathway Analysis}
```

**Poprawny układ:**

```latex
\section{Functional Pathway Analysis}
...
\section{Summary of Main Findings}
```

**Dodatkowo:**  
Nie umieszczać pustej ryciny informującej o braku istotnych terminów. Wystarczy tekst.

---

## [ ] P1-09. Ogólnie opisać wybór modeli do interpretowalności w Methods

**Zamienić nazwy modeli na regułę:**

```latex
Following completion of the comparative performance evaluation, the two
classifiers with the highest mean repetition-level AUROC were selected
for interpretability analysis. The highest-ranked classifier was treated
as the primary model, while the second-ranked classifier was included as
a comparison model.
```

**W Results dopiero ujawnić:**

```text
Random Forest and Naive Bayes were the two highest-ranked classifiers...
```

---

## [ ] P1-10. Doprecyzować permutation importance

**Lokalizacja:** Methods → `Held-out Permutation Importance`.

**Dodać:**

```latex
Each feature was permuted [N] times within every outer test fold.
```

**Dodać opis agregacji:**

```latex
For each gene, permutation importance was first averaged across repeated
permutations within an outer test fold and was then summarised across
all 50 outer test folds. The proportion of folds with positive
importance was additionally reported as a measure of fold-level
consistency.
```

**Doprecyzować niewybrane cechy:**

```latex
Predictors not retained by the fold-specific selector had an importance
of zero by construction because they did not influence the fitted
classifier.
```

---

## [ ] P1-11. Doprecyzować procedurę SHAP

**Minimalne informacje:**

1. Wyjaśniane było wyjście dla klasy UC.
2. Background dla model-agnostic SHAP pochodził wyłącznie z outer training fold.
3. Niewybrane cechy otrzymywały SHAP = 0.
4. Wyjaśnienia mapowano do pełnych 215 genów.
5. Pooled set zawierał 1440 wyjaśnień per model.
6. Podać limit obliczeń lub liczbę evaluations, jeśli używana.

**Zalecany tekst:**

```latex
SHAP values were calculated separately for each fitted outer-fold model
using observations from the corresponding held-out test fold. A
model-specific explainer was used when supported by the classifier
architecture; otherwise, a model-agnostic permutation explainer was
applied. Reference data for the model-agnostic explainer were drawn
exclusively from the corresponding outer training partition.

Explanations were calculated for the output associated with the positive
UC class. Fold-specific attributions were mapped back to the complete
set of 215 predictors, and features not retained in a given fold were
assigned an attribution of zero. Explanations from all 50 outer test
folds were pooled, yielding 1,440 held-out sample-level explanations per
classifier.
```

---

## [ ] P1-12. Zdefiniować callability operacyjnie

**Lokalizacja:** Methods → `Callability Audit of High-importance Predictors`.

**Należy dopisać dokładnie:**

- co oznacza callable gene,
- jak liczono callability per sample,
- czy używano DP/GQ,
- jak traktowano brak rekordu w variant-only VCF,
- jak traktowano assumed-reference,
- czy callability była liczona na pozycji, wariancie czy genie.

**Szablon:**

```latex
For each gene, callability was defined as [EXACT IMPLEMENTATION
DEFINITION]. Group-specific callability was calculated as the proportion
of samples satisfying this criterion. The signed UC-minus-JPT difference
and its absolute value were then reported.
```

---

## [ ] P1-13. Uogólnić pathway analysis w Methods

**Zamienić:**

```text
for the Random Forest explanations
```

na:

```text
for the highest-ranked classifier
```

**Dodać:** organism `hsapiens`, sources `GO:BP` i `REAC`, custom background 215 predictors, exact FDR method, threshold `adjusted p < 0.05`, minimum 2 contributing genes, g:Profiler version.

**Zalecany fragment:**

```latex
Over-representation analysis was performed for the highest-ranked
classifier using g:Profiler for Homo sapiens. Gene Ontology Biological
Process and Reactome terms were tested using the 215 modelled genes as a
custom background. Results were corrected for multiple testing using
false-discovery-rate adjustment. Terms with adjusted \(p<0.05\) and at
least two contributing genes were considered significant.
```

---

## [ ] P1-14. Wzmocnić ograniczenie dotyczące braku joint genotyping

**Lokalizacja:** Discussion → `Strengths and Limitations`.

**Dodać:**

```latex
The cohort-level VCF was created by merging independently called
variant-only VCF files rather than by joint genotyping of reference-
confidence GVCFs. Consequently, absence of a record could not be
interpreted as a jointly confirmed homozygous-reference genotype. The
assumed-reference approximation used during feature construction may
therefore have interacted with cohort-specific differences in coverage
and variant emission. This limitation is directly related to the
callability differences observed for several highly ranked predictors.
```

---

# P2 — poprawki merytoryczne i interpretacyjne

## [ ] P2-01. Ograniczyć twierdzenia o `prediction of UC`

**Lokalizacja:** tytuł, Introduction, Results, Conclusions.

**Preferowane sformułowanie:**

```text
internal discrimination between UC cases and JPT population-reference samples
```

zamiast:

```text
prediction of ulcerative colitis
```

**Rozważyć ostrożniejszy tytuł:**

```text
Comparative analysis of supervised learning algorithms for internal
discrimination of ulcerative colitis cases and population-reference
samples using aggregated rare genetic variants
```

---

## [ ] P2-02. Zastąpić `pathogenicity scores`

**Globalnie zamienić:**

```text
gene-level pathogenicity scores
```

na:

```text
GenePy-based gene-level burden scores
```

lub:

```text
gene-level GenePy scores
```

---

## [ ] P2-03. Osłabić twierdzenia o SHAP

**Nie używać:**

```text
SHAP ensures clinical plausibility
SHAP ensures biological plausibility
SHAP validates biological findings
```

**Zamiast tego:**

```text
SHAP was used to describe model-level feature contributions and to
support interpretation of the highest-performing classifiers.
```

---

## [ ] P2-04. Poprawić końcowy akapit Introduction

```latex
The remainder of this thesis is organised as follows. Chapters 1--3
provide the theoretical background, covering the clinical and genetic
context of UC, genomic data processing, and the machine-learning methods
used in the study. Chapter 4 describes the analytical workflow and
evaluation protocol. Chapter 5 presents the results, Chapter 6 discusses
their interpretation and limitations, and Chapter 7 summarises the main
conclusions.
```

---

## [ ] P2-05. Skrócić i uakademić fragment o personalised medicine

**Usunąć lub zastąpić fragment zaczynający się od:**

```text
Traditional healthcare usually follows a model that is based on
what works for most people...
```

**Proponowana wersja:**

```latex
Personalised medicine aims to integrate biological and clinical
heterogeneity into prevention, diagnosis, and treatment decisions.
In complex diseases such as UC, germline variation represents only one
component of this heterogeneity and should be considered together with
environmental, microbial, and clinical factors. The present study
therefore evaluates whether exome-derived gene-level features contain
internally reproducible classification signal without treating genomic
information as a complete representation of disease status.
```

---

## [ ] P2-06. Usunąć lub ograniczyć akapit o modelu Lalonde

**Rekomendacja:**  
Usunąć cały akapit albo zostawić jedno zdanie bez procentów:

```latex
This perspective is consistent with multifactorial models of health in
which biological, environmental, behavioural, and healthcare-related
factors jointly influence outcomes.
```

---

## [ ] P2-07. Poprawić zbyt mocne twierdzenia w Clinical Context

**Do przeglądu:**

- `UC is generally more prevalent than Crohn's disease`,
- genomika umożliwia `early prevention`,
- genomika umożliwia spersonalizowaną opiekę podczas flare,
- `the disease predominantly affects workers, specifically males aged 30–44`.

**Co zrobić:**

- dodać źródło lub usunąć generalizacje,
- ograniczyć polskie dane do populacji i okresu źródłowego,
- rozróżnić predisposition od diagnosis/monitoring.

---

## [ ] P2-08. Poprawić nazwę sekcji o regresji logistycznej

**Zamienić:**

```text
LASSO-Regularized Logistic Regression
```

na:

```text
Regularized Logistic Regression
```

---

## [ ] P2-09. Skrócić ogólny przegląd ML

**Do skrócenia:**

- Unsupervised Learning,
- Reinforcement Learning,
- szerokie przykłady AI niezwiązane z badaniem,
- przykłady typu AlphaFold, jeśli nie wspierają celu pracy.

**Cel:**  
Przenieść ciężar na supervised learning, wybrane algorytmy, nested CV, metryki, Wilcoxon i interpretowalność.

---

## [ ] P2-10. Poprawić opis AUROC

**Nie używać:**

```text
AUROC is invariant to class distribution
AUROC is the gold standard
```

**Zamienić na:**

```latex
AUROC provides a threshold-independent summary of ranking performance.
Because it does not describe performance at a single operating point,
it was interpreted together with class-specific and precision--recall
measures.
```

---

## [ ] P2-11. Poprawić interpretację testu Wilcoxona

```latex
The Wilcoxon signed-rank test is a non-parametric paired test based on
the ranks and signs of within-pair differences. Its \(p\)-value
quantifies the compatibility of the observed differences with the null
hypothesis under the assumptions of the test; it is not the probability
that the result occurred by chance.
```

---

## [ ] P2-12. Przenieść source confounding wyżej w Discussion

**Nowa zalecana kolejność ograniczeń:**

1. pełne skonfundowanie class label ze source study,
2. brak external validation,
3. brak joint genotyping i assumed-reference,
4. mały zbiór i repeated reuse,
5. restricted gene panel,
6. brak kalibracji i threshold optimisation.

---

## [ ] P2-13. Zastąpić `harmonisation` bardziej precyzyjnym opisem

**Zamienić:**

```text
Although preprocessing and harmonisation were designed to reduce...
```

na:

```text
Although common downstream processing was applied where possible...
```

---

## [ ] P2-14. Zaktualizować Future Work po wykonanej pathway analysis

```latex
Future work could evaluate rank-based enrichment methods and replicate
pathway-level patterns in larger independent cohorts. The absence of
significant over-representation in the present analysis suggests that
larger gene sets or alternative enrichment formulations may be required.
```

---

## [ ] P2-15. Osłabić deklarację reprodukowalności w Conclusions

**Zamiast:**

```text
a fully reproducible framework
```

użyć:

```text
a structured and partially reproducible analytical framework
```

lub po uzupełnieniu wersji:

```text
a reproducible downstream analytical framework based on the retained
cohort-level VCF and recorded configurations
```

---

# P3 — bibliografia, skład i język

## [ ] P3-01. Dodać cytowanie g:Profiler

**Lokalizacja:** tabela oprogramowania i Methods → pathway analysis.

---

## [ ] P3-02. Zweryfikować bibliografię 1000 Genomes

**Sprawdzić:** pełny tytuł, czasopismo, rok, wolumin, strony, DOI i autorów konsorcjum.

---

## [ ] P3-03. Ujednolicić format bibliografii

**Sprawdzić globalnie:** kapitalizację DNA, IBD, UC, WES; nazwy krajów; konsorcja; `et al.`; DOI/URL; duplikaty.

---

## [ ] P3-04. Zweryfikować licencję Figure 2.1

**Najbezpieczniej:** stworzyć własny schemat i podpisać:

```text
Author's representation based on the workflow described in [14].
```

---

## [ ] P3-05. Poprawić jakość Figure 1.1

Użyć SVG/PDF albo PNG minimum 300 dpi i sprawdzić czytelność napisów.

---

## [ ] P3-06. Ukryć czerwone ramki hyperlinków

```latex
\usepackage[hidelinks]{hyperref}
```

lub:

```latex
\hypersetup{
    colorlinks=true,
    linkcolor=black,
    citecolor=black,
    urlcolor=blue
}
```

---

## [ ] P3-07. Ujednolicić brytyjski angielski

**Globalne zamiany:**

```text
optimization      → optimisation
visualization     → visualisation
characterization  → characterisation
regularization    → regularisation
prioritization    → prioritisation
summarized        → summarised
analyzed          → analysed
modeling           → modelling
```

Nie zmieniać nazw funkcji, parametrów bibliotek, cytowanych tytułów ani nazw własnych.

---

## [ ] P3-08. Ujednolicić nazewnictwo JPT

**Preferowany termin:**

```text
JPT population-reference samples
```

**Unikać:** `healthy controls`, `JPT controls`, `control cohort`, chyba że jest to techniczna etykieta klasy.

Na osiach rycin rozważyć `JPT reference` zamiast `Control`.

---

## [ ] P3-09. Ujednolicić zapis nazw genów

**W tekście:**

```latex
\textit{NDUFAF2}
\textit{S1PR5}
\textit{ADCY3}
\textit{IL23R}
```

Nie kursywizować nazw białek, szlaków i narzędzi.

---

## [ ] P3-10. Ujednolicić zapis `machine learning`

- `machine learning` jako rzeczownik,
- `machine-learning` jako przymiotnik.

---

## [ ] P3-11. Ujednolicić terminologię cross-validation

**Preferowane:** nested cross-validation, repeated nested cross-validation, outer test fold, inner cross-validation, repetition-level estimate.

---

## [ ] P3-12. Sprawdzić puste strony i rozmieszczenie floatów

Szczególnie pathway analysis, diagram Methods, szerokie tabele i podpisy bez rycin.

---

# Globalne wyszukiwanie i zamiana

Przed ostateczną kompilacją wyszukać:

```text
[?]
[VERSION]
Not specified
planned comparisons
unrelated individuals
pathogenicity scores
healthy controls
harmonization
optimization
visualization
characterization
regularization
prioritization
accuracy were used
LASSO-Regularized
```

Każde wystąpienie zweryfikować kontekstowo.

---

# Końcowa kontrola jakości przed oddaniem

## [ ] Kompilacja

- brak `undefined citations`,
- brak `undefined references`,
- brak brakujących plików graficznych,
- brak placeholderów,
- brak przepełnionych tabel i niewidocznych rycin.

## [ ] Spójność liczb

- 144 osoby,
- 74 UC,
- 70 JPT,
- 488 genów panelu,
- 537 wariantów,
- 545 mapowań wariant–gen,
- 215 predyktorów,
- 10 powtórzeń,
- 5 outer folds,
- 3 inner folds,
- 50 outer test folds,
- 1440 pooled predictions per model.

## [ ] Spójność modeli

- Random Forest,
- Naive Bayes,
- XGBoost,
- Logistic Regression,
- k-NN,
- SVM.

## [ ] Spójność metryk

**Główna:** AUROC.  
**Uzupełniające:** balanced accuracy, sensitivity, specificity, precision, F1-score.  
**Opisowe:** ROC curve, precision–recall curve.

## [ ] Spójność interpretacji

- wynik dotyczy internal cohort discrimination,
- nie jest kliniczną walidacją,
- SHAP nie oznacza przyczynowości,
- pathway analysis nie dała istotnych terminów po FDR,
- callability ogranicza interpretację części genów,
- brak external validation pozostaje głównym ograniczeniem.

## [ ] Finalny eksport PDF

- wszystkie ryciny czytelne przy 100%,
- brak czerwonych ramek,
- poprawne numery stron,
- poprawny spis treści,
- poprawna lista skrótów,
- bibliografia bez braków,
- diagram metodologii widoczny.
