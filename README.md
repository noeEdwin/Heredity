
#🧬 Heredity – Inference with Probabilistic Models

This project uses a probabilistic model to infer the likelihood that each person in a family has a particular genetic trait (e.g., hearing impairment caused by mutations in the GJB2 gene), based on data about their parents and observed traits.

Using a Bayesian Network, the program:

    Computes the joint probability of various gene/trait combinations in the family.

    Aggregates these probabilities into final probability distributions for each individual.

    Normalizes these distributions to provide meaningful, interpretable results.

The model takes into account:

    Inheritance rules of dominant/recessive genes.

    Gene mutation probabilities.

    The observed presence or absence of a trait.

How It Works

The algorithm performs exhaustive enumeration over all possible gene/trait configurations. For each configuration, it:

    Calculates the joint probability of that scenario.

    Updates a running total of probabilities for each person.

    Normalizes the results so each distribution sums to 1.

This results in a complete probabilistic inference about each person’s gene copies and whether they exhibit the trait.


## How It Works

Given a CSV file describing a family (with columns for name, mother, father, and trait), the program:

1. Loads the family data.
2. Iterates over all possible gene and trait assignments.
3. Calculates the joint probability for each assignment, considering inheritance and mutation.
4. Aggregates and normalizes the probabilities for each person.

## Usage

To run the program, use:

```sh
python heredity.py family0.csv
