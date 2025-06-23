# Heredity

This project simulates the probability of gene and trait inheritance in a family, based on genetic information and known traits. It uses Bayesian inference to compute the likelihood that each person has a certain number of copies of a gene and whether they exhibit a particular trait.

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
