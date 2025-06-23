import csv
import itertools
import sys
import random

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    parents = []
    childs = []
    parents_probabilities = {}
    childs_probabilities = {}
    joint = 1

    for person in people:
        if check_parents(people, person):
            childs.append(person)
        else:
            parents.append(person)

    #giving probabilities to the parents
    for parent in parents:
        parents_probabilities[parent] = {}
        if parent in one_gene:
            parents_probabilities[parent]["gene"] = PROBS["gene"][1]
        elif parent in two_genes:
            parents_probabilities[parent]["gene"] = PROBS["gene"][2]
        else:
            parents_probabilities[parent]["gene"] = PROBS["gene"][0]

        parents_probabilities[parent]["trait"] = PROBS["trait"][check_for_trait(parents_probabilities[parent]["gene"])][parent in have_trait]

    for child in childs:
        childs_probabilities[child] = {}
        father = people[child]["father"]
        mother = people[child]["mother"]

        probability_father = passing_prob(father, one_gene, two_genes)
        probability_mother = passing_prob(mother, one_gene, two_genes)

        if child in one_gene:
            gene_prob = probability_mother * (1 - probability_father) + probability_father * (1-probability_mother)
            childs_probabilities[child]["gene"] = gene_prob
            childs_probabilities[child]["trait"] = PROBS["trait"][1][child in have_trait]

        elif child in two_genes:
            gene_prob = probability_father * probability_mother
            childs_probabilities[child]["gene"] = gene_prob
            childs_probabilities[child]["trait"] = PROBS["trait"][2][child in have_trait]

        else:
            gene_prob = (1 - probability_father) * (1-probability_mother)
            childs_probabilities[child]["gene"] = gene_prob
            childs_probabilities[child]["trait"] = PROBS["trait"][0][child in have_trait]

    for parent in list(parents_probabilities.keys()):
        joint *= parents_probabilities[parent]["gene"] * parents_probabilities[parent]["trait"]

    for child in list(childs_probabilities.keys()):
        joint *= childs_probabilities[child]["gene"] * childs_probabilities[child]["trait"]

    return joint

def passing_prob(person, one_gene, two_genes):
    if person in one_gene:
        return 0.5
    elif person in two_genes:
        return 1 - PROBS["mutation"]
    else:
        return PROBS["mutation"]

def check_for_trait(gene_probability):
    if gene_probability == 0.01:
        return 2
    elif gene_probability == 0.03:
        return 1
    else:
        return 0

def check_parents(people, person):
    return people[person]["father"] != None and people[person]["mother"] != None

def return_probability(field):
    return random.choices(list(PROBS[field].values()), list(PROBS[field].values()), k=1)[0]

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        # GENES
        if person in one_gene:
            probabilities[person]["gene"][1] += p
        elif person in two_genes:
            probabilities[person]["gene"][2] += p
        else:
            probabilities[person]["gene"][0] += p
        probabilities[person]["trait"][person in have_trait] += p

def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in list(probabilities.keys()):
        one_gene = probabilities[person]["gene"][1]
        two_gene = probabilities[person]["gene"][2]
        zero_gene = probabilities[person]["gene"][0]
        #normalizing distributions
        if zero_gene + two_gene + one_gene != 1:
            alpha = 1/(zero_gene + two_gene + one_gene)
            probabilities[person]["gene"][1] *= alpha
            probabilities[person]["gene"][2] *= alpha
            probabilities[person]["gene"][0] *= alpha

        have_trait_person = probabilities[person]["trait"][True]
        have_not_trait_person = probabilities[person]["trait"][False]

        if have_not_trait_person + have_trait_person != 1:
            alpha = 1/(have_trait_person+have_not_trait_person)
            probabilities[person]["trait"][True] *= alpha
            probabilities[person]["trait"][False] *= alpha

if __name__ == "__main__":
    main()
