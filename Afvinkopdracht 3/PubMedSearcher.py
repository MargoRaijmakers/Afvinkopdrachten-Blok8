# Margo Raijmakers
# 09-05-2021
# Afvinkopdracht 3: Text Mining

from itertools import chain, combinations
from Bio import Entrez, Medline
import random
import pubchempy as pcp


def file_reader(filename):
    """Deze functie leest het opgegeven bestand in en zet de inhoud in
    een lijst.

    :param filename: de naam van het bestand
    :return: een lijst met de inhoud van het bestand
    """
    content = []
    with open(filename, "r") as file:
        for line in file:
            content.append(line.strip())
    return content


def make_concepts_list(compounds, genes, molecular_effects):
    """Deze functie maakt van de compounds, genes en molecular_effects
    lijsten 1 lijst. Dit zijn de begrippen.

    :param compounds: de lijst met compounds
    :param genes: de lijst met genes
    :param molecular_effects: de lijst met molecular effects
    :return: de lijst met begrippen
    """
    return list(chain(compounds, genes, molecular_effects))


def get_synonyms(compounds):
    """Deze functie haalt de synoniemen op.

    :param compounds: een lijst van de compounds
    :return: een dictionary met als key de compound en als value de
    synoniemen
    """
    comp_syn_dict = {}
    for compound in compounds:
        try:
            comp_syn_dict[compound] = pcp.get_synonyms(compound, 'name')[0]\
                .get('Synonym')
        except IndexError:
            comp_syn_dict[compound] = []
    return comp_syn_dict


def get_most_combi(concepts, comp_syn_dict):
    """Deze functie haalt de ids op van de artikelen van de meest
    voorkomende combinaties.

    :param concepts: de lijst met begrippen
    :param comp_syn_dict: een dictionary met als key de compound en als
    value de synoniemen
    :return value: het aantal artikelen van de meest voorkomende combi's
    :return most_combi: een lijst met de meest voorkomende combi's
    :return ids_list: een lijst van ids van de artikelen van de meest
    voorkomende combi's
    """
    value = 0
    most_combi = []
    ids_list = []
    for concept1, concept2 in combinations(concepts, 2):
        # Maak de query
        query = make_query(comp_syn_dict, concept1, concept2)
        # Haal de artikelen op
        handle = Entrez.esearch(db="pubmed", term=query)
        record = Entrez.read(handle)
        handle.close()
        count = int(record["Count"])
        # Sla alleen de ID's op als de count hoger of gelijk is
        if count > value:
            value = count
            most_combi.clear()
            most_combi.append(concept1 + " en " + concept2)
            ids_list.clear()
            ids_list.append(record["IdList"])
        elif count == value:
            most_combi.append(concept1 + " en " + concept2)
            ids_list.append(record["IdList"])
        else:
            pass
    return value, most_combi, ids_list


def make_query(comp_syn_dict, concept1, concept2):
    """Deze functie maakt de query.

    :param comp_syn_dict: een dictionary met als key de compound en als
    value de synoniemen
    :param concept1: het eerste begrip
    :param concept2: het tweede begrip
    :return: de gemaakte query
    """
    query = ""
    # Voor het eerste begrip
    if comp_syn_dict.get(concept1) is not None \
            and len(comp_syn_dict.get(concept1)) != 0:
        query += "((" + concept1 + " [tiab] OR " + \
                 str(comp_syn_dict.get(concept1))\
                     .replace("', '", " [tiab] OR ")\
                     .replace("['", "").replace("']", " [tiab]") + ")"
    else:
        query += "(" + concept1 + " [tiab]"
    # Voor het tweede begrip
    if comp_syn_dict.get(concept2) is not None \
            and len(comp_syn_dict.get(concept2)) != 0:
        query += " AND ((" + concept2 + " [tiab] OR " + \
                 str(comp_syn_dict.get(concept2))\
                     .replace("', '", " [tiab] OR ")\
                     .replace("['", "").replace("']", " [tiab]") + "))"
    else:
        query += " AND " + concept2 + " [tiab])"
    return query


def choose_random(most_combi, ids_list):
    """Deze functie kiest random een van de meest voorkomende combi's.

    :param most_combi: een lijst met de meest voorkomende combi's
    :param ids_list: een lijst van ids van de artikelen van de meest
    voorkomende combi's
    :return chosen_combi: de gekozen combi
    :return chosen_ids: de ids van de artikelen die horen bij de gekozen
    combi
    """
    chosen_combi = random.choice(most_combi)
    chosen_ids = ids_list[most_combi.index(chosen_combi)]
    return chosen_combi, chosen_ids


def get_articles(chosen_ids):
    """Deze functie haalt de artikelen op van de gekozen meest
    voorkomende combi.

    :param chosen_ids: de ids van de artikelen die horen bij de gekozen
    combi
    :return: de artikelen
    """
    handle = Entrez.efetch(db="pubmed", id=chosen_ids, rettype="medline",
                           retmode="text")
    records = Medline.parse(handle)
    records = list(records)
    return records


def main():
    compounds_filename = "compounds.txt"
    genes_filename = "genes.txt"
    molecular_effects_filename = "molecular_effects.txt"
    compounds = file_reader(compounds_filename)
    genes = file_reader(genes_filename)
    molecular_effects = file_reader(molecular_effects_filename)
    # [0:4] zodat de dataset niet te groot is
    concepts = make_concepts_list(compounds[0:4], genes[0:4],
                                  molecular_effects[0:4])
    comp_syn_dict = get_synonyms(compounds[0:4])

    compounds = None
    molecular_effects = None
    genes = None

    Entrez.email = "margoraijmakers@gmail.com"
    value, most_combi, ids_list = get_most_combi(concepts, comp_syn_dict)
    chosen_combi, chosen_ids = choose_random(most_combi, ids_list)

    print("Meest voorkomende combinatie:", chosen_combi.replace("\"", "")
          .replace("AND", "en").replace(" [tiab]", ""))
    print("Aantal:", value)

    records = get_articles(chosen_ids)
    for record in records:
        print("Abstract:", record.get("AB"))
        print("Auteurs:", str(record.get("FAU")).replace("'", "")
              .replace("[", "").replace("]", ""))
        print("")


main()
